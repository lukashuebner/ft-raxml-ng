#include <algorithm>

#include "TreeInfo.hpp"
#include "ParallelContext.hpp"
#include "Checkpoint.hpp"
#include "io/file_io.hpp"
#include <cassert>
#include <signal.h>

using namespace std;

//extern "C" void (*pllmod_treeinfo_update_recovery_tree_benchmarked)(pllmod_treeinfo_t * treeinfo);

TreeInfo::TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
                    const IDVector& tip_msa_idmap,
                    const PartitionAssignment& part_assign,
                    std::function<PartitionAssignment(bool)> calculate_assignment_cb) :
                    redo_assignment_cb(calculate_assignment_cb)
{
  init_profiler();
  init(opts, tree, parted_msa, tip_msa_idmap, part_assign, std::vector<uintVector>());
  mini_checkpoint(true, true);
}

TreeInfo::TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
                    const IDVector& tip_msa_idmap,
                    const PartitionAssignment& part_assign,
                    const std::vector<uintVector>& site_weights)
{
  init_profiler();
  init(opts, tree, parted_msa, tip_msa_idmap, part_assign, site_weights);
  mini_checkpoint(true, true);
}

shared_ptr<ProfilerRegister> TreeInfo::_profiler_register = nullptr;
void TreeInfo::init_profiler() {
  _profiler_register = ProfilerRegister::getInstance();
  assert(_profiler_register != nullptr);
}

void TreeInfo::init(const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
                    const IDVector& tip_msa_idmap,
                    const PartitionAssignment& part_assign,
                    const std::vector<uintVector>& site_weights)
{
  _brlen_min = opts.brlen_min;
  _brlen_max = opts.brlen_max;
  _brlen_opt_method = opts.brlen_opt_method;
  _check_lh_impr = opts.safety_checks.isset(SafetyCheck::model_lh_impr);
  _partition_contributions.resize(parted_msa.part_count());

  // TODO Has this to be updated when the site assignment changes? --> Assume not
  _pll_treeinfo = pllmod_treeinfo_create(pll_utree_graph_clone(&tree.pll_utree_root()),
                                         tree.num_tips(),
                                         parted_msa.part_count(), opts.brlen_linkage);
  libpll_check_error("ERROR creating treeinfo structure");
  assert(_pll_treeinfo);

  if (ParallelContext::threads_per_group() > 1)
  {
    pllmod_treeinfo_set_parallel_context(_pll_treeinfo, (void *) nullptr,
                                         ParallelContext::parallel_reduce_cb);
  }

  init_partitions(opts, tree, parted_msa, tip_msa_idmap, part_assign, site_weights);
  partition_reinit_info = make_shared<partition_reinit_info_t>(opts, parted_msa, tip_msa_idmap, site_weights);

  // Enable profiling for tree updates
  pllmod_treeinfo_update_recovery_tree_set_benchmarked([](pllmod_treeinfo_t * treeinfo) {
    _profiler_register->profileFunction([&treeinfo]() {
      pllmod_treeinfo_update_recovery_tree(treeinfo);
    }, "UpdateTree");
  });
}

void TreeInfo::reinit_partitions(const PartitionAssignment& part_assign) {
  Tree tree;
  assert(_profiler_register != nullptr);
  _profiler_register->profileFunction([&]() {
    if (_pll_treeinfo->recovery_tree) {
      tree = Tree(*(_pll_treeinfo->recovery_tree));
    } else {
      tree = Tree(*(_pll_treeinfo->tree));
    }
  }, "CreateTree");

  // TODO: Maybe we can reuse some parts of the treeinfo structure? For this, we need to update
  // the tree in the treeinfo structure.
  // If you try: Check for bug in case of node failure before alpha optimization

  _profiler_register->profileFunction([&]() {
    for (unsigned int i = 0; i < _pll_treeinfo->partition_count; ++i)
    {
      if (_pll_treeinfo->partitions[i])
        pll_partition_destroy(_pll_treeinfo->partitions[i]);
    }
  }, "DestroyPartitions");

  _profiler_register->profileFunction([&]() {
    pllmod_treeinfo_reset_partitions(_pll_treeinfo);
  }, "TreeinfoResetPartitions");

  _profiler_register->profileFunction([&]() {
    pll_utree_graph_destroy(_pll_treeinfo->root, NULL);
  }, "UTreeGraphDestroy");

  _profiler_register->profileFunction([&]() {
    pllmod_treeinfo_destroy(_pll_treeinfo);
  }, "TreeinfoDestroy");

  //pllmod_treeinfo_set_tree(_pll_treeinfo, _pll_treeinfo->tree);
  _profiler_register->profileFunction([&]() {
    init(partition_reinit_info->opts,
                  tree,
                  partition_reinit_info->parted_msa,
                  partition_reinit_info->tip_msa_idmap,
                  part_assign,
                  partition_reinit_info->site_weights);
  }, "TreeinfoInit");

  _profiler_register->profileFunction([&]() {
    pllmod_treeinfo_update_partials_and_clvs(_pll_treeinfo);
  }, "TreeinfoUpdatePartialsAndCLVs");
  
  // pllmod_treeinfo_invalidate_all(_pll_treeinfo);
  // pllmod_treeinfo_compute_loglh(_pll_treeinfo, 0);

  LOG_DEBUG << "Restored the following tree:" << endl;
  LOG_DEBUG << to_newick_string_rooted(tree, 0) << endl;
}

void TreeInfo::init_partitions(const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
                               const IDVector& tip_msa_idmap,
                               const PartitionAssignment& part_assign,
                               const std::vector<uintVector>& site_weights)
{
  _parts_master.clear();
  int optimize_branches = opts.optimize_brlen ? PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE : 0;
  double total_weight = 0;

  for (size_t p = 0; p < parted_msa.part_count(); ++p)
  {
    const PartitionInfo& pinfo = parted_msa.part_info(p);
    const auto& weights = site_weights.empty() ? pinfo.msa().weights() : site_weights.at(p);
    int params_to_optimize = opts.optimize_model ? pinfo.model().params_to_optimize() : 0;
    params_to_optimize |= optimize_branches;

    _partition_contributions[p] = std::accumulate(weights.begin(), weights.end(), 0);
    total_weight += _partition_contributions[p];

    PartitionAssignment::const_iterator part_range = part_assign.find(p);
    if (part_range != part_assign.end())
    {
      /* create and init PLL partition structure */
      pll_partition_t * partition = create_pll_partition(opts, pinfo, tip_msa_idmap,
                                                         *part_range, weights);

      int retval = pllmod_treeinfo_init_partition(_pll_treeinfo, p, partition,
                                                  params_to_optimize,
                                                  pinfo.model().gamma_mode(),
                                                  pinfo.model().alpha(),
                                                  pinfo.model().ratecat_submodels().data(),
                                                  pinfo.model().submodel(0).rate_sym().data());

      if (!retval)
      {
        assert(pll_errno);
        libpll_check_error("ERROR adding treeinfo partition");
      }

      // set per-partition branch lengths or scalers
      if (opts.brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
      {
        assert (_pll_treeinfo->brlen_scalers);
        _pll_treeinfo->brlen_scalers[p] = pinfo.model().brlen_scaler();
      }
      else if (opts.brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED && !tree.partition_brlens().empty())
      {
        assert(_pll_treeinfo->branch_lengths[p]);
        memcpy(_pll_treeinfo->branch_lengths[p], tree.partition_brlens(p).data(),
               tree.num_branches() * sizeof(double));
      }

      if (part_range->master()) {
        _parts_master.insert(p);
      }
    }
    else
    {
      // this partition will be processed by other threads, but we still need to know
      // which parameters to optimize
      _pll_treeinfo->params_to_optimize[p] = params_to_optimize;
    }
  }

  // finalize partition contribution computation
  for (auto& c: _partition_contributions)
    c /= total_weight;
}

TreeInfo::~TreeInfo ()
{
  if (_pll_treeinfo)
  {
    for (unsigned int i = 0; i < _pll_treeinfo->partition_count; ++i)
    {
      if (_pll_treeinfo->partitions[i])
        pll_partition_destroy(_pll_treeinfo->partitions[i]);
    }

    pll_utree_graph_destroy(_pll_treeinfo->root, NULL);
    pllmod_treeinfo_destroy(_pll_treeinfo);
  }
}

void TreeInfo::assert_lh_eq(double lh1, double lh2) {
  assert(abs(lh1 - lh2) < -lh1 * RAXML_LOGLH_TOLERANCE);
}

void TreeInfo::assert_lh_improvement(double old_lh, double new_lh, const std::string& where)
{
  if (_check_lh_impr && !(old_lh - new_lh < -new_lh * RAXML_LOGLH_TOLERANCE))
  {
    throw runtime_error((where.empty() ? "" : "[" + where + "] ") +
                        "Worse log-likelihood after optimization!\n" +
                        "Old: " + to_string(old_lh) + "\n"
                        "New: " + to_string(new_lh) + "\n" +
                        "NOTE: You can disable this check with '--force model_lh_impr'");
  }
}


Tree TreeInfo::tree() const
{
  if (!_pll_treeinfo)
    return Tree();

  Tree tree(*_pll_treeinfo->tree);

  if (_pll_treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    // set per-partition branch lengths
    for (unsigned int i = 0; i < _pll_treeinfo->partition_count; ++i)
    {
      assert(_pll_treeinfo->branch_lengths[i]);
      doubleVector brlens(_pll_treeinfo->branch_lengths[i],
                          _pll_treeinfo->branch_lengths[i] + tree.num_branches());
      tree.add_partition_brlens(std::move(brlens));
    }

    // compute a weighted average of per-partition brlens
    tree.apply_avg_brlens(_partition_contributions);
  }

  return tree;
}

Tree TreeInfo::tree(size_t partition_id) const
{
  if (!_pll_treeinfo)
    return Tree();

  if (partition_id >= _pll_treeinfo->partition_count)
    throw out_of_range("Partition ID out of range");

  PllUTreeUniquePtr pll_utree(pllmod_treeinfo_get_partition_tree(_pll_treeinfo, partition_id));

  if (!pll_utree)
  {
    assert(pll_errno);
    libpll_check_error("treeinfo: cannot get partition tree");
  }

  return Tree(pll_utree);
}

void TreeInfo::tree(const Tree& tree)
{
  _pll_treeinfo->root = pll_utree_graph_clone(&tree.pll_utree_root());
}

double TreeInfo::loglh(bool incremental)
{
  double loglh = fault_tolerant_call("loglh computation", [this, incremental]() -> double {
      return pllmod_treeinfo_compute_loglh(_pll_treeinfo, incremental ? 1 : 0);
  }, false, false);

  #ifdef NON_FAILURE_TOLERANT_ASSERTS
  // Check ,if the loglh is equal across all PEs
  double my_loglh = loglh;
  ParallelContext::parallel_reduce(&loglh, 1, PLLMOD_COMMON_REDUCE_MAX);
  assert(my_loglh == loglh); // Should be exactly equal, without tolerance
  #endif

  return loglh;
}

void TreeInfo::model(size_t partition_id, const Model& model)
{
  if (partition_id >= _pll_treeinfo->partition_count)
    throw out_of_range("Partition ID out of range");

  if (!_pll_treeinfo->partitions[partition_id])
    return;

  assign(_pll_treeinfo->partitions[partition_id], model);
  _pll_treeinfo->alphas[partition_id] = model.alpha();
  if (_pll_treeinfo->brlen_scalers)
    _pll_treeinfo->brlen_scalers[partition_id] = model.brlen_scaler();
}

//#define DBG printf

double TreeInfo::optimize_branches(double lh_epsilon, double brlen_smooth_factor)
{
  /* update all CLVs and p-matrices before calling BLO */
  double new_loglh = loglh();

  if (_pll_treeinfo->params_to_optimize[0] & PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE)
  {
    int max_iters = brlen_smooth_factor * RAXML_BRLEN_SMOOTHINGS;
    new_loglh = fault_tolerant_call("branch length opt.",
      [this, lh_epsilon, brlen_smooth_factor, max_iters]() -> double {
        return -1 * pllmod_algo_opt_brlen_treeinfo(_pll_treeinfo,
                                                   _brlen_min,
                                                   _brlen_max,
                                                   lh_epsilon,
                                                   max_iters,
                                                   _brlen_opt_method,
                                                   PLLMOD_OPT_BRLEN_OPTIMIZE_ALL
                                                  );
    }, false, true);

    LOG_DEBUG << "\t - after brlen: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in branch length optimization");
    assert(isfinite(new_loglh));
  }

  /* optimize brlen scalers, if needed */
  if (_pll_treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED &&
      _pll_treeinfo->partition_count > 1)
  {
    new_loglh = fault_tolerant_call("branch length scalers opt.", [this]() -> double {
      return -1 * pllmod_algo_opt_brlen_scalers_treeinfo(_pll_treeinfo,
                                                         RAXML_BRLEN_SCALER_MIN,
                                                         RAXML_BRLEN_SCALER_MAX,
                                                         _brlen_min,
                                                         _brlen_max,
                                                         RAXML_PARAM_EPSILON);
    }, true, false); // TODO: Does this depend on the scaling mode?

    LOG_DEBUG << "\t - after brlen scalers: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in brlen scaler optimization");
    assert(isfinite(new_loglh));
  }

  return new_loglh;
}

void TreeInfo::update_to_new_assignment(bool rebalance) {
  assert(_profiler_register != nullptr);

  PartitionAssignment part_assign;
  _profiler_register->profileFunction([&]() {
    part_assign = redo_assignment_cb(rebalance);
  }, "RedoPartitionAssignment");

  reinit_partitions(part_assign);

  _profiler_register->profileFunction([&]() {  
    for (auto& m: CheckpointManager::all_models())
    {
      LOG_DEBUG << "Restoring model " << m.first << " to:" << endl;
      LOG_DEBUG << m.second << endl;

      // Locally restore the model, if we have a partition with this model assigned
      if (_pll_treeinfo->partitions[m.first]) {
        model(m.first, m.second);
      }
    }
  }, "RestoreModels");
}

//void TreeInfo::may_rebalance(bool force) {
//  if (ParallelContext::ranks_per_group() == 1) {
//    return;
//  }
//
//  double worked_for = _profiler_register->worked_for_ms();
//  ParallelContext::parallel_reduce(&worked_for, sizeof(worked_for), PLLMOD_COMMON_REDUCE_MAX);
//  if (force || worked_for > 10000) {
//    update_to_new_assignment(true);
//    _profiler_register->saveWorkByRank(true);
//  }
//}

void TreeInfo::print_tree() {
  LOG_DEBUG << to_newick_string_rooted(Tree(*(_pll_treeinfo->tree)), 0) << endl;
}

void TreeInfo::mini_checkpoint(bool save_models, bool save_tree) {
  if (save_models) {
    _profiler_register->profileFunction([this]() {
      CheckpointManager::update_models(*this);
    }, "UpdateModels");
  }
  // This operation is local and therefore cannot fail
  if (save_tree) {
    pllmod_treeinfo_update_recovery_tree_benchmarked(_pll_treeinfo);
  }
}

size_t failureCount = 0;
template<class F>
double TreeInfo::fault_tolerant_call(string parameter, F optimizer, bool changes_model, bool changes_tree) {
  //! Do not use loglh() inside this function, as loglh() will compute the log likelihood using
  //! fault_tolerant_call. Use pllmod_treeinfo_compute_loglh() directly instead (might fail)
  double new_loglh = 1;

  #ifdef NON_FAILURE_TOLERANT_ASSERTS
  //_profiler_register->startWorkTimer();
  double pre_fail_loglh = pllmod_treeinfo_compute_loglh(_pll_treeinfo, 0);
  //_profiler_register->endWorkTimer();
  #endif

  for (;;) {
    #ifndef NDEBUG
    string beforeFailureTree, beforeFailureModels;
    beforeFailureTree = to_newick_string_rooted(Tree(*(_pll_treeinfo->tree)), 0);
    beforeFailureModels = CheckpointManager::all_models_to_string();
    assert(beforeFailureTree != "" && beforeFailureModels != "");
    #endif

    try{
      assert(_profiler_register != nullptr);
      //_profiler_register->startWorkTimer();

      // Run optimization code during which rank failure might occur
      new_loglh = optimizer(); 
      // Check whether a rank failure occurred but not all ranks got notified yet.
      ParallelContext::check_for_rank_failure();
      
      //_profiler_register->endWorkTimer();

      mini_checkpoint(changes_model, changes_tree);
      break;
    } catch (ParallelContext::RankFailureException& e) {
      //_profiler_register->discardWorkTimer();
      LOG_ERROR << "Rank failure during " << parameter << endl;
      update_to_new_assignment();
      LOG_PROGR << "Restoration completed successfully." << endl;

      #ifndef NDEBUG
        assert(beforeFailureModels == CheckpointManager::all_models_to_string());
        if (!changes_tree) { // Tree operations are synchronized across all ranks, no need to restore
          assert(beforeFailureTree == to_newick_string_rooted(Tree(*(_pll_treeinfo->tree)), 0));
        }
        #ifdef NON_FAILURE_TOLERANT_ASSERTS
          //_profiler_register->startWorkTimer();
          double loglh = pllmod_treeinfo_compute_loglh(_pll_treeinfo, 0);
          //_profiler_register->endWorkTimer();
          if (parameter == "spr round" || parameter == "branch length opt.") {
            assert_lh_improvement(pre_fail_loglh, loglh);
          } else if (parameter != "branch length scalers") { // TODO: Check if this is expected
            assert_lh_eq(pre_fail_loglh, loglh);
          }
        #endif
      #endif
    }
  }

  #ifdef NON_FAILURE_TOLERANT_ASSERTS
  //_profiler_register->startWorkTimer();
  assert_lh_eq(new_loglh, pllmod_treeinfo_compute_loglh(_pll_treeinfo, 0));
  //_profiler_register->endWorkTimer();
  #endif
  assert(new_loglh <= 0);

  return new_loglh;
}


double TreeInfo::optimize_params(int params_to_optimize, double lh_epsilon)
{
  assert(!pll_errno);

  // TODO: Catch failures more elegantly
  mini_checkpoint(true, true);
  //may_rebalance();

  double
    cur_loglh = loglh(),
    new_loglh = cur_loglh;

  //ParallelContext::set_failure_prob(0.0005);

  /* optimize SUBSTITUTION RATES */
  if (params_to_optimize & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    new_loglh = fault_tolerant_call("substitution rates opt.", [this]() -> double {
      return -1 * pllmod_algo_opt_subst_rates_treeinfo(this->_pll_treeinfo,
                                                       0,
                                                       PLLMOD_OPT_MIN_SUBST_RATE,
                                                       PLLMOD_OPT_MAX_SUBST_RATE,
                                                       RAXML_BFGS_FACTOR,
                                                       RAXML_PARAM_EPSILON);
    }, true, false);

    LOG_DEBUG << "\t - after rates: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in substitution rates optimization");
    assert_lh_improvement(cur_loglh, new_loglh, "RATES");
    cur_loglh = new_loglh;
  }

  /* optimize BASE FREQS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREQUENCIES)
  {
    new_loglh = fault_tolerant_call("frequencies opt.", [this]() -> double {
      return -1 * pllmod_algo_opt_frequencies_treeinfo(this->_pll_treeinfo,
                                                       0,
                                                       PLLMOD_OPT_MIN_FREQ,
                                                       PLLMOD_OPT_MAX_FREQ,
                                                       RAXML_BFGS_FACTOR,
                                                       RAXML_PARAM_EPSILON);
    }, true, false);

    LOG_DEBUG << "\t - after freqs: logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in base frequencies optimization");
    assert_lh_improvement(cur_loglh, new_loglh, "FREQS");
    cur_loglh = new_loglh;
  }

  // TODO: co-optimization of PINV and ALPHA, mb with multiple starting points
  if (0 &&
      (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA) &&
      (params_to_optimize & PLLMOD_OPT_PARAM_PINV))
  {
    new_loglh = fault_tolerant_call("alpha & pinv opt.", [this]() -> double {
      return -1 * pllmod_algo_opt_alpha_pinv_treeinfo(_pll_treeinfo,
                                                      0,
                                                      PLLMOD_OPT_MIN_ALPHA,
                                                      PLLMOD_OPT_MAX_ALPHA,
                                                      PLLMOD_OPT_MIN_PINV,
                                                      PLLMOD_OPT_MAX_PINV,
                                                      RAXML_BFGS_FACTOR,
                                                      RAXML_PARAM_EPSILON);
    }, true, false);

    LOG_DEBUG << "\t - after a+i  : logLH = " << new_loglh << endl;

    libpll_check_error("ERROR in alpha/p-inv parameter optimization");
    assert_lh_improvement(cur_loglh, new_loglh, "ALPHA+PINV");
    cur_loglh = new_loglh;
  }
  else
  {
    /* optimize ALPHA */
    if (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA)
    {
      new_loglh = fault_tolerant_call("alpha opt.", [this]() -> double {
        return -1 * pllmod_algo_opt_onedim_treeinfo(this->_pll_treeinfo,
                                                    PLLMOD_OPT_PARAM_ALPHA,
                                                    PLLMOD_OPT_MIN_ALPHA,
                                                    PLLMOD_OPT_MAX_ALPHA,
                                                    RAXML_PARAM_EPSILON);
      }, true, false);

      LOG_DEBUG << "\t - after alpha: logLH = " << new_loglh << endl;

      libpll_check_error("ERROR in alpha parameter optimization");
      assert_lh_improvement(cur_loglh, new_loglh, "ALPHA");
      cur_loglh = new_loglh;
    }

    /* optimize PINV */
    if (params_to_optimize & PLLMOD_OPT_PARAM_PINV)
    {
      new_loglh = fault_tolerant_call("p-inv opt.", [this]() -> double {
        return -1 * pllmod_algo_opt_onedim_treeinfo(this->_pll_treeinfo,
                                                    PLLMOD_OPT_PARAM_PINV,
                                                    PLLMOD_OPT_MIN_PINV,
                                                    PLLMOD_OPT_MAX_PINV,
                                                    RAXML_PARAM_EPSILON);
      }, true, false);

      LOG_DEBUG << "\t - after p-inv: logLH = " << new_loglh << endl;

      libpll_check_error("ERROR in p-inv optimization");
      assert_lh_improvement(cur_loglh, new_loglh, "PINV");
      cur_loglh = new_loglh;
    }
  }

  /* optimize FREE RATES and WEIGHTS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREE_RATES)
  {
    cur_loglh = loglh();

    new_loglh = fault_tolerant_call("rates & weights opt.", [this]() -> double {
      return -1 * pllmod_algo_opt_rates_weights_treeinfo (this->_pll_treeinfo,
                                                          RAXML_FREERATE_MIN,
                                                          RAXML_FREERATE_MAX,
                                                          RAXML_BFGS_FACTOR,
                                                          RAXML_PARAM_EPSILON);
    }, true, false);
    
    /* normalize scalers and scale the branches accordingly */
    if (_pll_treeinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED &&
        _pll_treeinfo->partition_count > 1) {
      fault_tolerant_call("brlen scaler normalisation", [this]() -> double {
        pllmod_treeinfo_normalize_brlen_scalers(_pll_treeinfo);
        return pllmod_treeinfo_compute_loglh(_pll_treeinfo, 0);
      }, true, false); // TODO: Does this depend on the scaling mode?
    }

    LOG_DEBUG << "\t - after freeR: logLH = " << new_loglh << endl;
//    LOG_DEBUG << "\t - after freeR/crosscheck: logLH = " << loglh() << endl;

    libpll_check_error("ERROR in FreeRate rates/weights optimization");
    assert_lh_improvement(cur_loglh, new_loglh, "FREE RATES");
    cur_loglh = new_loglh;
  }

  if (params_to_optimize & PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE)
  {
    // Branch optimization has it's own failure mitigation
    new_loglh = optimize_branches(lh_epsilon, 0.25);

    assert_lh_improvement(cur_loglh, new_loglh, "BRLEN");
    cur_loglh = new_loglh;
  }

  mini_checkpoint(true, true);

  return new_loglh;
}

double TreeInfo::spr_round(spr_round_params& params)
{
  //may_rebalance();

  double loglh = fault_tolerant_call("spr round", [this, &params]() -> double {
    return pllmod_algo_spr_round(_pll_treeinfo, params.radius_min, params.radius_max,
                               params.ntopol_keep, params.thorough, _brlen_opt_method,
                               _brlen_min, _brlen_max, RAXML_BRLEN_SMOOTHINGS,
                               0.1,
                               params.subtree_cutoff > 0. ? &params.cutoff_info : nullptr,
                               params.subtree_cutoff);
  }, false, true);

  libpll_check_error("ERROR in SPR round");

  assert(isfinite(loglh) && loglh);

  return loglh;
}

void TreeInfo::set_topology_constraint(const Tree& cons_tree)
{
  if (!cons_tree.empty())
  {
    int retval = pllmod_treeinfo_set_constraint_tree(_pll_treeinfo, &cons_tree.pll_utree());
    if (!retval)
      libpll_check_error("ERROR: Cannot set topological constraint");
  }
}

// No MPI-parallelization support -> no fault tolerance
void TreeInfo::compute_ancestral(const AncestralStatesSharedPtr& ancestral,
                                 const PartitionAssignment& part_assign)
{
  pllmod_ancestral_t * pll_ancestral = pllmod_treeinfo_compute_ancestral(_pll_treeinfo);

  if (!pll_ancestral)
    libpll_check_error("Unable to compute ancestral states", true);

  assert(pll_ancestral->partition_count > 0 && pll_ancestral->partition_indices);

  if (ParallelContext::master_thread())
    assign_tree(*ancestral, *pll_ancestral);

  assign_probs(*ancestral, *pll_ancestral, part_assign);

  pllmod_treeinfo_destroy_ancestral(pll_ancestral);
}

void assign(PartitionedMSA& parted_msa, const TreeInfo& treeinfo)
{
  const pllmod_treeinfo_t& pll_treeinfo = treeinfo.pll_treeinfo();

  if (parted_msa.part_count() != pll_treeinfo.partition_count)
    throw runtime_error("Incompatible arguments");

  for (size_t p = 0; p < parted_msa.part_count(); ++p)
  {
    if (!pll_treeinfo.partitions[p])
      continue;

    Model model(parted_msa.model(p));
    assign(model, treeinfo, p);
    parted_msa.model(p, move(model));
  }
}

void assign(Model& model, const TreeInfo& treeinfo, size_t partition_id)
{
  const pllmod_treeinfo_t& pll_treeinfo = treeinfo.pll_treeinfo();

  if (partition_id >= pll_treeinfo.partition_count)
    throw out_of_range("Partition ID out of range");

  if (!pll_treeinfo.partitions[partition_id])
    return;

  assign(model, pll_treeinfo.partitions[partition_id]);
  model.alpha(pll_treeinfo.alphas[partition_id]);
  if (pll_treeinfo.brlen_scalers)
    model.brlen_scaler(pll_treeinfo.brlen_scalers[partition_id]);
}

void build_clv(ProbVector::const_iterator probs, size_t sites, WeightVector::const_iterator weights,
               pll_partition_t* partition, bool normalize, std::vector<double>& clv)
{
  const auto states = partition->states;
  auto clvp = clv.begin();

  for (size_t i = 0; i < sites; ++i)
  {
    if (weights[i] > 0)
    {
      double sum = 0.;
      for (size_t j = 0; j < states; ++j)
        sum += probs[j];

      for (size_t j = 0; j < states; ++j)
      {
        if (sum > 0.)
          clvp[j] =  normalize ? probs[j] / sum : probs[j];
        else
          clvp[j] = 1.0;
      }

      clvp += states;
    }

    /* NB: clv has to be padded, but msa arrays are not! */
    probs += states;
  }

  assert(clvp == clv.end());
}

void set_partition_tips(const Options& opts, const MSA& msa, const IDVector& tip_msa_idmap,
                        const PartitionRange& part_region,
                        pll_partition_t* partition, const pll_state_t * charmap)
{
  /* get "true" sequence offset considering that MSA can be partially loaded */
  auto seq_offset = msa.get_local_offset(part_region.start);

//  printf("\n\n rank %lu, GLOBAL OFFSET %lu, LOCAL OFFSET %lu \n\n", ParallelContext::proc_id(), part_region.start, seq_offset);

  /* set pattern weights */
  if (!msa.weights().empty())
    pll_set_pattern_weights(partition, msa.weights().data() + seq_offset);

  if (opts.use_prob_msa && msa.probabilistic())
  {
    assert(!(partition->attributes & PLL_ATTRIB_PATTERN_TIP));
    assert(partition->states == msa.states());

    auto normalize = !msa.normalized();
    auto weights_start = msa.weights().cbegin() + seq_offset;

    // we need a libpll function for that!
    auto clv_size = partition->sites * partition->states;
    std::vector<double> tmp_clv(clv_size);
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      auto prob_start = msa.probs(seq_id, seq_offset);
      build_clv(prob_start, partition->sites, weights_start, partition, normalize, tmp_clv);
      pll_set_tip_clv(partition, tip_id, tmp_clv.data(), PLL_FALSE);
    }
  }
  else
  {
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      pll_set_tip_states(partition, tip_id, charmap, msa.at(seq_id).c_str() + seq_offset);
    }
  }
}

void set_partition_tips(const Options& opts, const MSA& msa, const IDVector& tip_msa_idmap,
                        const PartitionRange& part_region,
                        pll_partition_t* partition, const pll_state_t * charmap,
                        const WeightVector& weights)
{
  assert(!weights.empty());

  const auto pstart = msa.get_local_offset(part_region.start);
  const auto plen = part_region.length;
  const auto pend = pstart + plen;

  /* compress weights array by removing all zero entries */
  uintVector comp_weights;
  for (size_t j = pstart; j < pend; ++j)
  {
    if (weights[j] > 0)
      comp_weights.push_back(weights[j]);
  }

  /* now set tip sequences, ignoring all columns with zero weights */
  if (opts.use_prob_msa && msa.probabilistic())
  {
    assert(!(partition->attributes & PLL_ATTRIB_PATTERN_TIP));
    assert(partition->states == msa.states());

    auto normalize = !msa.normalized();
    auto weights_start = msa.weights().cbegin() + pstart;

    // we need a libpll function for that!
    auto clv_size = plen * partition->states;
    std::vector<double> tmp_clv(clv_size);
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      auto prob_start = msa.probs(seq_id, pstart);
      build_clv(prob_start, plen, weights_start, partition, normalize, tmp_clv);
      pll_set_tip_clv(partition, tip_id, tmp_clv.data(), PLL_FALSE);
    }
  }
  else
  {
    std::vector<char> bs_seq(plen);
    for (size_t tip_id = 0; tip_id < partition->tips; ++tip_id)
    {
      auto seq_id = tip_msa_idmap.empty() ? tip_id : tip_msa_idmap[tip_id];
      const char * full_seq = msa.at(seq_id).c_str();
      size_t pos = 0;
      for (size_t j = pstart; j < pend; ++j)
      {
        if (weights[j] > 0)
          bs_seq[pos++] = full_seq[j];
      }
      assert(pos == comp_weights.size());

      pll_set_tip_states(partition, tip_id, charmap, bs_seq.data());
    }
  }

  pll_set_pattern_weights(partition, comp_weights.data());
}

pll_partition_t* create_pll_partition(const Options& opts, const PartitionInfo& pinfo,
                                      const IDVector& tip_msa_idmap,
                                      const PartitionRange& part_region, const uintVector& weights)
{
  const MSA& msa = pinfo.msa();
  const Model& model = pinfo.model();
  const auto pstart = msa.get_local_offset(part_region.start);

//  printf("\n\n rank %lu, GLOBAL OFFSET %lu, LOCAL OFFSET %lu \n\n", ParallelContext::proc_id(), part_region.start, pstart);

  /* part_length doesn't include columns with zero weight */
  const size_t part_length = weights.empty() ? part_region.length :
                             std::count_if(weights.begin() + pstart,
                                           weights.begin() + pstart + part_region.length,
                                           [](uintVector::value_type w) -> bool
                                             { return w > 0; }
                                           );

  unsigned int attrs = opts.simd_arch;

  if (opts.use_rate_scalers && model.num_ratecats() > 1)
  {
    attrs |= PLL_ATTRIB_RATE_SCALERS;
  }

  if (opts.use_repeats)
  {
    assert(!(opts.use_prob_msa));
    attrs |= PLL_ATTRIB_SITE_REPEATS;
  }
  else if (opts.use_tip_inner)
  {
    assert(!(opts.use_prob_msa));
    // 1) SSE3 tip-inner kernels are not implemented so far, so generic version will be faster
    // 2) same for state-rich models
    if (opts.simd_arch != PLL_ATTRIB_ARCH_SSE && model.num_states() <= 20)
    {
      // TODO: use proper auto-tuning
      const unsigned long min_len_ti = model.num_states() > 4 ? 40 : 100;
      if ((unsigned long) part_length > min_len_ti)
        attrs |= PLL_ATTRIB_PATTERN_TIP;
    }
  }

  // NOTE: if partition is split among multiple threads, asc. bias correction must be applied only once!
  if (model.ascbias_type() == AscBiasCorrection::lewis ||
      (model.ascbias_type() != AscBiasCorrection::none && part_region.master()))
  {
    attrs |=  PLL_ATTRIB_AB_FLAG;
    attrs |= (unsigned int) model.ascbias_type();
  }

  BasicTree tree(msa.size());
  pll_partition_t * partition = pll_partition_create(
      tree.num_tips(),         /* number of tip sequences */
      tree.num_inner(),        /* number of CLV buffers */
      model.num_states(),      /* number of states in the data */
      part_length,             /* number of alignment sites/patterns */
      model.num_submodels(),   /* number of different substitution models (LG4 = 4) */
      tree.num_branches(),     /* number of probability matrices */
      model.num_ratecats(),    /* number of (GAMMA) rate categories */
      tree.num_inner(),        /* number of scaling buffers */
      attrs                    /* list of flags (SSE3/AVX, TIP-INNER special cases etc.) */
  );

  libpll_check_error("ERROR creating pll_partition");
  assert(partition);

  if (part_region.master() && !model.ascbias_weights().empty())
    pll_set_asc_state_weights(partition, model.ascbias_weights().data());

  if (part_length == part_region.length)
    set_partition_tips(opts, msa, tip_msa_idmap, part_region, partition, model.charmap());
  else
    set_partition_tips(opts, msa, tip_msa_idmap, part_region, partition, model.charmap(), weights);

  assign(partition, model);

  return partition;
}
