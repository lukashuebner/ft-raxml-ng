/*
    Copyright (C) 2017-2019 Alexey Kozlov, Alexandros Stamatakis, Diego Darriba, Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#include <algorithm>
#include <chrono>

#include <memory>

#include "version.h"
#include "common.h"
#include "MSA.hpp"
#include "Options.hpp"
#include "CommandLineParser.hpp"
#include "Optimizer.hpp"
#include "PartitionInfo.hpp"
#include "PartitionedMSAView.hpp"
#include "TreeInfo.hpp"
#include "io/file_io.hpp"
#include "io/binary_io.hpp"
#include "ParallelContext.hpp"
#include "loadbalance/LoadBalancer.hpp"
#include "loadbalance/CoarseLoadBalancer.hpp"
#include "bootstrap/BootstrapGenerator.hpp"
#include "bootstrap/BootstopCheck.hpp"
#include "bootstrap/TransferBootstrapTree.hpp"
#include "bootstrap/ConsensusTree.hpp"
#include "autotune/ResourceEstimator.hpp"
#include "ICScoreCalculator.hpp"
#include "topology/RFDistCalculator.hpp"

#ifdef _RAXML_TERRAPHAST
#include "terraces/TerraceWrapper.hpp"
#endif

using namespace std;

struct RaxmlWorker;

enum class RaxmlRunPhase
{
  start,
  mlsearch,
  bootstrap,
  finish
};

struct RaxmlInstance
{
  Options opts;
  shared_ptr<PartitionedMSA> parted_msa;
  unique_ptr<PartitionedMSA> parted_msa_parsimony;
  map<BranchSupportMetric, shared_ptr<SupportTree> > support_trees;
  shared_ptr<ConsensusTree> consens_tree;

  // When a rank failure occurs, the work has to be redistributed from within
  // the algorithm (TreeInfo for now). This requires callbacks.
  function<PartitionAssignmentList(shared_ptr<vector<double>>)> load_balancer_cb;

  TreeList start_trees;
  BootstrapReplicateList bs_reps;
  TreeList bs_start_trees;

  /* IDs of the trees that have been already inferred (eg after resuming from a checkpoint) */
  IDSet done_ml_trees;
  IDSet done_bs_trees;

  // load balancing
  PartitionAssignmentList proc_part_assign;
  unique_ptr<LoadBalancer> load_balancer;
  unique_ptr<CoarseLoadBalancer> coarse_load_balancer;

  // mini-checkpointing
  // All the IDs of the rank which are the master of at least one model and therefore
  // need to broadcast model parameters.
  shared_ptr<IDVector> ranks_which_are_part_masters; 

  // bootstopping convergence test, only autoMRE is supported for now
  unique_ptr<BootstopCheckMRE> bootstop_checker;
  bool bs_converged;
  RaxmlRunPhase run_phase;

  // mapping taxon name -> tip_id/clv_id in the tree
  NameIdMap tip_id_map;

  // mapping tip_id in the tree (array index) -> sequence index in MSA
  IDVector tip_msa_idmap;

 // unique_ptr<TerraceWrapper> terrace_wrapper;

//  unique_ptr<RandomGenerator> starttree_seed_gen;
//  unique_ptr<RandomGenerator> bootstrap_seed_gen;

  unique_ptr<NewickStream> start_tree_stream;

  /* this is just a dummy random tree used for convenience, e,g, if we need tip labels or
   * just 'any' valid tree for the alignment at hand */
  Tree random_tree;

  /* topological constraint */
  Tree constraint_tree;

  MLTree ml_tree;

  unique_ptr<RFDistCalculator> dist_calculator;
  AncestralStatesSharedPtr ancestral_states;

  vector<RaxmlWorker> workers;
  RaxmlWorker& get_worker() { return workers.at(ParallelContext::local_group_id()); }

  RaxmlInstance() : bs_converged(false), run_phase(RaxmlRunPhase::start) {}
};

struct RaxmlWorker
{
  RaxmlWorker(RaxmlInstance& inst, unsigned int id) :
    instance(inst), worker_id(id) {}

  RaxmlInstance& instance;

  unsigned int worker_id;
//  TreeList start_trees;
//  BootstrapReplicateList bs_reps;
//  TreeList bs_start_trees;

  IDVector start_trees;
  IDVector bs_trees;
  PartitionAssignmentList proc_part_assign;

  size_t total_num_searches() const { return start_trees.size() + bs_trees.size(); }
};

void print_banner()
{
  LOG_INFO << endl << "RAxML-NG v. " << RAXML_VERSION << " released on " << RAXML_DATE <<
      " by The Exelixis Lab." << endl;
  LOG_INFO << "Developed by: Alexey M. Kozlov and Alexandros Stamatakis." << endl;
  LOG_INFO << "Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, "
              "Sarah Lutteropp, Ben Bettisworth." << endl;
  LOG_INFO << "Latest version: https://github.com/amkozlov/raxml-ng" << endl;
  LOG_INFO << "Questions/problems/suggestions? "
              "Please visit: https://groups.google.com/forum/#!forum/raxml" << endl << endl;
}

void init_part_info(RaxmlInstance& instance)
{
  auto& opts = instance.opts;

  instance.parted_msa = std::make_shared<PartitionedMSA>();
  auto& parted_msa = *instance.parted_msa;

  if (!sysutil_file_exists(opts.msa_file))
  {
    throw runtime_error("Alignment file not found: " + opts.msa_file);
  }

  /* check if we have a binary input file */
  if (opts.msa_format == FileFormat::binary ||
      (opts.msa_format == FileFormat::autodetect && RBAStream::rba_file(opts.msa_file)))
  {
    opts.msa_format = FileFormat::binary;

    if (!opts.model_file.empty())
    {
      LOG_WARN <<
          "WARNING: The model you specified on the command line (" << opts.model_file <<
                    ") will be ignored " << endl <<
          "         since the binary MSA file already contains a model definition." << endl <<
          "         If you want to change the model, please re-run RAxML-NG "  << endl <<
          "         with the original PHYLIP/FASTA alignment and --redo option."
          << endl << endl;
    }

    LOG_INFO_TS << "Loading binary alignment from file: " << opts.msa_file << endl;

    auto rba_elem = opts.use_rba_partload ? RBAStream::RBAElement::metadata : RBAStream::RBAElement::all;
    RBAStream bs(opts.msa_file);
    bs >> RBAStream::RBAOutput(parted_msa, rba_elem, nullptr);

    // binary probMSAs are not supported yet
    instance.opts.use_prob_msa = false;

    LOG_INFO_TS << "Alignment comprises " << parted_msa.taxon_count() << " taxa, " <<
        parted_msa.part_count() << " partitions and " <<
        parted_msa.total_length() << " patterns\n" << endl;

    LOG_INFO << parted_msa;

    LOG_INFO << endl;
  }
  /* check if model is a file */
  else if (sysutil_file_exists(opts.model_file))
  {
    // read partition definitions from file
    try
    {
      RaxmlPartitionStream partfile(opts.model_file, ios::in);
      partfile >> parted_msa;
    }
    catch(exception& e)
    {
      throw runtime_error("Failed to read partition file:\n" + string(e.what()));
    }
  }
  else if (!opts.model_file.empty())
  {
    // create and init single pseudo-partition
    parted_msa.emplace_part_info("noname", opts.data_type, opts.model_file);
  }
  else
    throw runtime_error("Please specify an evolutionary model with --model switch");

  assert(parted_msa.part_count() > 0);

  /* make sure that linked branch length mode is set for unpartitioned alignments */
  if (parted_msa.part_count() == 1)
    opts.brlen_linkage = PLLMOD_COMMON_BRLEN_LINKED;

  /* in the scaled brlen mode, use ML optimization of brlen scalers by default */
  if (opts.brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
  {
    for (auto& pinfo: parted_msa.part_list())
      pinfo.model().set_param_mode_default(PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER, ParamValue::ML);
  }

  int freerate_count = 0;

  for (const auto& pinfo: parted_msa.part_list())
  {
    LOG_DEBUG << "|" << pinfo.name() << "|   |" << pinfo.model().to_string() << "|   |" <<
        pinfo.range_string() << "|" << endl;

    if (pinfo.model().ratehet_mode() == PLLMOD_UTIL_MIXTYPE_FREE)
      freerate_count++;
  }

  if (parted_msa.part_count() > 1 && freerate_count > 0 &&
      opts.brlen_linkage == PLLMOD_COMMON_BRLEN_LINKED)
  {
    throw runtime_error("LG4X and FreeRate models are not supported in linked branch length mode.\n"
        "Please use the '--brlen scaled' option to switch into proportional branch length mode.");
  }
}

void print_reduced_msa(const RaxmlInstance& instance, const PartitionedMSAView& reduced_msa_view)
{
  // save reduced MSA and partition files
  auto reduced_msa_fname = instance.opts.output_fname("reduced.phy");
  PhylipStream ps(reduced_msa_fname);

  ps << reduced_msa_view;

  LOG_INFO << "\nNOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) "
              "\nNOTE: was saved to: ";
  LOG_INFO << sysutil_realpath(reduced_msa_fname) << endl;

  // save reduced partition file
  if (sysutil_file_exists(instance.opts.model_file))
  {
    auto reduced_part_fname = instance.opts.output_fname("reduced.partition");
    RaxmlPartitionStream ps(reduced_part_fname, ios::out);

    ps << reduced_msa_view;

    LOG_INFO << "\nNOTE: The corresponding reduced partition file was saved to:\n";
    LOG_INFO << sysutil_realpath(reduced_part_fname) << endl;
  }
}

bool check_msa_global(const MSA& msa)
{
  bool msa_valid = true;

  /* check taxa count */
  if (msa.size() < 4)
  {
    LOG_ERROR << "\nERROR: Your alignment contains less than 4 sequences! " << endl;
    msa_valid = false;
  }

  /* check for duplicate taxon names */
  unsigned long stats_mask = PLLMOD_MSA_STATS_DUP_TAXA;

  pllmod_msa_stats_t * stats = pllmod_msa_compute_stats(msa.pll_msa(),
                                                        4,
                                                        pll_map_nt, // map is not used here
                                                        NULL,
                                                        stats_mask);

  libpll_check_error("ERROR computing MSA stats");
  assert(stats);

  if (stats->dup_taxa_pairs_count > 0)
  {
    LOG_ERROR << endl;
    for (unsigned long c = 0; c < stats->dup_taxa_pairs_count; ++c)
    {
      auto id1 = stats->dup_taxa_pairs[c*2];
      auto id2 = stats->dup_taxa_pairs[c*2+1];
      LOG_ERROR << "ERROR: Sequences " << id1+1 << " and "
                << id2+1 << " have identical name: "
                << msa.label(id1) << endl;
    }
    LOG_ERROR << "\nERROR: Duplicate sequence names found: "
              << stats->dup_taxa_pairs_count << endl;

    msa_valid = false;
  }

  pllmod_msa_destroy_stats(stats);

  return msa_valid;
}

bool check_msa(RaxmlInstance& instance)
{
  LOG_VERB_TS << "Checking the alignment...\n";

  const auto& opts = instance.opts;
  auto& parted_msa = *instance.parted_msa;
  const auto& full_msa = parted_msa.full_msa();
  const auto pll_msa = full_msa.pll_msa();

  bool msa_valid = true;
  bool msa_corrected = false;
  PartitionedMSAView parted_msa_view(instance.parted_msa);

  vector<pair<size_t,size_t> > dup_taxa;
  vector<pair<size_t,size_t> > dup_seqs;
  std::set<size_t> gap_seqs;

  /* check taxon names for invalid characters */
  if (opts.safety_checks.isset(SafetyCheck::msa_names))
  {
    const string invalid_chars = "(),;:' \t\n";
    for (const auto& taxon: parted_msa.taxon_names())
    {
      if (taxon.find_first_of(invalid_chars) != std::string::npos)
      {
        size_t i = 0;
        auto fixed_name = taxon;
        while ((i = fixed_name.find_first_of(invalid_chars, i)) != std::string::npos)
          fixed_name[i++] = '_';
        parted_msa_view.map_taxon_name(taxon, fixed_name);
      }
    }

    msa_valid &= parted_msa_view.taxon_name_map().empty();
  }

  /* check for duplicate sequences */
  if (opts.safety_checks.isset(SafetyCheck::msa_dups))
  {
    unsigned long stats_mask = PLLMOD_MSA_STATS_DUP_SEQS;

    pllmod_msa_stats_t * stats = pllmod_msa_compute_stats(pll_msa,
                                                          4,
                                                          pll_map_nt, // map is not used here
                                                          NULL,
                                                          stats_mask);

    libpll_check_error("ERROR computing MSA stats");
    assert(stats);

    for (unsigned long c = 0; c < stats->dup_seqs_pairs_count; ++c)
    {
      dup_seqs.emplace_back(stats->dup_seqs_pairs[c*2],
                            stats->dup_seqs_pairs[c*2+1]);
    }

    pllmod_msa_destroy_stats(stats);
  }

  size_t total_gap_cols = 0;
  size_t part_num = 0;
  for (auto& pinfo: parted_msa.part_list())
  {
    /* check for invalid MSA characters */
    pllmod_msa_errors_t * errs = pllmod_msa_check(pinfo.msa().pll_msa(),
                                                  pinfo.model().charmap());

    if (errs)
    {
      if (errs->invalid_char_count > 0)
      {
        msa_valid = false;
        LOG_ERROR << endl;
        for (unsigned long c = 0; c < errs->invalid_char_count; ++c)
        {
          auto global_pos = parted_msa.full_msa_site(part_num, errs->invalid_char_pos[c]);
          LOG_ERROR << "ERROR: Invalid character in sequence " <<  errs->invalid_char_seq[c]+1
                    << " at position " <<  global_pos+1  << ": " << errs->invalid_chars[c] << endl;
        }
        part_num++;
        continue;
      }
      pllmod_msa_destroy_errors(errs);
    }
    else
      libpll_check_error("MSA check failed");


    /* Check for all-gap columns and sequences */
    if (opts.safety_checks.isset(SafetyCheck::msa_allgaps))
    {
      unsigned long stats_mask = PLLMOD_MSA_STATS_GAP_SEQS | PLLMOD_MSA_STATS_GAP_COLS;

      pllmod_msa_stats_t * stats = pinfo.compute_stats(stats_mask);

      if (stats->gap_cols_count > 0)
      {
        total_gap_cols += stats->gap_cols_count;
        std::vector<size_t> gap_cols(stats->gap_cols, stats->gap_cols + stats->gap_cols_count);
        pinfo.msa().remove_sites(gap_cols);
  //      parted_msa_view.exclude_sites(part_num, gap_cols);
      }

      std::set<size_t> cur_gap_seq(stats->gap_seqs, stats->gap_seqs + stats->gap_seqs_count);

      if (!part_num)
      {
        gap_seqs = cur_gap_seq;
      }
      else
      {
        for(auto it = gap_seqs.begin(); it != gap_seqs.end();)
        {
          if(cur_gap_seq.find(*it) == cur_gap_seq.end())
            it = gap_seqs.erase(it);
          else
            ++it;
        }
      }

      pllmod_msa_destroy_stats(stats);
    }

    part_num++;
  }

  if (total_gap_cols > 0)
  {
    LOG_WARN << "\nWARNING: Fully undetermined columns found: " << total_gap_cols << endl;
    msa_corrected = true;
  }

  if (!gap_seqs.empty())
  {
   LOG_WARN << endl;
   for (auto c : gap_seqs)
   {
     parted_msa_view.exclude_taxon(c);
     LOG_VERB << "WARNING: Sequence #" << c+1 << " (" << parted_msa.taxon_names().at(c)
              << ") contains only gaps!" << endl;
   }
   LOG_WARN << "WARNING: Fully undetermined sequences found: " << gap_seqs.size() << endl;
  }

  if (!dup_seqs.empty())
  {
    size_t dup_count = 0;
    LOG_WARN << endl;
    for (const auto& p: dup_seqs)
    {
      /* ignore gap-only sequences */
      if (gap_seqs.count(p.first) || gap_seqs.count(p.second))
        continue;

      ++dup_count;
      parted_msa_view.exclude_taxon(p.second);
      LOG_WARN << "WARNING: Sequences " << parted_msa.taxon_names().at(p.first) << " and " <<
          parted_msa.taxon_names().at(p.second) << " are exactly identical!" << endl;
    }
    if (dup_count > 0)
      LOG_WARN << "WARNING: Duplicate sequences found: " << dup_count << endl;
  }

  if (!instance.opts.nofiles_mode && (msa_corrected || !parted_msa_view.identity()))
  {
    print_reduced_msa(instance, parted_msa_view);
  }

  if (!parted_msa_view.taxon_name_map().empty())
  {
    LOG_ERROR << endl;
    for (auto it: parted_msa_view.taxon_name_map())
      LOG_ERROR << "ERROR: Following taxon name contains invalid characters: " << it.first << endl;

    LOG_ERROR << endl;
    LOG_INFO << "NOTE: Following symbols are not allowed in taxa names to ensure Newick compatibility:\n"
                "NOTE: \" \" (space), \";\" (semicolon), \":\" (colon), \",\" (comma), "
                       "\"()\" (parentheses), \"'\" (quote). " << endl;
    LOG_INFO << "NOTE: Please either correct the names manually, or use the reduced alignment file\n"
                "NOTE: generated by RAxML-NG (see above).";
    LOG_INFO << endl;
  }

  return msa_valid;
}

size_t total_free_params(const RaxmlInstance& instance)
{
  const auto& parted_msa = *instance.parted_msa;
  size_t free_params = parted_msa.total_free_model_params();
  size_t num_parts = parted_msa.part_count();
  auto tree = BasicTree(parted_msa.taxon_count());
  auto num_branches = tree.num_branches();
  auto brlen_linkage = instance.opts.brlen_linkage;

  if (brlen_linkage == PLLMOD_COMMON_BRLEN_LINKED)
    free_params += num_branches;
  else if (brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
    free_params += num_branches + num_parts - 1;
  else if (brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
    free_params += num_branches * num_parts;

  return free_params;
}

void check_models(const RaxmlInstance& instance)
{
  const auto& opts = instance.opts;
  bool zero_freqs = false;

  for (const auto& pinfo: instance.parted_msa->part_list())
  {
    auto stats = pinfo.stats();
    auto model = pinfo.model();

    // check for non-recommended model combinations
    if (opts.safety_checks.isset(SafetyCheck::model_lg4_freqs))
    {
      if ((model.name() == "LG4X" || model.name() == "LG4M") &&
          model.param_mode(PLLMOD_OPT_PARAM_FREQUENCIES) != ParamValue::model)
      {
        throw runtime_error("Partition \"" + pinfo.name() +
                            "\": You specified LG4M or LG4X model with shared stationary based frequencies (" +
                            model.to_string(false) + ").\n"
                            "Please be warned, that this is against the idea of LG4 models and hence it's not recommended!" + "\n"
                            "If you know what you're doing, you can add --force command line switch to disable this safety check.");
      }
    }

    // check for zero state frequencies
    if (opts.safety_checks.isset(SafetyCheck::model_zero_freqs))
    {
      if (model.param_mode(PLLMOD_OPT_PARAM_FREQUENCIES) == ParamValue::empirical)
      {
        const auto& freqs = stats.emp_base_freqs;
        for (unsigned int i = 0; i < freqs.size(); ++i)
        {
          if (freqs[i] < PLL_EIGEN_MINFREQ)
          {
            if (!zero_freqs)
            {
              LOG_WARN << endl;
              zero_freqs = true;
            }

            LOG_WARN << "WARNING: State " << to_string(i) <<
                (instance.parted_msa->part_count() > 1 ? " in partition " + pinfo.name() : "") <<
                " has very low frequency (" << FMT_PREC9(freqs[i]) << ")!" << endl;

            LOG_VERB << "Base frequencies: ";
            for (unsigned int j = 0; j < freqs.size(); ++j)
              LOG_VERB << FMT_PREC9(freqs[j]) <<  " ";
            LOG_VERB << endl;
          }
        }
      }
    }

    // check for user-defined state frequencies which do not sum up to one
    if (opts.safety_checks.isset(SafetyCheck::model_invalid_freqs))
    {
      if (model.param_mode(PLLMOD_OPT_PARAM_FREQUENCIES) == ParamValue::user)
      {
        const auto& freqs = model.base_freqs(0);
        double sum = 0.;
        for (unsigned int i = 0; i < freqs.size(); ++i)
          sum += freqs[i];

        if (fabs(sum - 1.0) > 0.01)
        {
          LOG_ERROR << "\nBase frequencies: ";
          for (unsigned int j = 0; j < freqs.size(); ++j)
            LOG_ERROR << FMT_PREC9(freqs[j]) <<  " ";
          LOG_ERROR << endl;

          throw runtime_error("User-specified stationary base frequencies"
                              " in partition " + pinfo.name() + " do not sum up to 1.0!\n"
                              "Please provide normalized frequencies.");
        }
      }
    }

    if (model.num_submodels() > 1 &&
        (model.param_mode(PLLMOD_OPT_PARAM_FREQUENCIES) == ParamValue::ML ||
         model.param_mode(PLLMOD_OPT_PARAM_SUBST_RATES) == ParamValue::ML))
    {
      throw runtime_error("Invalid model " + model.to_string(false) + " in partition " + pinfo.name() + ":\n"
                          "Mixture models with ML estimates of rates/frequencies are not supported yet!");
    }

    // check partitions which contain invariant sites and have ascertainment bias enabled
    if (opts.safety_checks.isset(SafetyCheck::model_asc_bias))
    {
      if (model.ascbias_type() != AscBiasCorrection::none && stats.inv_count() > 0)
      {
        throw runtime_error("You enabled ascertainment bias correction for partition " +
                             pinfo.name() + ", but it contains " +
                             to_string(stats.inv_count()) + " invariant sites.\n"
                            "This is not allowed! Please either remove invariant sites or "
                            "disable ascertainment bias correction.");
      }
    }
  }

  if (zero_freqs)
  {
    LOG_WARN << endl << "WARNING: Some states have very low frequencies, "
        "which might lead to numerical issues!" << endl;
  }

  /* Check for extreme cases of overfitting (K >= n) */
  if (opts.safety_checks.isset(SafetyCheck::model_overfit))
  {
    if (instance.parted_msa->part_count() > 1)
    {
      size_t model_free_params = instance.parted_msa->total_free_model_params();
      size_t free_params = total_free_params(instance);
      size_t sample_size = instance.parted_msa->total_sites();
      string errmsg = "Number of free parameters (K=" + to_string(free_params) +
          ") is larger than alignment size (n=" + to_string(sample_size) + ").\n" +
          "       This might lead to overfitting and compromise tree inference results!\n" +
          "       Please consider revising your partitioning scheme, conducting formal model selection\n" +
          "       and/or using linked/scaled branch lengths across partitions.\n" +
          "NOTE:  You can disable this check by adding the --force option.\n";

      if (free_params >= sample_size)
      {
        if (model_free_params >= sample_size ||
            instance.opts.brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
        {
          throw runtime_error(errmsg);
        }
        else
          LOG_WARN << endl << "WARNING: " << errmsg << endl;
      }
    }
  }
}

void check_tree(const PartitionedMSA& msa, const Tree& tree)
{
  auto missing_taxa = 0;
  auto duplicate_taxa = 0;

  if (msa.taxon_count() > tree.num_tips())
    throw runtime_error("Alignment file contains more sequences than expected");
  else if (msa.taxon_count() != tree.num_tips())
    throw runtime_error("Some taxa are missing from the alignment file");

  unordered_set<string> tree_labels;
  unordered_set<string> msa_labels(msa.taxon_names().cbegin(), msa.taxon_names().cend());

  for (const auto& tip: tree.tip_labels())
  {
    if (!tree_labels.insert(tip.second).second)
    {
      LOG_ERROR << "ERROR: Taxon name appears more than once in the tree: " << tip.second << endl;
      duplicate_taxa++;
    }

    if (msa_labels.count(tip.second) == 0)
    {
      LOG_ERROR << "ERROR: Taxon name not found in the alignment: " << tip.second << endl;
      missing_taxa++;
    }
  }

  if (duplicate_taxa > 0)
    throw runtime_error("Tree contains duplicate taxon names (see above)!");

  if (missing_taxa > 0)
    throw runtime_error("Please check that sequence labels in the alignment and in the tree file are identical!");

  /* check for negative branch length */
  for (const auto& branch: tree.topology())
  {
    if (branch.length < 0.)
      throw runtime_error("Tree file contains negative branch lengths!");
  }
}

/* This function is called after MPI init, but before thread init, MSA loading etc. */
void check_options_early(Options& opts)
{
  auto num_procs = opts.num_ranks * opts.num_threads;

  if (!opts.constraint_tree_file.empty() &&
      (opts.start_trees.count(StartingTree::parsimony) > 0 ||
       opts.start_trees.count(StartingTree::user)))
  {
    throw runtime_error(" User and parsimony starting trees are not supported in combination with "
                        "constrained tree inference.\n"
                        "       Please use random starting trees instead.");
  }

  if (opts.num_workers > num_procs)
  {
    throw OptionException("The specified number of parallel tree searches (" +
                          to_string(opts.num_workers) +
                          ") is higher than the number of available threads (" +
                          to_string(num_procs) + ")");
  }

  /* check for unsupported coarse-grained topology */
  if (opts.coarse() && opts.num_ranks > opts.num_workers)
  {
    throw runtime_error("Unsupported parallelization topology!\n"
                        "NOTE:  Multiple MPI ranks per worker are not allowed in coarse-grained mode.\n");
  }

  if (opts.command == Command::ancestral)
  {
    if (opts.use_pattern_compression)
      throw runtime_error("Pattern compression is not supported in ancestral state reconstruction mode!");
    if (opts.use_repeats)
      throw runtime_error("Site repeats are not supported in ancestral state reconstruction mode!");
    if (opts.use_rate_scalers)
      throw runtime_error("Per-rate scalers are not supported in ancestral state reconstruction mode!");
    if (opts.num_ranks > 1)
      throw runtime_error("MPI parallelization is not supported in ancestral state reconstruction mode!");
  }

  /* autodetect if we can use partial RBA loading */
  opts.use_rba_partload &= (opts.num_ranks > 1 && !opts.coarse());                // only useful for fine-grain MPI runs
  opts.use_rba_partload &= (!opts.start_trees.count(StartingTree::parsimony));    // does not work with parsimony
  opts.use_rba_partload &= (opts.command == Command::search ||                    // currently doesn't work with bootstrap
                            opts.command == Command::evaluate ||
                            opts.command == Command::ancestral);

  LOG_DEBUG << "RBA partial loading: " << (opts.use_rba_partload ? "ON" : "OFF") << endl;
}

void check_options(RaxmlInstance& instance)
{
  const auto& opts = instance.opts;

  /* check that all outgroup taxa are present in the alignment */
  if (!opts.outgroup_taxa.empty())
  {
    NameList missing_taxa;
    for (const auto& ot: opts.outgroup_taxa)
    {
      if (!instance.parted_msa->taxon_id_map().count(ot))
        missing_taxa.push_back(ot);
    }

    if (!missing_taxa.empty())
    {
      LOG_ERROR << "ERROR: Following taxa were specified as an outgroup "
                                                     "but are missing from the alignment:" << endl;
      for (const auto& mt: missing_taxa)
        LOG_ERROR << mt << endl;
      LOG_ERROR << endl;
      throw runtime_error("Outgroup taxon not found.");
    }
  }

  /* check that we have enough patterns per thread */
  if (opts.safety_checks.isset(SafetyCheck::perf_threads))
  {
    if (ParallelContext::master_rank() && ParallelContext::num_procs() > 1)
    {
      StaticResourceEstimator resEstimator(*instance.parted_msa, instance.opts);
      auto res = resEstimator.estimate();
      if (ParallelContext::threads_per_group() > res.num_threads_response)
      {
        LOG_WARN << endl;
        LOG_WARN << "WARNING: You might be using too many threads (" << ParallelContext::num_procs()
                 <<  ") for your alignment with "
                 << (opts.use_pattern_compression ?
                        to_string(instance.parted_msa->total_patterns()) + " unique patterns." :
                        to_string(instance.parted_msa->total_sites()) + " alignment sites.")
                 << endl;
        LOG_WARN << "NOTE:    For the optimal throughput, please consider using fewer threads " << endl;
        LOG_WARN << "NOTE:    and parallelize across starting trees/bootstrap replicates." << endl;
        LOG_WARN << "NOTE:    As a general rule-of-thumb, please assign at least 200-1000 "
            "alignment patterns per thread." << endl << endl;

        if (ParallelContext::threads_per_group() > 2 * res.num_threads_response)
        {
          throw runtime_error("Too few patterns per thread! "
                              "RAxML-NG will terminate now to avoid wasting resources.\n"
                              "NOTE:  Please reduce the number of threads (see guidelines above).\n"
                              "NOTE:  This check can be disabled with the '--force' option.");
        }
      }
    }
  }

  if (instance.parted_msa->taxon_count() > RAXML_RATESCALERS_TAXA &&
      !instance.opts.use_rate_scalers && opts.command != Command::ancestral)
  {
    LOG_INFO << "\nNOTE: Per-rate scalers were automatically enabled to prevent numerical issues "
        "on taxa-rich alignments." << endl;
    LOG_INFO << "NOTE: You can use --force switch to skip this check "
        "and fall back to per-site scalers." << endl << endl;
    instance.opts.use_rate_scalers = true;
  }

  /* make sure we do not check for convergence too often in coarse-grained parallelization mode */
  instance.opts.bootstop_interval = std::max(opts.bootstop_interval, opts.num_workers*2);
}

void load_msa(RaxmlInstance& instance)
{
  const auto& opts = instance.opts;
  auto& parted_msa = *instance.parted_msa;

  LOG_INFO_TS << "Reading alignment from file: " << opts.msa_file << endl;

  /* load MSA */
  auto msa = msa_load_from_file(opts.msa_file, opts.msa_format);

  LOG_INFO_TS << "Loaded alignment with " << msa.size() << " taxa and " <<
      msa.num_sites() << " sites" << endl;

  if (msa.probabilistic() && opts.use_prob_msa)
  {
    instance.opts.use_pattern_compression = false;
    instance.opts.use_tip_inner = false;

    if (parted_msa.part_count() > 1)
      throw runtime_error("Partitioned probabilistic alignments are not supported yet, sorry...");
  }
  else
    instance.opts.use_prob_msa = false;

  if (!check_msa_global(msa))
    throw runtime_error("Alignment check failed (see details above)!");

  parted_msa.full_msa(std::move(msa));

  LOG_VERB_TS << "Extracting partitions... " << endl;

  parted_msa.split_msa();

  /* check alignment */
  if (!check_msa(instance))
    throw runtime_error("Alignment check failed (see details above)!");

  if (opts.use_pattern_compression)
  {
    LOG_VERB_TS << "Compressing alignment patterns... " << endl;
    parted_msa.compress_patterns();

    // temp workaround: since MSA pattern compression calls rand(), it will change all random
    // numbers generated afterwards. so just reset seed to the initial value to ensure that
    // starting trees, BS replicates etc. are the same regardless whether pat.comp is ON or OFF
    srand(opts.random_seed);
  }

//  if (parted_msa.part_count() > 1)
//    instance.terrace_wrapper.reset(new TerraceWrapper(parted_msa));

  parted_msa.set_model_empirical_params();

  check_models(instance);

  LOG_INFO << endl;

  LOG_INFO << "Alignment comprises " << parted_msa.part_count() << " partitions and "
           << parted_msa.total_length() << (opts.use_pattern_compression ? " patterns" : " sites")
           << endl << endl;

  LOG_INFO << parted_msa;

  LOG_INFO << endl;

  if (ParallelContext::master_rank() &&
      !instance.opts.use_prob_msa && !instance.opts.binary_msa_file().empty())
  {
    auto binary_msa_fname = instance.opts.binary_msa_file();
    if (sysutil_file_exists(binary_msa_fname) && !opts.redo_mode &&
        opts.command != Command::parse)
    {
      LOG_INFO << "NOTE: Binary MSA file already exists: " << binary_msa_fname << endl << endl;
    }
    else if (opts.command != Command::check)
    {
      RBAStream bs(binary_msa_fname);
      bs << parted_msa;
      LOG_INFO << "NOTE: Binary MSA file created: " << binary_msa_fname << endl << endl;
    }
  }
}

void load_parted_msa(RaxmlInstance& instance)
{
  init_part_info(instance);

  assert(instance.parted_msa);

  if (instance.opts.msa_format != FileFormat::binary)
    load_msa(instance);

  // use MSA sequences IDs as "normalized" tip IDs in all trees
  instance.tip_id_map = instance.parted_msa->taxon_id_map();
}

void prepare_tree(const RaxmlInstance& instance, Tree& tree)
{
  /* fix missing branch lengths */
  tree.fix_missing_brlens();

  /* make sure tip indices are consistent between MSA and pll_tree */
  assert(!instance.parted_msa->taxon_id_map().empty());
  tree.reset_tip_ids(instance.tip_id_map);
}

Tree generate_tree(const RaxmlInstance& instance, StartingTree type)
{
  Tree tree;

  const auto& opts = instance.opts;
  const auto& parted_msa = *instance.parted_msa;
  const auto  tree_rand_seed = rand();

  switch (type)
  {
    case StartingTree::user:
    {
      assert(instance.start_tree_stream);

      /* parse the unrooted binary tree in newick format, and store the number
         of tip nodes in tip_nodes_count */
      *instance.start_tree_stream >> tree;

      LOG_DEBUG << "Loaded user starting tree with " << tree.num_tips() << " taxa from: "
                           << opts.tree_file << endl;

      check_tree(parted_msa, tree);

      break;
    }
    case StartingTree::random:
      /* no starting tree provided, generate a random one */

      LOG_DEBUG << "Generating a random starting tree with " << parted_msa.taxon_count()
                << " taxa" << endl;

      if (instance.constraint_tree.empty())
        tree = Tree::buildRandom(parted_msa.taxon_names(), tree_rand_seed);
      else
        tree = Tree::buildRandomConstrained(parted_msa.taxon_names(), tree_rand_seed,
                                            instance.constraint_tree);

      break;
    case StartingTree::parsimony:
    {
      LOG_DEBUG << "Generating a parsimony starting tree with " << parted_msa.taxon_count()
                << " taxa" << endl;

      unsigned int score;
      unsigned int attrs = opts.simd_arch;

      // TODO: check if there is any reason not to use tip-inner
      attrs |= PLL_ATTRIB_PATTERN_TIP;

      const PartitionedMSA& pars_msa = instance.parted_msa_parsimony ?
                                    *instance.parted_msa_parsimony.get() : *instance.parted_msa;
      tree = Tree::buildParsimony(pars_msa, tree_rand_seed, attrs, &score);

      LOG_DEBUG << "Parsimony score of the starting tree: " << score << endl;

      break;
    }
    default:
      sysutil_fatal("Unknown starting tree type: %d\n", type);
  }

  assert(!tree.empty());

  prepare_tree(instance, tree);

  return tree;
}

void load_start_trees(RaxmlInstance& instance)
{
  NewickStream ts(instance.opts.start_tree_file(), std::ios::in);
  size_t i = 0;
  while (ts.peek() != EOF)
  {
    Tree tree;
    ts >> tree;
    i++;

    prepare_tree(instance, tree);
    instance.start_trees.emplace_back(tree);
  }

  if (instance.opts.start_trees.count(StartingTree::user) > 0)
  {
    // in case of user starting trees, we do not know num_searches
    // until we read trees from the file. that's why we update num_searches here.
    assert(i >= instance.opts.num_searches);
    instance.opts.num_searches = i;
  }
  else
    assert(i == instance.opts.num_searches);
}

void load_checkpoint(RaxmlInstance& instance, CheckpointManager& cm)
{
  /* init checkpoint and set to the manager */
  cm.init_checkpoints(instance.random_tree, instance.parted_msa->models());

  auto& ckpfile = cm.checkp_file();

  if (!instance.opts.redo_mode)
  {
    if (cm.read())
    {
      // read start trees from file to avoid re-generation
      // NOTE: doesn't work for constrained tree search
      if (sysutil_file_exists(instance.opts.start_tree_file()) &&
          instance.opts.num_searches > 0 &&
          instance.opts.constraint_tree_file.empty())
      {
        load_start_trees(instance);
      }
    }

    /* collect all trees inferred so far at master rank */
    cm.gather_ml_trees();
    cm.gather_bs_trees();

    // NB: consider BS trees from the previous run when performing bootstopping test
    if (instance.bootstop_checker)
    {
      auto bs_tree = instance.random_tree;
      for (auto it: ckpfile.bs_trees)
      {
        bs_tree.topology(it.second.second);

        instance.bootstop_checker->add_bootstrap_tree(bs_tree);
      }
    }

    /* determine if we are in ML tree search or in bootstrapping phase */
    instance.run_phase = (ckpfile.ml_trees.size() >= instance.opts.num_searches) ?
        RaxmlRunPhase::bootstrap : RaxmlRunPhase::mlsearch;

    if (ParallelContext::num_ranks() > 1)
    {
      ParallelContext::mpi_broadcast(instance.run_phase);

      /* gather in-progress tree ids */
      auto& in_work_trees = instance.run_phase == RaxmlRunPhase::bootstrap ?
          instance.done_bs_trees : instance.done_ml_trees;
      for (auto& c: ckpfile.checkp_list)
      {
        if (c.tree_index > 0)
          in_work_trees.insert(c.tree_index);
      }

      auto worker_cb = [&in_work_trees](void * buf, size_t buf_size) -> int
          {
            return (int) BinaryStream::serialize((char*) buf, buf_size, in_work_trees);
          };

      /* receive callback -> master rank */
      auto master_cb = [&in_work_trees](void * buf, size_t buf_size)
         {
           BinaryStream bs((char*) buf, buf_size);

           IDSet recv_trees;
           bs >> recv_trees;
           in_work_trees.insert(recv_trees.cbegin(), recv_trees.cend());
         };

      ParallelContext::mpi_gather_custom(worker_cb, master_cb);

      if (ParallelContext::master())
      {
        for (const auto& p: ckpfile.ml_trees)
          instance.done_ml_trees.insert(p.first);

        for (const auto& p: ckpfile.bs_trees)
          instance.done_bs_trees.insert(p.first);
      }

      ParallelContext::global_mpi_barrier();

      // broadcast done_start_trees + done_bs_trees
      ParallelContext::mpi_broadcast(instance.done_ml_trees);
      ParallelContext::mpi_broadcast(instance.done_bs_trees);
    }

    if (!instance.done_ml_trees.empty())
    {
      LOG_INFO_TS << "NOTE: Resuming execution from checkpoint " <<
          "(logLH: " << cm.checkpoint().loglh() <<
          ", ML trees: " << ckpfile.ml_trees.size() <<
          ", bootstraps: " << ckpfile.bs_trees.size() <<
          ")"
          << endl;
    }

//    printf("DONE: %lu %lu\n", ParallelContext::rank_id(), instance.done_ml_trees.size());
  }
}

void load_constraint(RaxmlInstance& instance)
{
  const auto& parted_msa = *instance.parted_msa;
  const auto& opts = instance.opts;

  if (!instance.opts.constraint_tree_file.empty())
  {
    if (!sysutil_file_exists(opts.constraint_tree_file))
      throw runtime_error("Constraint tree file not found: " + opts.constraint_tree_file);

    NewickStream nw_cons(instance.opts.constraint_tree_file, std::ios::in);
    Tree& cons_tree = instance.constraint_tree;
    nw_cons >> cons_tree;

    LOG_INFO_TS << "Loaded " <<
        (cons_tree.num_tips() == parted_msa.taxon_count() ? "" : "non-") <<
        "comprehensive constraint tree with " << cons_tree.num_tips() << " taxa" << endl;

    // check if taxa names are consistent between contraint tree and MSA
    {
      NameList missing_taxa;
      for (const auto& l: cons_tree.tip_labels())
      {
        if (!parted_msa.taxon_id_map().count(l.second))
          missing_taxa.push_back(l.second);;
      }

      if (!missing_taxa.empty())
      {
        stringstream ss;
        ss << "Following " << missing_taxa.size() <<
            " taxa present in the constraint tree can not be found in the alignment: " << endl;
        for (const auto& taxon: missing_taxa)
          ss << taxon << endl;
        throw runtime_error(ss.str());
      }
    }

//    /* branch lengths in constraint tree can be misleading, so just reset them to default value */
//    cons_tree.reset_brlens();

    if (cons_tree.num_tips() < parted_msa.taxon_count())
    {
      // incomplete constraint tree -> adjust tip IDs such that all taxa in the constraint tree
      // go before the remaining free taxa
      instance.tip_id_map.clear();
      instance.tip_msa_idmap.resize(parted_msa.taxon_count());
      auto cons_name_map = cons_tree.tip_ids();
      size_t seq_id = 0;
      size_t cons_tip_id = 0;
      size_t free_tip_id = cons_tree.num_tips();
      for (const auto& name: parted_msa.taxon_names())
      {
        auto tip_id = cons_name_map.count(name) ? cons_tip_id++ : free_tip_id++;
        instance.tip_id_map[name] = tip_id;
        instance.tip_msa_idmap[tip_id] = seq_id++;
      }
      assert(cons_tip_id == cons_tree.num_tips());
      assert(free_tip_id == instance.tip_id_map.size());
      assert(instance.tip_id_map.size() == parted_msa.taxon_count());
    }
    else if (cons_tree.binary() && !instance.opts.force_mode)
    {
      throw runtime_error("You provided a comprehensive, fully-resolved tree as a topological constraint.\n"
          "Since this is almost certainly not what you intended, RAxML-NG will now exit...");
    }

    /* make sure tip indices are consistent between MSA and pll_tree */
    cons_tree.reset_tip_ids(instance.tip_id_map);

//    pll_utree_show_ascii(&cons_tree.pll_utree_root(), PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_BRANCH_LENGTH |
//                                     PLL_UTREE_SHOW_CLV_INDEX );
  }
}

void build_parsimony_msa(RaxmlInstance& instance)
{
  // create 1 partition per datatype
  const PartitionedMSA& orig_msa = *instance.parted_msa;

  instance.parted_msa_parsimony.reset(new PartitionedMSA(orig_msa.taxon_names()));
  PartitionedMSA& pars_msa = *instance.parted_msa_parsimony.get();

  NameIdMap datatype_pinfo_map;
  for (const auto& pinfo: orig_msa.part_list())
  {
    const auto& model = pinfo.model();
    auto data_type_name = model.data_type_name();

    auto iter = datatype_pinfo_map.find(data_type_name);
    if (iter == datatype_pinfo_map.end())
    {
      pars_msa.emplace_part_info(data_type_name, model.data_type(), model.to_string());
      auto& pars_pinfo = pars_msa.part_list().back();
      pars_pinfo.msa(MSA(pinfo.msa().num_sites()));
      datatype_pinfo_map[data_type_name] = pars_msa.part_count()-1;
    }
    else
    {
      auto& msa = pars_msa.part_list().at(iter->second).msa();
      msa.num_sites(msa.num_sites() + pinfo.msa().num_sites());
    }
  }

  // set_per-datatype MSA
  for (size_t j = 0; j < orig_msa.taxon_count(); ++j)
  {
    for (auto& pars_pinfo: pars_msa.part_list())
    {
      auto pars_datatype = pars_pinfo.model().data_type_name();
      std::string sequence;
      sequence.resize(pars_pinfo.msa().num_sites());
      size_t offset = 0;

      for (const auto& pinfo: orig_msa.part_list())
      {
        // different datatype -> skip for now
        if (pinfo.model().data_type_name() != pars_datatype)
          continue;

        const auto w = pinfo.msa().weights();
        const auto s = pinfo.msa().at(j);

        if (w.empty())
        {
          for (size_t k = 0; k < s.size(); ++k)
            sequence[offset++] = s[k];
        }
        else
        {
          for (size_t k = 0; k < w.size(); ++k)
          {
            auto wk = w[k];
            while(wk-- > 0)
              sequence[offset++] = s[k];
          }
        }
      }

      assert(offset == sequence.size());

      pars_pinfo.msa().append(sequence);
    }
  }

  // compress patterns
  if (instance.opts.use_pattern_compression)
  {
    for (auto& pinfo: pars_msa.part_list())
    {
      pinfo.compress_patterns();
    }
  }
}

void build_start_trees(RaxmlInstance& instance, size_t skip_trees)
{
  auto& opts = instance.opts;
  const auto& parted_msa = *instance.parted_msa;

  /* all start trees were already generated/loaded -> return */
  if (skip_trees + instance.start_trees.size() >= instance.opts.num_searches)
    return;

  for (auto& st_tree: opts.start_trees)
  {
    auto st_tree_type = st_tree.first;
    auto& st_tree_count = st_tree.second;
    switch (st_tree_type)
    {
      case StartingTree::user:
        LOG_INFO_TS << "Loading user starting tree(s) from: " << opts.tree_file << endl;
        if (!sysutil_file_exists(opts.tree_file))
          throw runtime_error("File not found: " + opts.tree_file);
        instance.start_tree_stream.reset(new NewickStream(opts.tree_file, std::ios::in));
        break;
      case StartingTree::random:
        LOG_INFO_TS << "Generating " << st_tree_count << " random starting tree(s) with "
                    << parted_msa.taxon_count() << " taxa" << endl;
        break;
      case StartingTree::parsimony:
        if (parted_msa.part_count() > 1)
        {
          LOG_DEBUG_TS << "Generating MSA partitioned by data type for parsimony computation" << endl;
          build_parsimony_msa(instance);
        }
        LOG_INFO_TS << "Generating " << st_tree_count << " parsimony starting tree(s) with "
                    << parted_msa.taxon_count() << " taxa" << endl;
        break;
      default:
        assert(0);
    }

    for (size_t i = 0; i < st_tree_count; ++i)
    {
      auto tree = generate_tree(instance, st_tree_type);

      // TODO use universal starting tree generator
      if (st_tree_type == StartingTree::user)
      {
        if (instance.start_tree_stream->peek() != EOF)
        {
          st_tree_count++;
          opts.num_searches++;
        }
      }

      if (skip_trees > 0)
      {
        skip_trees--;
        continue;
      }

      instance.start_trees.emplace_back(tree);
    }
  }

  // free memory used for parsimony MSA
  instance.parted_msa_parsimony.release();

  if (::ParallelContext::master_rank())
  {
    NewickStream nw_start(opts.start_tree_file());
    for (auto const& tree: instance.start_trees)
      nw_start << tree;
  }
}

// This will build a function which will recalculate the load balancing based on the current
// number of ranks. This can be helpful if the number of ranks changed due to a rank failure.
std::function<PartitionAssignmentList(shared_ptr<vector<double>>)> build_load_balancer(RaxmlInstance& instance) {
  return [&instance] (shared_ptr<vector<double>> work_by_partition = nullptr) -> PartitionAssignmentList {
    PartitionAssignment part_sizes;

    /* init list of partition sizes */
    for (size_t i = 0; i < instance.parted_msa->part_count(); i++)
    {
      auto const& pinfo = instance.parted_msa->part_list()[i];
      if (work_by_partition == nullptr) {
        part_sizes.assign_sites(i, 0, pinfo.length(), pinfo.model().clv_entry_size());
      } else {
        assert(work_by_partition->size() == instance.parted_msa->part_count());
        assert(work_by_partition->size() > i);
        assert(work_by_partition->at(i) > 0);
        assert(pinfo.length() > 0);
        double weight = work_by_partition->at(i) / pinfo.length();
        part_sizes.assign_sites(i, 0, pinfo.length(), weight);
      }
    }

    return instance.load_balancer->get_all_assignments(part_sizes, ParallelContext::threads_per_group());
  };
}

void balance_load(RaxmlInstance& instance, std::function<PartitionAssignmentList()> load_balancer)
{
  instance.proc_part_assign = load_balancer();
  LOG_INFO_TS << "Data distribution: " << PartitionAssignmentStats(instance.proc_part_assign) << endl;
  LOG_VERB << endl << instance.proc_part_assign;
}

// TODO: Failure mitigation in case of bootstrap replicates (different load balancer)
PartitionAssignmentList balance_load(RaxmlInstance& instance, WeightVectorList part_site_weights)
{
  /* This function is used to re-distribute sites across processes for each bootstrap replicate.
   * Since during bootstrapping alignment sites are sampled with replacement, some sites will be
   * absent from BS alignment. Therefore, site distribution computed for original alignment can
   * be suboptimal for BS replicates. Here, we recompute the site distribution, ignoring all sites
   * that are not present in BS replicate (i.e., have weight of 0 in part_site_weights).
   * */

  PartitionAssignmentList assign_list;
  PartitionAssignment part_sizes;
  WeightVectorList comp_pos_map(part_site_weights.size());

  /* init list of partition sizes */
  size_t i = 0;
  for (auto const& weights: part_site_weights)
  {
    /* build mapping from compressed indices to the original/uncompressed ones */
    comp_pos_map[i].reserve(weights.size());
    for (size_t s = 0; s < weights.size(); ++s)
    {
      if (weights[s] > 0)
      {
        comp_pos_map[i].push_back(s);
      }
    }

    LOG_DEBUG << "Partition #" << i << ": " << comp_pos_map[i].size() << endl;

    /* add compressed partition length to the */
    part_sizes.assign_sites(i, 0, comp_pos_map[i].size(),
                            instance.parted_msa->model(i).clv_entry_size());
    ++i;
  }

  assign_list = instance.load_balancer->get_all_assignments(part_sizes,
                                                            ParallelContext::threads_per_group());

  LOG_VERB_TS << "Data distribution: " << PartitionAssignmentStats(assign_list) << endl;
  LOG_DEBUG << endl << assign_list;

  // translate partition range coordinates: compressed -> uncompressed
  for (auto& part_assign: assign_list)
  {
    for (auto& part_range: part_assign)
    {
      const auto& pos_map = comp_pos_map[part_range.part_id];
      const auto comp_start = part_range.start;
      const auto comp_end = comp_start + part_range.length - 1;

      part_range.start = (comp_start > 0) ? pos_map[comp_start] : 0;
      part_range.length = pos_map[comp_end] - part_range.start + 1;
    }
  }

//  LOG_VERB_TS << "(uncompressed) Data distribution: " << PartitionAssignmentStats(instance.proc_part_assign) << endl;
//  LOG_DEBUG << endl << instance.proc_part_assign;

  return assign_list;
}

void balance_load_coarse(RaxmlInstance& instance, const CheckpointFile& ckpfile)
{
  auto num_workers = ParallelContext::num_groups();

  CoarseAssignment todo_start_trees, todo_bs_trees;
  for (size_t i = 1; i <= instance.start_trees.size(); ++i)
  {
    if (!instance.done_ml_trees.count(i))
      todo_start_trees.push_back(i);
  }

  for (size_t i = 1; i <= instance.bs_start_trees.size(); ++i)
  {
    if (!instance.done_bs_trees.count(i))
      todo_bs_trees.push_back(i);
  }

  /* distribute ML and BS tree searches */
  CoarseAssignmentList start_tree_assign = instance.coarse_load_balancer->get_all_assignments(todo_start_trees, num_workers);
  CoarseAssignmentList bs_tree_assign = instance.coarse_load_balancer->get_all_assignments(todo_bs_trees, num_workers);

  assert(instance.workers.size() == ckpfile.checkp_list.size());
  for (size_t i = 0; i < instance.workers.size(); ++i)
  {
    auto& wrk = instance.workers[i];
    wrk.start_trees = start_tree_assign.at(wrk.worker_id);
    wrk.bs_trees = bs_tree_assign.at(wrk.worker_id);

    /* add current tree from a checkpoint */
    auto& ckp = ckpfile.checkp_list[i];
    auto& in_work_trees = instance.run_phase == RaxmlRunPhase::bootstrap ? wrk.bs_trees : wrk.start_trees;
    if (ckp.tree_index > 0)
      in_work_trees.insert(in_work_trees.begin(), ckp.tree_index);
  }

  LOG_INFO_TS << "Data distribution: max. searches per worker: "
              << instance.workers.at(0).total_num_searches() << endl;
}

void generate_bootstraps(RaxmlInstance& instance, const CheckpointFile& checkp)
{
  if (instance.opts.command == Command::bootstrap || instance.opts.command == Command::all ||
      instance.opts.command == Command::bsmsa)
  {
    assert(instance.parted_msa);

    /* generate replicate alignments */
    BootstrapGenerator bg;
    for (size_t b = 0; b < instance.opts.num_bootstraps; ++b)
    {
      auto seed = rand();

//      /* check if this BS was already computed in the previous run and saved in checkpoint */
//      if (b < checkp.bs_trees.size())
//        continue;

      instance.bs_reps.emplace_back(bg.generate(*instance.parted_msa, seed));
    }

    /* generate starting trees for bootstrap searches */
    for (size_t b = 0; b < instance.opts.num_bootstraps; ++b)
    {
      auto tree = generate_tree(instance, StartingTree::random);

//      if (b < checkp.bs_trees.size())
//        continue;

      instance.bs_start_trees.emplace_back(move(tree));
    }
  }
}

void init_ancestral(RaxmlInstance& instance)
{
  if (instance.opts.command == Command::ancestral)
  {
    const auto& parted_msa = *instance.parted_msa;
    const Tree& tree = instance.start_trees.at(0);

    instance.ancestral_states = make_shared<AncestralStates>(tree.num_inner(), parted_msa);
  }
}

void reroot_tree_with_outgroup(const Options& opts, Tree& tree, bool add_root_node)
{
  if (!opts.outgroup_taxa.empty())
  {
    try
    {
      tree.reroot(opts.outgroup_taxa, add_root_node);
    }
    catch (std::runtime_error& e)
    {
      if (pll_errno == PLLMOD_TREE_ERROR_POLYPHYL_OUTGROUP)
        LOG_WARN << "WARNING: " << e.what() << endl << endl;
      else
        throw e;
    }
  }
}

void postprocess_tree(const Options& opts, Tree& tree)
{
  reroot_tree_with_outgroup(opts, tree, true);
  // TODO: collapse short branches
  // TODO: regraft previously removed duplicate seqs etc.
}

void draw_bootstrap_support(RaxmlInstance& instance, Tree& ref_tree,
                            const TreeTopologyList& bs_trees)
{
  reroot_tree_with_outgroup(instance.opts, ref_tree, false);

  for (auto metric: instance.opts.bs_metrics)
  {
      shared_ptr<SupportTree> sup_tree;
      bool support_in_pct = false;

      if (metric == BranchSupportMetric::fbp)
      {
        sup_tree = make_shared<BootstrapTree>(ref_tree);
        support_in_pct = true;
      }
      else if (metric == BranchSupportMetric::tbe)
      {
        sup_tree = make_shared<TransferBootstrapTree>(ref_tree, instance.opts.tbe_naive);
        support_in_pct = false;
      }
      else
        assert(0);

      Tree tree = ref_tree;
      for (auto bs: bs_trees)
      {
        tree.topology(bs);
        sup_tree->add_replicate_tree(tree);
      }
      sup_tree->draw_support(support_in_pct);

      instance.support_trees[metric] = sup_tree;
  }
}

bool check_bootstop(const RaxmlInstance& instance, const TreeTopologyList& bs_trees,
                    bool print = false)
{
  if (!instance.bootstop_checker)
    return false;

  const auto& opts       = instance.opts;
  auto& bootstop_checker = instance.bootstop_checker;

  if (!bootstop_checker->max_bs_trees())
    bootstop_checker->max_bs_trees(bs_trees.size());

  if (print)
  {
    LOG_INFO << "Performing bootstrap convergence assessment using autoMRE criterion"
             << endl << endl;

    // # Trees     Avg WRF in %    # Perms: wrf <= 2.00 %
    LOG_INFO << " # trees       "
             << " avg WRF      "
             << " avg WRF in %      "
             << " # perms: wrf <= " << setprecision(2) << opts.bootstop_cutoff * 100 << " %    "
             << " converged?  " << endl;
  }

  assert(!instance.random_tree.empty());

  Tree bs_tree = instance.random_tree;
  size_t bs_num = 0;
  bool converged = false;
  for (auto it: bs_trees)
  {
    bs_tree.topology(it);

    bootstop_checker->add_bootstrap_tree(bs_tree);

    bs_num++;

    if (bs_num % opts.bootstop_interval == 0 || bs_num == bs_trees.size())
    {
      converged = bootstop_checker->converged(rand());

      if (print)
      {
        LOG_INFO << setw(8) << bs_num << " "
                 << setw(14) << setprecision(3) << bootstop_checker->avg_wrf() << "   "
                 << setw(16) << setprecision(3) << bootstop_checker->avg_pct() << "   "
                 << setw(26) << bootstop_checker->num_better() << "        "
                 << (converged ? "YES" : "NO") << endl;
      }

      if (converged)
        break;
    }
  }

  if (print)
  {
    LOG_INFO << "Bootstopping test " << (converged ? "converged" : "did not converge")
             << " after " <<  bootstop_checker->num_bs_trees() << " trees" << endl << endl;
  }

  return converged;
}

TreeTopologyList read_newick_trees(Tree& ref_tree, const std::string& fname,
                                   const std::string& tree_kind)
{
  NameIdMap ref_tip_ids;
  TreeTopologyList trees;
  unsigned int bs_num = 0;

  if (!sysutil_file_exists(fname))
    throw runtime_error("File not found: " + fname);

  NewickStream boots(fname, std::ios::in);
  auto tree_kind_cap = tree_kind;
  tree_kind_cap[0] = toupper(tree_kind_cap[0]);

  LOG_INFO << "Reading " << tree_kind << " trees from file: " << fname << endl;

  while (boots.peek() != EOF)
  {
    Tree tree;
    boots >> tree;

    if (trees.empty())
    {
      if (ref_tree.empty())
        ref_tree = tree;
      ref_tip_ids = ref_tree.tip_ids();
    }

    assert(!ref_tip_ids.empty());

    if (!tree.binary())
    {
      LOG_DEBUG << "REF #branches: " << ref_tree.num_branches()
                << ", BS #branches: " << tree.num_branches() << endl;
      throw runtime_error(tree_kind_cap + " tree #" + to_string(bs_num+1) +
                          " contains multifurcations!");
    }

    try
    {
      tree.reset_tip_ids(ref_tip_ids);
    }
    catch (out_of_range& e)
    {
      throw runtime_error(tree_kind_cap + " tree #" + to_string(bs_num+1) +
                          " contains incompatible taxon name(s)!");
    }
    catch (invalid_argument& e)
    {
      throw runtime_error(tree_kind_cap + " tree #" + to_string(bs_num+1) +
                          " has wrong number of tips: " + to_string(tree.num_tips()));
    }
    trees.push_back(tree.topology());
    bs_num++;
  }

  LOG_INFO << "Loaded " << trees.size() << " trees with "
           << ref_tree.num_tips() << " taxa." << endl << endl;

  return trees;
}

TreeTopologyList read_bootstrap_trees(const RaxmlInstance& instance, Tree& ref_tree)
{
  auto bs_trees = read_newick_trees(ref_tree, instance.opts.bootstrap_trees_file(), "bootstrap");

  if (bs_trees.size() < 2)
  {
    throw runtime_error("You must provide a file with multiple bootstrap trees!");
  }

  return bs_trees;
}

void read_multiple_tree_files(RaxmlInstance& instance)
{
  const auto& opts = instance.opts;

  vector<string> fname_list;
  if (sysutil_file_exists(opts.tree_file))
    fname_list.push_back(opts.tree_file);
  else
    fname_list = split_string(opts.tree_file, ',');

  for (const auto& fname: fname_list)
  {
    Tree ref_tree = instance.random_tree;
    auto topos = read_newick_trees(ref_tree, fname, "input");
    for (const auto& t: topos)
    {
      ref_tree.topology(t);
      instance.start_trees.emplace_back(ref_tree);
    }
  }
}

void command_bootstop(RaxmlInstance& instance)
{
  auto bs_trees = read_bootstrap_trees(instance, instance.random_tree);

  check_bootstop(instance, bs_trees, true);
}

void command_support(RaxmlInstance& instance)
{
  const auto& opts = instance.opts;

  LOG_INFO << "Reading reference tree from file: " << opts.tree_file << endl;

  if (!sysutil_file_exists(opts.tree_file))
    throw runtime_error("File not found: " + opts.tree_file);

  Tree ref_tree;
  NewickStream refs(opts.tree_file, std::ios::in);
  refs >> ref_tree;

  LOG_INFO << "Reference tree size: " << to_string(ref_tree.num_tips()) << endl << endl;

  /* read all bootstrap trees from a Newick file */
  auto bs_trees = read_bootstrap_trees(instance, ref_tree);

  draw_bootstrap_support(instance, ref_tree, bs_trees);
  check_bootstop(instance, bs_trees, true);
}

void command_rfdist(RaxmlInstance& instance)
{
  const auto& opts = instance.opts;

  if (opts.start_trees.count(StartingTree::random) +
      opts.start_trees.count(StartingTree::parsimony) > 0)
  {
    /* generate random/parsimony trees -> we need an MSA for this */
    assert(!opts.msa_file.empty());
    load_parted_msa(instance);
    build_start_trees(instance, 0);
  }
  else
  {
    /* load trees from Newick file(s) */
    read_multiple_tree_files(instance);
  }

  if (instance.start_trees.size() < 2)
    throw runtime_error("Cannot compute RF distances since tree file contains fewer than 2 trees!");

  instance.dist_calculator.reset(new RFDistCalculator(instance.start_trees));
}

void command_consense(RaxmlInstance& instance)
{
  const auto& opts = instance.opts;

  /* load trees from Newick file(s) */
  read_multiple_tree_files(instance);

  if (instance.start_trees.size() < 2)
    throw runtime_error("Cannot consensus tree since tree file contains fewer than 2 trees!");

  instance.consens_tree.reset(new ConsensusTree(instance.start_trees, opts.consense_cutoff));
  if (instance.consens_tree)
    instance.consens_tree->draw_support();
  else
    runtime_error("Cannot create consensus tree!");
}

void command_bsmsa(RaxmlInstance& instance, const CheckpointFile& checkp)
{
  load_parted_msa(instance);
  generate_bootstraps(instance, checkp);
}

void check_terrace(const RaxmlInstance& instance, const Tree& tree)
{
#ifdef _RAXML_TERRAPHAST
  const auto& parted_msa = *instance.parted_msa;

  if (parted_msa.part_count() > 1 && instance.opts.brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    try
    {
      TerraceWrapper terrace_wrapper(parted_msa, tree);
      auto terrace_size = terrace_wrapper.terrace_size();
      if (terrace_size > 1)
      {
        LOG_WARN << "WARNING: Best-found ML tree lies on a terrace of size: "
                 << terrace_size << endl << endl;

        if (!instance.opts.terrace_file().empty())
        {
          ofstream fs(instance.opts.terrace_file());
          terrace_wrapper.print_terrace_compressed(fs);
          LOG_INFO << "Tree terrace (in compressed Newick format) was saved to: "
              << sysutil_realpath(instance.opts.terrace_file()) << endl;

          if (terrace_size <= instance.opts.terrace_maxsize)
          {
            auto nwk_fname = instance.opts.terrace_file() + "Newick";
            ofstream fsn(nwk_fname);
            terrace_wrapper.print_terrace_newick(fsn);
            LOG_INFO << "Tree terrace (in multi-line Newick format) was saved to: "
                << sysutil_realpath(nwk_fname) << endl;
          }

          LOG_INFO << endl;
        }
      }
      else
      {
        LOG_INFO << "NOTE: Tree does not lie on a phylogenetic terrace." << endl << endl;
      }
    }
    catch (terraces::no_usable_root_error& e)
    {
      if (instance.opts.command == Command::terrace)
      {
        LOG_ERROR << "ERROR: Cannot check for phylogenetic terraces "
            "since no comprehensive taxon is found." << endl << endl;
      }
      else
      {
        LOG_VERB << "NOTE: Cannot check for phylogenetic terraces "
            "since no comprehensive taxon is found." << endl << endl;
      }
    }
    catch (std::runtime_error& e)
    {
      LOG_ERROR << "ERROR: Unexpected terraphast error: " << e.what() << endl << endl;
    }
  }
#else
  RAXML_UNUSED(instance);
  RAXML_UNUSED(tree);
#endif
}

void save_ml_trees(const Options& opts, const CheckpointFile& checkp)
{
  NewickStream nw(opts.ml_trees_file(), std::ios::out);
  for (auto topol: checkp.ml_trees)
  {
    Tree ml_tree = checkp.tree();
    ml_tree.topology(topol.second.second);
    postprocess_tree(opts, ml_tree);
    nw << ml_tree;
  }
}

void print_ic_scores(const RaxmlInstance& instance, double loglh)
{
  const auto& parted_msa = *instance.parted_msa;

  size_t free_params = total_free_params(instance);
  size_t sample_size = parted_msa.total_sites();

  ICScoreCalculator ic_calc(free_params, sample_size);
  auto ic_scores = ic_calc.all(loglh);

  LOG_INFO << "AIC score: " << ic_scores[InformationCriterion::aic] << " / ";
  LOG_INFO << "AICc score: " << ic_scores[InformationCriterion::aicc] << " / ";
  LOG_INFO << "BIC score: " << ic_scores[InformationCriterion::bic] << endl;
  LOG_INFO << "Free parameters (model + branch lengths): " << free_params << endl << endl;

  if (free_params >= sample_size)
  {
    LOG_WARN << "WARNING: Number of free parameters (K=" << free_params << ") "
             << "is larger than alignment size (n=" << sample_size << ").\n"
             << "         This might lead to overfitting and compromise tree inference results!\n"
            << endl << endl;
  }
}

void print_final_output(const RaxmlInstance& instance, const CheckpointFile& checkp)
{
  auto const& opts = instance.opts;
  const auto& parted_msa = *instance.parted_msa;

  if (opts.command == Command::search || opts.command == Command::all ||
      opts.command == Command::evaluate || opts.command == Command::ancestral)
  {
    auto model_log_lvl = parted_msa.part_count() > 1 ? LogLevel::verbose : LogLevel::info;
    const auto& ml_models = instance.ml_tree.models;

    assert(ml_models.size() == parted_msa.part_count());

    RAXML_LOG(model_log_lvl) << "Optimized model parameters:" << endl;

    for (size_t p = 0; p < parted_msa.part_count(); ++p)
    {
      RAXML_LOG(model_log_lvl) << "\n   Partition " << p << ": " <<
          parted_msa.part_info(p).name().c_str() << endl;
      RAXML_LOG(model_log_lvl) << ml_models.at(p);
    }

    RAXML_LOG(model_log_lvl) << endl;
  }

  if (opts.command == Command::search || opts.command == Command::all ||
      opts.command == Command::evaluate)
  {
    auto best_loglh = instance.ml_tree.loglh;

    LOG_INFO << endl;
    LOG_RESULT << "Final LogLikelihood: " << FMT_LH(best_loglh) << endl;
    LOG_INFO << endl;

    print_ic_scores(instance, best_loglh);

    Tree best_tree = instance.ml_tree.tree;

    check_terrace(instance, best_tree);

    postprocess_tree(opts, best_tree);

//    pll_utree_show_ascii(&best_tree.pll_utree_root(),
//                         PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_BRANCH_LENGTH | PLL_UTREE_SHOW_CLV_INDEX );
//    printf("\n\n");

    if (!opts.best_tree_file().empty())
    {
      NewickStream nw_result(opts.best_tree_file());
      nw_result << best_tree;

      LOG_INFO << "Best ML tree saved to: " << sysutil_realpath(opts.best_tree_file()) << endl;
    }

    if (opts.brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED && !opts.partition_trees_file().empty())
    {
      NewickStream nw_result(opts.partition_trees_file());

      for (size_t p = 0; p < parted_msa.part_count(); ++p)
      {
        best_tree.apply_partition_brlens(p);
        nw_result << best_tree;
      }

      LOG_INFO << "Best per-partition ML trees saved to: " <<
          sysutil_realpath(opts.partition_trees_file()) << endl;
    }

    if (checkp.ml_trees.size() > 1 && !opts.ml_trees_file().empty())
    {
      save_ml_trees(opts, checkp);

      LOG_INFO << "All ML trees saved to: " << sysutil_realpath(opts.ml_trees_file()) << endl;
    }
  }

  if (opts.command == Command::all || opts.command == Command::support)
  {
    assert(!instance.support_trees.empty());

    for (const auto& it: instance.support_trees)
    {
      postprocess_tree(instance.opts, *it.second);

      auto sup_file = opts.support_tree_file(it.first);
      if (!sup_file.empty())
      {
        NewickStream nw(sup_file, std::ios::out);
        nw << *it.second;

        std::string metric_name = "";
        if (it.first == BranchSupportMetric::fbp)
          metric_name = "Felsenstein bootstrap (FBP)";
        else if (it.first == BranchSupportMetric::tbe)
          metric_name = "Transfer bootstrap (TBE)";

        LOG_INFO << "Best ML tree with " << metric_name << " support values saved to: " <<
            sysutil_realpath(sup_file) << endl;
      }
    }
  }

  if (opts.command == Command::consense)
  {
    assert(instance.consens_tree);

    auto cons_file = opts.cons_tree_file();
    if (!cons_file.empty())
    {
      NewickStream nw(cons_file, std::ios::out);
      nw.brlens(false);
      nw << *instance.consens_tree;

      LOG_INFO << opts.consense_type_name() << " consensus tree saved to: " <<
          sysutil_realpath(cons_file) << endl;
    }
  }

  if (opts.command == Command::search || opts.command == Command::all ||
      opts.command == Command::evaluate)
  {
    if (!opts.best_model_file().empty())
    {
      RaxmlPartitionStream model_stream(opts.best_model_file(), true);
      model_stream.print_model_params(true);
      model_stream << fixed << setprecision(logger().precision(LogElement::model));
      model_stream << parted_msa;

      LOG_INFO << "Optimized model saved to: " << sysutil_realpath(opts.best_model_file()) << endl;
    }
  }

  if (opts.command == Command::bootstrap || opts.command == Command::all)
  {
    // TODO now only master process writes the output, this will have to change with
    // coarse-grained parallelization scheme (parallel start trees/bootstraps)
    if (!opts.bootstrap_trees_file().empty())
    {
  //    NewickStream nw(opts.bootstrap_trees_file(), std::ios::out | std::ios::app);
      NewickStream nw(opts.bootstrap_trees_file(), std::ios::out);

      for (auto topol: checkp.bs_trees)
      {
        Tree bs_tree = checkp.tree();
        bs_tree.topology(topol.second.second);
        postprocess_tree(opts, bs_tree);
        nw << bs_tree;
      }

      LOG_INFO << "Bootstrap trees saved to: " << sysutil_realpath(opts.bootstrap_trees_file()) << endl;
    }
  }

  if (opts.command == Command::bsmsa)
  {
    if (!opts.bootstrap_msa_file(1).empty())
    {
      PartitionedMSAView bs_msa_view(instance.parted_msa);

      bool print_part_file = instance.parted_msa->part_count() > 1;

      size_t bsnum = 0;
      for (const auto& bsrep: instance.bs_reps)
      {
        bsnum++;
        PhylipStream ps(opts.bootstrap_msa_file(bsnum));

        bs_msa_view.site_weights(bsrep.site_weights);
        ps << bs_msa_view;
      }

      LOG_INFO << "Bootstrap replicate MSAs saved to: "
               << sysutil_realpath(opts.bootstrap_msa_file(1))
               << "  ... " << endl
               << "                                   "
               << sysutil_realpath(opts.bootstrap_msa_file(opts.num_bootstraps)) << endl;

      if (print_part_file)
      {
        RaxmlPartitionStream ps(opts.bootstrap_partition_file(), ios::out);

        ps << bs_msa_view;

        LOG_INFO << endl;
        LOG_INFO << "Partition file for (all) bootstrap replicate MSAs saved to: "
                 << sysutil_realpath(opts.bootstrap_partition_file()) << endl << endl;

        LOG_INFO << "IMPORTANT: You MUST use the aforementioned adjusted partitioned file" << endl
                 << "           when running tree searches on bootstrap replicate MSAs!" << endl;
      }
    }
  }

  if (opts.command == Command::rfdist)
  {
    assert(instance.dist_calculator);

    const auto& rfcalc = *instance.dist_calculator;

    LOG_RESULT << "Average absolute RF distance in this tree set: " << rfcalc.avg_rf() << endl;
    LOG_RESULT << "Average relative RF distance in this tree set: " << rfcalc.avg_rrf() << endl;
    LOG_RESULT << "Number of unique topologies in this tree set: "  << rfcalc.num_uniq_trees() << endl;

    if (!opts.rfdist_file().empty())
    {
      fstream fs(opts.rfdist_file(), ios::out);
      fs << rfcalc;

      LOG_INFO << "\nPairwise RF distances saved to: " << sysutil_realpath(opts.rfdist_file()) << endl;
    }
  }

  if (opts.command == Command::ancestral)
  {
    assert(instance.ancestral_states);

    if (!opts.asr_probs_file().empty())
    {
      AncestralProbStream as(opts.asr_probs_file());
      as.precision(logger().precision(LogElement::other));
      as << *instance.ancestral_states;

      LOG_INFO << "Marginal ancestral probabilities saved to: " << sysutil_realpath(opts.asr_probs_file()) << endl;
    }

    if (!opts.asr_states_file().empty())
    {
      AncestralStateStream as(opts.asr_states_file());
      as << *instance.ancestral_states;

      LOG_INFO << "Reconstructed ancestral sequences saved to: " << sysutil_realpath(opts.asr_states_file()) << endl;
    }

    if (!opts.asr_tree_file().empty())
    {
      NewickStream nw_result(opts.asr_tree_file());
      nw_result << *instance.ancestral_states;

      LOG_INFO << "Node-labeled tree saved to: " << sysutil_realpath(opts.asr_tree_file()) << endl;
    }
  }

  if (!opts.log_file().empty())
      LOG_INFO << "\nExecution log saved to: " << sysutil_realpath(opts.log_file()) << endl;

  LOG_INFO << "\nAnalysis started: " << global_timer().start_time();
  LOG_INFO << " / finished: " << global_timer().current_time() << std::endl;
  LOG_INFO << "\nElapsed time: " << FMT_PREC3(global_timer().elapsed_seconds()) << " seconds";
  if (checkp.elapsed_seconds > 0.)
  {
    LOG_INFO << " (this run) / ";
    LOG_INFO << FMT_PREC3(checkp.elapsed_seconds + global_timer().elapsed_seconds()) <<
        " seconds (total with restarts)";
  }

  LOG_INFO << endl << endl;
}

void print_resources(const RaxmlInstance& instance)
{
  StaticResourceEstimator resEstimator(*instance.parted_msa, instance.opts);
  auto res = resEstimator.estimate();

  LOG_VERB << "* Per-taxon CLV size (elements)                : "
      << res.taxon_clv_size << endl;
  LOG_INFO << "* Estimated memory requirements                : " <<
      (size_t) (((float) res.total_mem_size) / (1024 * 1024) + 1) << " MB" << endl << endl;
  LOG_INFO << "* Recommended number of threads / MPI processes: "
           << res.num_threads_balanced << endl;
  LOG_VERB << "* Maximum     number of threads / MPI processes: "
           << res.num_threads_response << endl;
  LOG_VERB << "* Minimum     number of threads / MPI processes: "
           << res.num_threads_throughput << endl;

  LOG_INFO << endl << "Please note that numbers given above are rough estimates only. " << endl <<
      "Actual memory consumption and parallel performance on your system may differ!"
      << endl << endl;
}

void init_parallel_buffers(const RaxmlInstance& instance)
{
  auto const& parted_msa = *instance.parted_msa;
  auto const& opts = instance.opts;

  // we need 2 doubles for each partition AND threads to perform parallel reduction,
  // so resize the buffer accordingly
  const size_t reduce_buffer_size = std::max(1024lu, 2 * sizeof(double) *
                                     parted_msa.part_count() * ParallelContext::num_threads());

  size_t worker_buf_size = 0;
  if (ParallelContext::num_ranks() > 1)
  {
    auto model_size = BinaryStream::serialized_size(parted_msa.models());
    auto tree_size = BinaryStream::serialized_size(instance.random_tree.topology());

    // buffer needs enough space to store serialized model parameters
    worker_buf_size = model_size;

    // for coarse-grained, add extra space to store ML/BS trees sent from workers to master
    if (ParallelContext::num_groups() > 1)
      worker_buf_size += (opts.bootstop_interval / ParallelContext::num_groups() + 1) * tree_size;

    // add some reserve
    worker_buf_size *= 1.2;
  }

  LOG_INFO << "Parallel reduction/worker buffer size: " << reduce_buffer_size/1024 <<  " KB  / "
            <<  worker_buf_size/1024 << " KB\n\n";

  ParallelContext::resize_buffers(reduce_buffer_size, worker_buf_size);
}

void load_assignment_data_for_this_rank(RaxmlInstance& instance) {
  // TODO: Fault-tolerance in case of coarse-grained parallelization
  // doesn't work with coarse-grained parallelization!
  assert(ParallelContext::num_groups() == 1);

  // collect PartitionAssignments from all worker threads
  PartitionAssignment local_part_ranges;
  assert(instance.opts.num_threads == 1);
  for (size_t i = 0; i < instance.opts.num_threads; ++i)
  {
    auto thread_ranges = instance.proc_part_assign.at(ParallelContext::local_proc_id() + i);
    for (auto& r: thread_ranges) {
      local_part_ranges.assign_sites(r.part_id, r.start, r.length);
    }
  }

  LOG_DEBUG << "Loading MSA segments from RBA file " << instance.opts.msa_file << "..." << endl;

  RBAStream bs(instance.opts.binary_msa_file());
  bs >> RBAStream::RBAOutput(*instance.parted_msa, RBAStream::RBAElement::seqdata, &local_part_ranges);
}

shared_ptr<vector<double>> compute_work_by_partition(RaxmlInstance& instance) {
  auto work_by_rank = ProfilerRegister::getInstance()->work_by_rank();
  assert(work_by_rank->size() == ParallelContext::num_ranks());

  auto work_by_partition = make_shared<vector<double>>(instance.parted_msa->part_count());

  assert(ParallelContext::num_procs() == instance.proc_part_assign.size());
  for (size_t proc = 0; proc < ParallelContext::num_procs(); proc++) {
    assert(ParallelContext::rank_id(ParallelContext::proc_id() == ParallelContext::rank_id()));
    size_t rank_id = ParallelContext::rank_id(proc);
    auto& assignment_list = instance.proc_part_assign.at(proc);

    // Determine the work per site for this rank
    assert(ParallelContext::num_procs() == ParallelContext::num_ranks());
    size_t nSites = 0;
    for (auto& assignment: assignment_list) {
      nSites += assignment.length;
    }
    assert(nSites > 0);
    assert(work_by_rank->at(rank_id) > 0);
    double work_per_site = work_by_rank->at(rank_id) / nSites;
    assert(work_per_site > 0);

    for (auto& assignment: assignment_list) {
      work_by_partition->at(assignment.part_id) += work_per_site * assignment.length;
    }
  }

  return work_by_partition;
}

bool firstLoad = true;
void thread_infer_ml(RaxmlInstance& instance, CheckpointManager& cm)
{
  auto& worker = instance.get_worker();
  Checkpoint& checkp = cm.checkpoint();
  auto const& master_msa = *instance.parted_msa;
  auto const& opts = instance.opts;

  unique_ptr<TreeInfo> treeinfo;

  auto gather_ml_trees = [&instance, &cm](unsigned int& batch_id) -> void
    {
      if (instance.opts.coarse() && ParallelContext::num_ranks() > 1)
      {
        ParallelContext::global_barrier();

        if (ParallelContext::group_master_thread())
          cm.gather_ml_trees();

        ParallelContext::global_thread_barrier();
        batch_id++;
      }
    };

  if (opts.command == Command::evaluate)
  {
    LOG_INFO << "\nEvaluating " << opts.num_searches <<
        " trees" << endl;
  }
  else
  {
    LOG_INFO << "\nStarting ML tree search with " << opts.num_searches <<
        " distinct starting trees" << endl;
  }

  (instance.start_trees.size() > 1 ? LOG_RESULT : LOG_INFO) << endl;

  unsigned int batch_id = (instance.done_ml_trees.size() / opts.bootstop_interval) + 1;

  auto compute_part_masters = [&instance]() {
    assert(ParallelContext::rank_id(ParallelContext::proc_id()) == ParallelContext::rank_id());
    
    if (instance.ranks_which_are_part_masters == nullptr) {
      instance.ranks_which_are_part_masters = make_shared<IDVector>();
    } else {
      instance.ranks_which_are_part_masters->clear();
    }
    instance.ranks_which_are_part_masters->reserve(ParallelContext::num_ranks());
    assert(instance.ranks_which_are_part_masters->empty());

    assert(ParallelContext::num_procs() == instance.proc_part_assign.size());
    for (size_t proc = 0; proc < ParallelContext::num_procs(); proc++) {
      auto& assignment = instance.proc_part_assign.at(proc);
      // assignment is the assignment to one PE (rank or thread)
      for (auto& range: assignment) {
        // Is the processor this assignment belongs to master for this range?
        size_t rank_id = ParallelContext::rank_id(proc);
        if (range.master()) {
          // Add each rank only once
          if (instance.ranks_which_are_part_masters->empty() ||
              instance.ranks_which_are_part_masters->back() != rank_id) {
            instance.ranks_which_are_part_masters->push_back(rank_id);
          }
          break;
        } 
      }
    }
    CheckpointManager::set_model_masters(instance.ranks_which_are_part_masters);
  };

  auto redo_assignment_cb = [&instance, &compute_part_masters] (bool rebalance = false) {
    auto profiler_register = ProfilerRegister::getInstance();

    LOG_DEBUG << "Rebalancing load:" << endl;
    profiler_register->profileFunction([&]() {
      if (rebalance) {
        LOG_DEBUG << "Using local work as weights to rebalance the load." << endl;
        auto work_by_partition = compute_work_by_partition(instance);
        balance_load(instance, bind(instance.load_balancer_cb, work_by_partition));
      } else {
        balance_load(instance, bind(instance.load_balancer_cb, nullptr));
      }
    }, "BalanceLoad");

    profiler_register->profileFunction([&]() {
      compute_part_masters();
    }, "ComputePartMasters");

    if (firstLoad) {
      profiler_register->profileFunction([&]() {
        load_assignment_data_for_this_rank(instance);
      }, "LoadAssignmentDataFirstLoad");
      firstLoad= false;
    } else {
      profiler_register->profileFunction([&]() {
        load_assignment_data_for_this_rank(instance);
      }, "LoadAssignmentData");
    }

    return instance.proc_part_assign.at(ParallelContext::local_proc_id());
  };

  auto ckp_tree_index = instance.run_phase == RaxmlRunPhase::mlsearch ? checkp.tree_index : 0;
  ParallelContext::thread_barrier();
  for (auto start_tree_num: worker.start_trees)
  {
    // Do this inside this loop as a node failure may have caused a reassignment
    auto const& part_assign = instance.proc_part_assign.at(ParallelContext::proc_id());
    compute_part_masters(); // TODO: Is this needed here?

    const auto& tree = instance.start_trees.at(start_tree_num-1);
    assert(!tree.empty());

    if (ckp_tree_index == start_tree_num)
    {
      // restore search state from checkpoint (tree + model params)
      treeinfo.reset(new TreeInfo(opts, checkp.tree, master_msa,
                                  instance.tip_msa_idmap, part_assign,
                                  redo_assignment_cb));
      assign_models(*treeinfo, checkp);
    }
    else
    {
      if (ParallelContext::group_master_thread())
        checkp.tree_index = start_tree_num;
      treeinfo.reset(new TreeInfo(opts, tree, master_msa, instance.tip_msa_idmap,
                                  part_assign, redo_assignment_cb));
    }

    treeinfo->set_topology_constraint(instance.constraint_tree);

    auto log_level = instance.start_trees.size() > 1 ? LogLevel::result : LogLevel::info;
    Optimizer optimizer(opts);
    if (opts.command == Command::evaluate || opts.command == Command::ancestral)
    {
      // check if we have anything to optimize
      if (opts.optimize_brlen || opts.optimize_model)
      {
        LOG_INFO_TS << "Tree #" << start_tree_num <<
            ", initial LogLikelihood: " << FMT_LH(treeinfo->loglh()) << endl;
        LOG_PROGR << endl;
        optimizer.evaluate(*treeinfo, cm);
      }
      else
      {
        double loglh = treeinfo->loglh();
        if (ParallelContext::master_thread())
          cm.search_state().loglh = loglh;

        cm.update_and_write(*treeinfo);
      }

      LOG_PROGR << endl;
      LOG_WORKER_TS(log_level) << "Tree #" << start_tree_num <<
                           ", final logLikelihood: " << FMT_LH(checkp.loglh()) << endl;
      LOG_PROGR << endl;
    }
    else
    {
      optimizer.optimize_topology(*treeinfo, cm);
      LOG_PROGR << endl;
      LOG_WORKER_TS(log_level) << "ML tree search #" << start_tree_num <<
                           ", logLikelihood: " << FMT_LH(checkp.loglh()) << endl;
      LOG_PROGR << endl;
    }

    cm.save_ml_tree();
    cm.reset_search_state();

    // coarse: collect ML trees from MPI workers
    if (start_tree_num > batch_id * opts.bootstop_interval)
      gather_ml_trees(batch_id);
  }

  gather_ml_trees(batch_id);

  // No failure mitigation as MPI parallelization is not supported when searching
  // for ancestral states.
  if (opts.command == Command::ancestral)
  {
    auto const& part_assign = instance.proc_part_assign.at(ParallelContext::local_proc_id());
    assert(!opts.use_pattern_compression);
    treeinfo->compute_ancestral(instance.ancestral_states, part_assign);
    ParallelContext::thread_barrier();
  }
}

void thread_infer_bootstrap(RaxmlInstance& instance, CheckpointManager& cm)
{
  auto const& opts = instance.opts;
  auto const& master_msa = *instance.parted_msa;
  auto& worker = instance.get_worker();
  Checkpoint& checkp = cm.checkpoint();

  unique_ptr<TreeInfo> treeinfo;

  auto gather_bs_trees = [&instance, &opts, &cm](unsigned int& batch_start, unsigned int& batch_end) -> void
    {
      ParallelContext::global_thread_barrier();

      if (ParallelContext::group_master_thread())
        cm.gather_bs_trees();

      /* check bootstrapping convergence */
      if (instance.bootstop_checker)
      {
        if (ParallelContext::master())
        {
          Tree tree = instance.random_tree;
          for (unsigned int  i = batch_start; i < batch_end; ++i)
          {
            tree.topology(cm.checkp_file().bs_trees.at(i+1).second);

            instance.bootstop_checker->add_bootstrap_tree(tree);
          }

          instance.bs_converged = instance.bootstop_checker->converged(opts.random_seed);

          if (instance.bs_converged)
          {
            auto num_bs_trees = cm.checkp_file().bs_trees.size();
            LOG_INFO_TS << "Bootstrapping converged after " << num_bs_trees << " replicates." << endl;
          }
        }

        if (ParallelContext::master_thread())
          ParallelContext::mpi_broadcast(&instance.bs_converged, sizeof(bool));
      }

      ParallelContext::global_thread_barrier();

      batch_start = batch_end;
      batch_end = std::min(opts.num_bootstraps, batch_end+opts.bootstop_interval);
    };

  if (!instance.bs_reps.empty())
  {
    if (opts.command == Command::all)
    {
      LOG_INFO << endl;
      LOG_INFO_TS << "ML tree search completed, best tree logLH: " <<
          FMT_LH(cm.checkp_file().ml_trees.best_score()) << endl << endl;
    }

    LOG_INFO_TS << "Starting bootstrapping analysis with " << opts.num_bootstraps
             << " replicates." << endl << endl;
  }

  /* infer bootstrap trees if needed */
  unsigned int bs_batch_start = instance.bootstop_checker ?
                                instance.bootstop_checker->num_bs_trees() : cm.checkp_file().bs_trees.size();

  ParallelContext::global_master_broadcast(&bs_batch_start, sizeof(unsigned int));

  auto bs_batch_offset = bs_batch_start % opts.bootstop_interval;
  unsigned int  bs_batch_end = bs_batch_start - bs_batch_offset + opts.bootstop_interval;
  auto bs_num = worker.bs_trees.cbegin();

  auto ckp_tree_index = instance.run_phase == RaxmlRunPhase::bootstrap ? checkp.tree_index : 0;

  ParallelContext::global_thread_barrier();

  while (!instance.bs_converged && bs_num != worker.bs_trees.cend())
  {
    const auto& bs_start_tree = instance.bs_start_trees.at(*bs_num - 1);
    auto bs_rep = instance.bs_reps.at(*bs_num - 1);

    // rebalance sites
    if (ParallelContext::group_master_thread())
    {
      worker.proc_part_assign = balance_load(instance, bs_rep.site_weights);
    }
    ParallelContext::thread_barrier();

    auto const& bs_part_assign = worker.proc_part_assign.at(ParallelContext::local_proc_id());

    if (ckp_tree_index == *bs_num)
    {
      // restore search state from checkpoint (tree + model params)
      treeinfo.reset(new TreeInfo(opts, checkp.tree, master_msa, instance.tip_msa_idmap,
                                  bs_part_assign, bs_rep.site_weights));
      assign_models(*treeinfo, checkp);
    }
    else
    {
      if (ParallelContext::group_master_thread())
        checkp.tree_index = *bs_num;
      treeinfo.reset(new TreeInfo(opts, bs_start_tree, master_msa, instance.tip_msa_idmap,
                                  bs_part_assign, bs_rep.site_weights));
    }

    treeinfo->set_topology_constraint(instance.constraint_tree);

    Optimizer optimizer(opts);
    optimizer.optimize_topology(*treeinfo, cm);

    LOG_PROGR << endl;
    LOG_WORKER_TS(LogLevel::info) << "Bootstrap tree #" << *bs_num <<
                                     ", logLikelihood: " << FMT_LH(checkp.loglh()) << endl;
    LOG_PROGR << endl;

    cm.save_bs_tree();
    cm.reset_search_state();

    bs_num++;

    if (bs_num == worker.bs_trees.cend() || *bs_num > bs_batch_end)
      gather_bs_trees(bs_batch_start, bs_batch_end);

    ParallelContext::thread_barrier();
  }

  /* special case: if this worker has no bsreps to infer, it still must synchronize! */
  if (worker.bs_trees.empty())
    gather_bs_trees(bs_batch_start, bs_batch_end);
}


void thread_main(RaxmlInstance& instance, CheckpointManager& cm)
{
  /* wait until master thread prepares all global data */
//  printf("WORKER: %u, LOCAL_THREAD: %u\n", ParallelContext::group_id(), ParallelContext::local_proc_id());
  ParallelContext::global_barrier();

  auto const& opts = instance.opts;

  if ((opts.command == Command::search || opts.command == Command::all ||
      opts.command == Command::evaluate || opts.command == Command::ancestral) &&
      !instance.start_trees.empty())
  {
    thread_infer_ml(instance, cm);
    ParallelContext::global_barrier();
  }

  if ((opts.command == Command::bootstrap || opts.command == Command::all))
  {
    thread_infer_bootstrap(instance, cm);
    ParallelContext::global_barrier();
  }

  (instance.start_trees.size() > 1 ? LOG_RESULT : LOG_INFO) << endl;
}

void master_main(RaxmlInstance& instance, CheckpointManager& cm)
{
  auto const& opts = instance.opts;

  /* init workers */
  assert(opts.num_workers > 0);
  for (size_t i = 0; i < ParallelContext::num_local_groups(); ++i)
  {
    const auto& grp = ParallelContext::thread_group(i);
    instance.workers.emplace_back(instance, grp.group_id);
  }

  /* if resuming from a checkpoint, use binary MSA (if exists) */
  if (!opts.redo_mode &&
      sysutil_file_exists(opts.checkp_file()) &&
      sysutil_file_exists(opts.binary_msa_file()) &&
      RBAStream::rba_file(opts.binary_msa_file(), true))
  {
    instance.opts.msa_file = opts.binary_msa_file();
    instance.opts.msa_format = FileFormat::binary;
  }

  load_parted_msa(instance);
  assert(instance.parted_msa);
  auto& parted_msa = *instance.parted_msa;

  load_constraint(instance);

  check_options(instance);

  instance.bs_converged = false;

  /* init template tree */
  srand(instance.opts.random_seed);
  instance.random_tree = generate_tree(instance, StartingTree::random);

  init_parallel_buffers(instance);

  /* load checkpoint */
  load_checkpoint(instance, cm);

  // Also initialize mini-checkpoint's datastructures
  CheckpointManager::init_models(instance.parted_msa->models());

  /* load/create starting tree if not already loaded from checkpoint */
  if (instance.start_trees.size() < opts.num_searches)
  {
    if (ParallelContext::master_rank() || !opts.constraint_tree_file.empty() ||
        opts.start_tree_file().empty())
    {
      /* only master MPI rank generates starting trees (doesn't work with constrained search) */
      build_start_trees(instance, 0);
      ParallelContext::global_mpi_barrier();
    }
    else
    {
      /* non-master ranks load starting trees from a file */
      ParallelContext::global_mpi_barrier();
      load_start_trees(instance);
    }
  }

  LOG_VERB << endl << "Initial model parameters:" << endl;
  for (size_t p = 0; p < parted_msa.part_count(); ++p)
  {
    LOG_VERB << "   Partition: " << parted_msa.part_info(p).name() << endl <<
        parted_msa.model(p) << endl;
  }

  /* run load balancing algorithm */
  instance.load_balancer_cb = build_load_balancer(instance);
  balance_load(instance, bind(instance.load_balancer_cb, nullptr));

  /* lazy-load part of the alignment assigned to the current MPI rank */
  if (opts.msa_format == FileFormat::binary && opts.use_rba_partload)
  {
    load_assignment_data_for_this_rank(instance);
  }

  // TEMP WORKAROUND: here we reset random seed once again to make sure that BS replicates
  // are not affected by the number of ML search starting trees that has been generated before
  srand(instance.opts.random_seed);

  /* generate bootstrap replicates */
  generate_bootstraps(instance, cm.checkp_file());

  balance_load_coarse(instance, cm.checkp_file());

  init_ancestral(instance);

  if (ParallelContext::master_rank())
    instance.opts.remove_result_files();

  thread_main(instance, cm);

  if (ParallelContext::master_rank())
  {
    instance.ml_tree = cm.checkp_file().best_tree();

    if (opts.command == Command::all)
    {
      auto& checkp = cm.checkp_file();

      TreeTopologyList bs_trees;
      for (auto& t: checkp.bs_trees)
        bs_trees.push_back(t.second.second);

      draw_bootstrap_support(instance, instance.ml_tree.tree, bs_trees);
    }

    const auto& ml_models = instance.ml_tree.models;
    assert(ml_models.size() == parted_msa.part_count());
    for (size_t p = 0; p < parted_msa.part_count(); ++p)
    {
      parted_msa.model(p, ml_models.at(p));
    }
  }
}

int clean_exit(int retval)
{
  // We have checked for rank failures after finishing our computations, we therefore can ignore
  // all futher rank failures, as they won't influence the result.
  try {
    ParallelContext::finalize(retval != EXIT_SUCCESS);
  } catch (ParallelContext::RankFailureException& e) {
    LOG_WARN << "A rank failure occurred when finalizing the MPI subsystem. This should not have influenced the computation." << endl;
  } catch (ParallelContext::UnrecoverableRankFailureException& e) {
    LOG_WARN << "A rank failure occurred when finalizing the MPI subsystem. This should not have influenced the computation." << endl;
  }
  return retval;
}

int internal_main(int argc, char** argv, void* comm)
{
  int retval = EXIT_SUCCESS;

  RaxmlInstance instance;
  auto& opts = instance.opts;

  ParallelContext::init_mpi(argc, argv, comm);

  opts.num_ranks = ParallelContext::num_ranks();

  logger().add_log_stream(&cout);

  CommandLineParser cmdline;
  try
  {
    cmdline.parse_options(argc, argv, opts);
  }
  catch (OptionException &e)
  {
    LOG_INFO << "ERROR: " << e.message() << std::endl;
    return clean_exit(EXIT_FAILURE);
  }

  /* handle trivial commands first */
  switch (opts.command)
  {
    case Command::help:
      print_banner();
      cmdline.print_help();
      return clean_exit(EXIT_SUCCESS);
      break;
    case Command::version:
      print_banner();
      return clean_exit(EXIT_SUCCESS);
      break;
    case Command::evaluate:
    case Command::search:
    case Command::bootstrap:
    case Command::all:
    case Command::support:
    case Command::start:
    case Command::terrace:
    case Command::bsmsa:
    case Command::rfdist:
    case Command::consense:
    case Command::ancestral:
      if (!opts.redo_mode && opts.result_files_exist())
      {
        LOG_ERROR << endl << "ERROR: Result files for the run with prefix `" <<
                            (opts.outfile_prefix.empty() ? opts.msa_file : opts.outfile_prefix) <<
                            "` already exist!\n" <<
                            "Please either choose a new prefix, remove old files, or add "
                            "--redo command line switch to overwrite them." << endl << endl;
        return clean_exit(EXIT_FAILURE);
      }
      break;
    case Command::bsconverge:
    default:
      break;
  }

  /* now get to the real stuff */
  try
  {
    // make sure all MPI ranks use the same random seed
    ParallelContext::mpi_broadcast(&opts.random_seed, sizeof(long));
    srand(opts.random_seed);

    logger().log_level(instance.opts.log_level);
    logger().precision(instance.opts.precision);

    /* only master process writes the log file */
    if (ParallelContext::master() && !instance.opts.log_file().empty())
    {
      auto mode = !instance.opts.redo_mode && sysutil_file_exists(instance.opts.checkp_file()) ?
          ios::app : ios::out;
      logger().set_log_filename(opts.log_file(), mode);
    }

    print_banner();
    LOG_INFO << opts;

    check_options_early(opts);

    if (opts.redo_mode)
    {
      LOG_WARN << "WARNING: Running in REDO mode: existing checkpoints are ignored, "
          "and all result files will be overwritten!" << endl << endl;
    }

    if (opts.force_mode)
    {
      LOG_WARN << "WARNING: Running in FORCE mode: "
               << (opts.safety_checks.isnone() ? "all" : "some")
               << " safety checks are disabled!"
               << endl << endl;
    }

    /* init bootstopping */
    switch (opts.bootstop_criterion)
    {
      case BootstopCriterion::autoMRE:
        instance.bootstop_checker.reset(new BootstopCheckMRE(opts.num_bootstraps,
                                                             opts.bootstop_cutoff,
                                                             opts.bootstop_permutations));
        break;
      case BootstopCriterion::none:
        break;
      default:
        throw runtime_error("Only autoMRE bootstopping criterion is supported for now, sorry!");
    }

    CheckpointManager cm(opts);

    switch (opts.command)
    {
      case Command::evaluate:
      case Command::search:
      case Command::bootstrap:
      case Command::all:
      case Command::ancestral:
      {
        /* init load balancer */
        switch(opts.load_balance_method)
        {
          case LoadBalancing::naive:
            instance.load_balancer.reset(new SimpleLoadBalancer());
            break;
          case LoadBalancing::kassian:
            instance.load_balancer.reset(new KassianLoadBalancer());
            break;
          case LoadBalancing::benoit:
            instance.load_balancer.reset(new BenoitLoadBalancer());
            break;
          default:
            assert(0);
        }

        // use naive coarse-grained load balancer for now
        instance.coarse_load_balancer.reset(new SimpleCoarseLoadBalancer());

        ParallelContext::init_pthreads(opts, std::bind(thread_main,
                                                       std::ref(instance),
                                                       std::ref(cm)));

        auto profiler_register = ProfilerRegister::createInstance("recovery_timings.csv");
        profiler_register->registerProfiler("recalculate-assignment");
        profiler_register->registerProfiler("reload-sites");
        profiler_register->registerProfiler("restore-models");
        profiler_register->registerProfiler("mini-checkpoints");

        while (42) {
          try {
            master_main(instance, cm);
            // Some rank failures might hide from some ranks until the next collective operation
            ParallelContext::check_for_rank_failure();
            break; // on success
          } catch (ParallelContext::RankFailureException& e) {
            LOG_ERROR << endl << "####################" << endl;
            LOG_PROGR << e.what() << " Restarting from checkpoint" << endl;
            LOG_PROGR << "####################" << endl << endl;
            instance.opts.redo_mode = false; // Prevent program from starting over after each rank failure
            srand(opts.random_seed); // Reset random number generator for enhanced reproducability
          } catch (ParallelContext::UnrecoverableRankFailureException& e) {
            LOG_ERROR << endl << "####################" << endl;
            LOG_ERROR << "MPI Communicator is invalid, cannot proceed." << endl;
            LOG_ERROR << "####################" << endl;
            break;
          }
        }
        ProfilerRegister::getInstance()->writeStats(ParallelContext::rankToProcessorName);
        break;
      }
      case Command::support:
        command_support(instance);
        break;
      case Command::bsconverge:
        command_bootstop(instance);
        break;
#ifdef _RAXML_TERRAPHAST
      case Command::terrace:
      {
        load_parted_msa(instance);
        assert(!opts.tree_file.empty());
        LOG_INFO << "Loading tree from: " << opts.tree_file << endl << endl;
        if (!sysutil_file_exists(opts.tree_file))
          throw runtime_error("File not found: " + opts.tree_file);
        instance.start_tree_stream.reset(new NewickStream(opts.tree_file, std::ios::in));
        Tree tree = generate_tree(instance, StartingTree::user);
        check_terrace(instance, tree);
        break;
      }
#endif
      case Command::check:
        opts.use_pattern_compression = false;
        /* fall through */
      case Command::parse:
      {
        load_parted_msa(instance);
        if (!opts.tree_file.empty())
        {
          LOG_INFO << "Loading tree from: " << opts.tree_file << endl << endl;
          if (!sysutil_file_exists(opts.tree_file))
            throw runtime_error("File not found: " + opts.tree_file);
          instance.start_tree_stream.reset(new NewickStream(opts.tree_file, std::ios::in));
          Tree tree = generate_tree(instance, StartingTree::user);
        }
        if (opts.command == Command::parse)
          print_resources(instance);

        LOG_INFO << "Alignment can be successfully read by RAxML-NG." << endl << endl;
        break;
      }
      case Command::start:
      {
        load_parted_msa(instance);
        build_start_trees(instance, 0);
        if (!opts.start_tree_file().empty())
        {
          LOG_INFO << "\nAll starting trees saved to: " <<
              sysutil_realpath(opts.start_tree_file()) << endl << endl;
        }
        else
        {
          LOG_INFO << "\nStarting trees have been successfully generated." << endl << endl;
        }
        break;
      }
      case Command::bsmsa:
      {
        command_bsmsa(instance, cm.checkp_file());
        break;
      }
      case Command::rfdist:
      {
        command_rfdist(instance);
        break;
      }
      case Command::consense:
      {
        command_consense(instance);
        break;
      }
      case Command::none:
      default:
        LOG_ERROR << "Unknown command!" << endl;
        retval = EXIT_FAILURE;
    }

    /* finalize */
    if (ParallelContext::master_rank())
      print_final_output(instance, cm.checkp_file());

    /* analysis finished successfully, remove checkpoint file */
    if (ParallelContext::group_master_rank())
      cm.remove();
  }
  catch(exception& e)
  {
    LOG_ERROR << endl << "ERROR: " << e.what() << endl << endl;
    retval = EXIT_FAILURE;
  }

  return clean_exit(retval);
}


#ifdef _RAXML_BUILD_AS_LIB

extern "C" int dll_main(int argc, char** argv, void* comm)
{
  return internal_main(argc, argv, comm);
}

#else

int main(int argc, char** argv)
{
  return internal_main(argc, argv, 0);
}

#endif
