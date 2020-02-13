#ifndef RAXML_TREEINFO_HPP_
#define RAXML_TREEINFO_HPP_

#include "common.h"
#include "Tree.hpp"
#include "Options.hpp"
#include "AncestralStates.hpp"
#include "loadbalance/PartitionAssignment.hpp"
#include <memory>
#include <string>
#include <functional>

struct spr_round_params
{
  bool thorough;
  int radius_min;
  int radius_max;
  int ntopol_keep;
  double subtree_cutoff;
  cutoff_info_t cutoff_info;

  void reset_cutoff_info(double loglh)
  {
    cutoff_info.lh_dec_count = 0;
    cutoff_info.lh_dec_sum = 0.;
    cutoff_info.lh_cutoff = loglh / -1000.0;
  }
};

class TreeInfo
{
public:
  TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
            std::function<PartitionAssignment()> redo_assignment_cb);
  TreeInfo (const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
            const std::vector<uintVector>& site_weights);
  virtual
  ~TreeInfo ();

  const pllmod_treeinfo_t& pll_treeinfo() const { return *_pll_treeinfo; }
  const pll_unode_t& pll_utree_root() const { assert(_pll_treeinfo); return *_pll_treeinfo->root; }

  Tree tree() const;
  Tree tree(size_t partition_id) const;
  void tree(const Tree& tree);

  /* in parallel mode, partition can be share among multiple threads and TreeInfo objects;
   * this method returns list of partition IDs for which this thread is designated as "master"
   * and thus responsible for e.g. sending model parameters to the main thread. */
  const IDSet& parts_master() const { return _parts_master; }

  void model(size_t partition_id, const Model& model);

  void set_topology_constraint(const Tree& cons_tree);

  double loglh(bool incremental = false);
  double optimize_params(int params_to_optimize, double lh_epsilon);
  double optimize_params_all(double lh_epsilon)
  { return optimize_params(PLLMOD_OPT_PARAM_ALL, lh_epsilon); } ;
  double optimize_model(double lh_epsilon)
  { return optimize_params(PLLMOD_OPT_PARAM_ALL & ~PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE, lh_epsilon); } ;
  double optimize_branches(double lh_epsilon, double brlen_smooth_factor);
  double spr_round(spr_round_params& params);
  void compute_ancestral(const AncestralStatesSharedPtr& ancestral,
                         const PartitionAssignment& part_assign);

private:
  pllmod_treeinfo_t * _pll_treeinfo;
  IDSet _parts_master;
  int _brlen_opt_method;
  double _brlen_min;
  double _brlen_max;
  bool _check_lh_impr;
  doubleVector _partition_contributions;
  struct __partition_reinit_info_t {
    __partition_reinit_info_t(const Options &opts,
                              const PartitionedMSA& parted_msa,
                              const IDVector& tip_msa_idmap,
                              const std::vector<uintVector>& site_weights) :
      parted_msa(parted_msa), site_weights(site_weights), opts(opts), tip_msa_idmap(tip_msa_idmap) {}

    const PartitionedMSA& parted_msa;
    // This has to be copied, as the object passed to the constructor won't life long enough to be
    // around for the reinitialization.
    const std::vector<uintVector> site_weights;
    const Options &opts;
    const IDVector& tip_msa_idmap;
  };
  typedef __partition_reinit_info_t partition_reinit_info_t;
  std::shared_ptr<partition_reinit_info_t> partition_reinit_info = nullptr;

  void init(const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
            const IDVector& tip_msa_idmap, const PartitionAssignment& part_assign,
            const std::vector<uintVector>& site_weights);
  
  // Will use the stored partition info to reinitialize partitions
  void reinit_partitions(const PartitionAssignment& part_assign); 

  void init_partitions(const Options &opts, const Tree& tree, const PartitionedMSA& parted_msa,
                               const IDVector& tip_msa_idmap,
                               const PartitionAssignment& part_assign,
                               const std::vector<uintVector>& site_weights);

  void assert_lh_improvement(double old_lh, double new_lh, const std::string& where = "");

  // Mini checkpoints are used to restore to in case of rank failure. This creates one.
  void mini_checkpoint();

  // The callback function which can be used by this object to recalculate the sites to rank assignment
  std::function<PartitionAssignment()> redo_assignment_cb;

  // Calculate a new site to rank assignment and update the treeinfo structure to it. 
  void update_to_new_assignment();

  // Wrapper to call an optimization function until a pass succeeds without a rank failure.
  // On failure, the ParallelContext is updated to include only non-failed ranks, the models and tree stored at the
  // last mini_checkpoint() are restored and the optimization routine is invoked again.
  double fault_tolerant_optimization(std::string parameter, const std::function<double()> optimizer);
};

void assign(PartitionedMSA& parted_msa, const TreeInfo& treeinfo);
void assign(Model& model, const TreeInfo& treeinfo, size_t partition_id);


pll_partition_t* create_pll_partition(const Options& opts, const PartitionInfo& pinfo,
                                      const IDVector& tip_msa_idmap,
                                      const PartitionRange& part_region, const uintVector& weights);

#endif /* RAXML_TREEINFO_HPP_ */
