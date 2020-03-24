#include <stdio.h>

#include "Checkpoint.hpp"
#include "io/binary_io.hpp"
#include "io/file_io.hpp"

using namespace std;

void Checkpoint::reset_search_state()
{
  search_state = SearchState();
};

// TODO: fix!
Tree CheckpointFile::tree() const
{
  return checkp_list.at(0).tree;
}

MLTree CheckpointFile::best_tree() const
{
  MLTree result;
  auto best = ml_trees.best()->second;
  result.loglh = best.first;
  result.tree = tree();
  result.tree.topology(best.second);
  result.models = best_models;
  return result;
}

CheckpointManager::CheckpointManager(const Options& opts) :
    _active(true), _ckp_fname(opts.checkp_file())
{
  _checkp_file.opts = opts;
}

const Checkpoint& CheckpointManager::checkpoint(size_t ckp_id) const
{
  return _checkp_file.checkp_list.at(ckp_id);
}

Checkpoint& CheckpointManager::checkpoint(size_t ckp_id)
{
  return _checkp_file.checkp_list.at(ckp_id);
}

void CheckpointManager::init_checkpoints(const Tree& tree, const ModelCRefMap& models)
{
  /* create one checkpoint per *local* worker */
  for (size_t i = 0; i < ParallelContext::num_local_groups(); ++i)
    _checkp_file.checkp_list.emplace_back();

  for (auto& ckp: _checkp_file.checkp_list)
  {
    ckp.tree = tree;
    for (auto it: models)
      ckp.models[it.first] = it.second;
  }

  for (auto it: models)
    _checkp_file.best_models[it.first] = it.second;
}

void CheckpointManager::write(const std::string& ckp_fname) const
{
  if (ParallelContext::master() ||
      (ParallelContext::group_master() && ParallelContext::num_groups() > 1))
  {
  //  printf("write ckp rank=%lu\n", ParallelContext::rank_id());

    backup();

    BinaryFileStream fs(ckp_fname, std::ios::out);

    fs << _checkp_file;

    remove_backup();
  }
}

bool CheckpointManager::read(const std::string& ckp_fname)
{
  if (sysutil_file_exists(ckp_fname))
  {
    try
    {
      BinaryFileStream fs(ckp_fname, std::ios::in);

      fs >> _checkp_file;

      return true;
    }
    catch (runtime_error& e)
    {
      reset_search_state();
      LOG_DEBUG << "Error reading checkpoint: " << e.what() << endl;
      return false;
    }
  }
  else
    return false;
}

void CheckpointManager::remove()
{
  if (sysutil_file_exists(_ckp_fname))
    std::remove(_ckp_fname.c_str());
}

void CheckpointManager::backup() const
{
  if (sysutil_file_exists(_ckp_fname))
    std::rename(_ckp_fname.c_str(), backup_fname().c_str());
}

void CheckpointManager::remove_backup() const
{
  if (sysutil_file_exists(backup_fname()))
    std::remove(backup_fname().c_str());
}

SearchState& CheckpointManager::search_state()
{
  if (_active)
    return checkpoint().search_state;
  else
  {
    return _empty_search_state;
  }
};

void CheckpointManager::reset_search_state()
{
  ParallelContext::thread_barrier();

  if (ParallelContext::group_master_thread())
  {
    if (_active)
      checkpoint().reset_search_state();
    else
      _empty_search_state = SearchState();
  }

  ParallelContext::thread_barrier();
};

void CheckpointManager::save_ml_tree()
{
  if (ParallelContext::group_master())
  {
    /* we will modify a global data in _checkp_file -> define critical section */
    ParallelContext::UniqueLock lock;

    Checkpoint& ckp = checkpoint();

//    printf("WORKER %u: save ML tree # %u index loglh = %lf\n",
//           ParallelContext::group_id(), index, ckp.loglh());

    auto& ml_trees = _checkp_file.ml_trees;
    if (ml_trees.empty() || ckp.loglh() > ml_trees.best_score())
      _checkp_file.best_models = ckp.models;

    ml_trees.insert(ckp.tree_index, ScoredTopology(ckp.loglh(), ckp.tree.topology()));

    if (_active)
      write();
  }
}

void CheckpointManager::save_bs_tree()
{
  if (ParallelContext::group_master())
  {
    /* we will modify a global data in _checkp_file -> define critical section */
    ParallelContext::UniqueLock lock;

    Checkpoint& ckp = checkpoint();

//    printf("WORKER %u: save BS tree # %u index loglh = %lf\n",
//           ParallelContext::group_id(), index, ckp.loglh());

    _checkp_file.bs_trees.insert(ckp.tree_index, ScoredTopology(ckp.loglh(), ckp.tree.topology()));

    if (_active)
      write();
  }
}

bool CheckpointManager::_models_initialized = false;
ModelMap CheckpointManager::_working_models;
ModelMap CheckpointManager::_tmp_models;

// The models have to be initalized once from instance.parted_msa.models(), as assing(model, treeinfo_partiton)
// won't copy over everything (e.g. rate heterogeneity).
void CheckpointManager::init_models(const ModelCRefMap& models) {
  // Do not assert that the models were not initialized before, as this may happen if we
  // restart from a checkpoint after rank failure.
  assert(models.size() > 0);
  for (auto m: models) {
    _working_models[m.first] = m.second;
    _tmp_models[m.first] = m.second;
  }
  _models_initialized = true;
}

shared_ptr<IDVector> CheckpointManager::_model_master_ranks = nullptr;
void CheckpointManager::set_model_masters(shared_ptr<IDVector> model_master_ranks) {
  assert(model_master_ranks->size() > 0);
  assert(model_master_ranks->size() <= all_models().size());
  assert(model_master_ranks->size() <= ParallelContext::num_ranks());
  _model_master_ranks = model_master_ranks;
}

const ModelMap& CheckpointManager::all_models() {
  assert(_models_initialized);
  assert(_working_models.size() > 0);
  return _working_models; 
}

void CheckpointManager::print_models() {
  for (auto& model: all_models()) {
    ParallelContext::log(model.second.to_string(true));
  }
}

string CheckpointManager::all_models_to_string() {
  string s = "";
  for (auto& model: all_models()) {
    s.append("<model " + to_string(model.first) + ">");
    s.append(model.second.to_string(true));
    s.append("brlen_scaler: " + to_string(model.second.brlen_scaler()));
  }
  return s;
}

// Only checks the working copy of the models, not the receiving copy
void CheckpointManager::assert_models_are_the_same_on_all_ranks() {
  #ifdef NON_FAILURE_TOLERANT_ASSERTS
  auto serialize = [](void * buf, size_t buf_size) -> int
  {
    BinaryStream bs((char*) buf, buf_size);
    bs << all_models().size();
    for (auto& model: all_models())
    {
      bs << model.first << model.second;
    }
    return (int) bs.pos();
  };

  auto deserialize = [](void * buf, size_t buf_size)
  {
    BinaryStream bs((char*) buf, buf_size);
    auto model_count = bs.get<size_t>();
    for (size_t m = 0; m < model_count; ++m)
    {
      size_t part_id;
      bs >> part_id;

      Model model = all_models().at(part_id);
      bs >> model;
      assert(model.to_string(true) == all_models().at(part_id).to_string(true));
    }
  };

  ParallelContext::global_broadcast_custom(serialize, deserialize, sizeof(Model) * all_models().size(), 0);
  #endif
}

// This method will pass on all RankFailureExceptions to the caller which should handle them.
// In case this method throws such an exception, it guarantees, that none of the models have
// been updated on any rank.
void CheckpointManager::update_models(const TreeInfo& treeinfo) {
  assert(ParallelContext::ranks_per_group() >= 1);
  if (ParallelContext::ranks_per_group() == 1) {
    for (auto p: treeinfo.parts_master()) {
      assign(_working_models[p], treeinfo, p);
    }
  } else {
    assert(_models_initialized);
    assert(_model_master_ranks != nullptr);
    assert(_working_models.size() == _tmp_models.size());
    IDSet modelsToSend;

    // Do not clear models. The assign does NOT update all aspects of the model. For example rate
    // heterogeneity will be lost. Models therefore need to be updated using assign.
    // _all_models.clear()

    // First, collect the model ids of those models that this rank is the "master" for, i.e. responsible
    // for sending the model parameters to the other ranks.
    assert(modelsToSend.empty());
    for (auto p: treeinfo.parts_master())
    {
      /* we will modify a global map -> define critical section */
      ParallelContext::GroupLock lock;

      assign(_tmp_models[p], treeinfo, p);
      modelsToSend.insert(p);
    }
    assert(modelsToSend.size() == treeinfo.parts_master().size());

    // Serialize all the models this rank has to send out into a bytestream which
    // will be broadcasted to all other ranks.
    auto serialize = [&modelsToSend](void * buf, size_t buf_size) -> int
      {
        assert(modelsToSend.size() > 0);
        BinaryStream bs((char*) buf, buf_size);
        bs << modelsToSend.size();
        for (size_t p: modelsToSend)
        {
          bs << p << _tmp_models.at(p);
          // ParallelContext::log("Sent model " + to_string(p) + ": " + _tmp_models.at(p).to_string(true));
        }
        return (int) bs.pos();
      };

    // Deserialize all models received from one rank.
    auto deserialize = [&modelsToSend](void * buf, size_t buf_size)
      {
        BinaryStream bs((char*) buf, buf_size);
        auto model_count = bs.get<size_t>();
        assert(model_count > 0);
        for (size_t m = 0; m < model_count; ++m)
        {
          size_t part_id;
          bs >> part_id;
          assert(part_id < _tmp_models.size());
          assert(modelsToSend.find(part_id) == modelsToSend.end());

          bs >> _tmp_models[part_id];
          // ParallelContext::log("Rcvd model " + to_string(part_id) + ": " + _tmp_models.at(part_id).to_string(true));
        }
      };

    // Now, perform one broadcast from every rank which is master of at least one model

    assert(_model_master_ranks->size() >= 1 &&
          _model_master_ranks->size() <= ParallelContext::num_ranks() &&
          _model_master_ranks->size() <= _tmp_models.size()
    );

    for (size_t rank: *_model_master_ranks) {
      // The sender needs to know large of a buffer to allocate.
      // The receivers will get a message from the sender describing the length of the encoding.
      size_t sizeOfBuffer;
      if (ParallelContext::rank_id() == rank) {
        assert(modelsToSend.size() > 0);
        sizeOfBuffer = sizeof(modelsToSend.size());
        for (size_t modelId: modelsToSend)  {
          size_t lengthOfModelEncoding = _tmp_models.at(modelId).encodingLength();
          assert(lengthOfModelEncoding > 0);
          sizeOfBuffer += sizeof(size_t) + lengthOfModelEncoding; // part_id + model
        }
        assert(sizeOfBuffer > 0);
      } else {
        sizeOfBuffer = 0;
      }
      try {
        ParallelContext::global_broadcast_custom(serialize, deserialize, sizeOfBuffer, rank);
      } catch (ParallelContext::RankFailureException& e){
        LOG_ERROR << "Rank failure during mini-checkpointing" << endl;
        throw e;
      }
    }

    // If no rank failrue have occurred, locally update the working copy of the models from the received models.
    // If a rank fails after this barrier, the received models will still be the state the search is restored to,
    // as every rank still alive has copied them over to their working copy.
    ParallelContext::check_for_rank_failure();
    assert(_tmp_models.size() == _working_models.size());
    assert(_working_models.size() > 0);
    for (size_t model_id = 0; model_id < _tmp_models.size(); model_id++) {
      _working_models[model_id] = _tmp_models[model_id];
      assert(_working_models[model_id].to_string(true) == _tmp_models[model_id].to_string(true));
      // ParallelContext::log("Copied model " + to_string(model_id) + ": " + _working_models[model_id].to_string(true));
    }

    assert_models_are_the_same_on_all_ranks();
  }
}

void CheckpointManager::update_and_write(const TreeInfo& treeinfo)
{
  ProfilerRegister::getInstance()->writeStats(ParallelContext::rankToProcessorName);

  Checkpoint& ckp = checkpoint();

  for (auto& model: _working_models)
  {
    /* we will modify a global map -> define critical section */
    ParallelContext::GroupLock lock;

    ckp.models.at(model.first) = model.second;
  }

  if (ParallelContext::group_master())
  {
    assign_tree(ckp, treeinfo);
    if (_active)
      write();
  }
}

void CheckpointManager::gather_ml_trees()
{
  if (ParallelContext::num_ranks() == 1 || ParallelContext::num_groups() == 1)
    return;

  /* send callback -> worker ranks */
  auto worker_cb = [this](void * buf, size_t buf_size) -> int
      {
        BinaryStream bs((char*) buf, buf_size);

        bs << _checkp_file.ml_trees;

        bs << _checkp_file.best_models;

//        printf("after worker: %u\n", bs.pos());

        // clear this batch of ML trees from the worker, since they will now be stored by master
        _checkp_file.ml_trees.clear();

        return (int) bs.pos();
      };

  /* receive callback -> master rank */
  auto master_cb = [this](void * buf, size_t buf_size)
     {
       double old_score = _checkp_file.ml_trees.best_score();
       BinaryStream bs((char*) buf, buf_size);

       bs >>  _checkp_file.ml_trees;

       if (_checkp_file.ml_trees.best_score() > old_score)
         bs >> _checkp_file.best_models;
     };

  ParallelContext::mpi_gather_custom(worker_cb, master_cb);
}

void CheckpointManager::gather_bs_trees()
{
  if (ParallelContext::num_ranks() == 1 || ParallelContext::num_groups() == 1)
    return;

  /* send callback -> worker ranks */
  auto worker_cb = [this](void * buf, size_t buf_size) -> int
      {
        BinaryStream bs((char*) buf, buf_size);

        bs << _checkp_file.bs_trees;

//        printf("after worker: %u\n", bs.pos());

        // clear this batch of BS trees from the worker, since they will now be stored by master
        _checkp_file.bs_trees.clear();

        return (int) bs.pos();
      };

  /* receive callback -> master rank */
  auto master_cb = [this](void * buf, size_t buf_size)
     {
       BinaryStream bs((char*) buf, buf_size);

       bs >>  _checkp_file.bs_trees;
     };

  ParallelContext::mpi_gather_custom(worker_cb, master_cb);
}

BasicBinaryStream& operator<<(BasicBinaryStream& stream, const Checkpoint& ckp)
{
  stream << ckp.search_state;

  stream << ckp.tree_index;

  stream << ckp.tree.topology();

  stream << ckp.models;

  return stream;
}

BasicBinaryStream& operator>>(BasicBinaryStream& stream, Checkpoint& ckp)
{
  stream >> ckp.search_state;

  stream >> ckp.tree_index;

  ckp.tree.topology(stream.get<TreeTopology>());

  stream >> ckp.models;

  return stream;
}

BasicBinaryStream& operator<<(BasicBinaryStream& stream, const CheckpointFile& ckpfile)
{
  stream << ckpfile.version;

  // NB: accumulated runtime from past runs + current elapsed time
  stream << ckpfile.elapsed_seconds + global_timer().elapsed_seconds();

  stream << ckpfile.opts;

  stream << ckpfile.checkp_list;

  stream << ckpfile.best_models;

  stream << ckpfile.ml_trees;

  stream << ckpfile.bs_trees;

  return stream;
}

BasicBinaryStream& operator>>(BasicBinaryStream& stream, CheckpointFile& ckpfile)
{
  stream >> ckpfile.version;

  if (ckpfile.version < RAXML_CKP_MIN_SUPPORTED_VERSION || ckpfile.version > RAXML_CKP_VERSION)
  {
    throw runtime_error("Unsupported checkpoint file version!");
  }

  stream >> ckpfile.elapsed_seconds;

  stream >> ckpfile.opts;

  {
    // we should take special care in case number of workers has been changed after restart:
    // - if #workers increased, "extra" workers will start tree search from scratch
    // - if #workers decreased, "extra" checkpoints in file will be ignored (not optimal, but simpler)
    size_t num_ckp_in_file = stream.get<size_t>();
    size_t num_ckp_to_load = std::min(num_ckp_in_file, ckpfile.checkp_list.size());
    auto dummy_ckp = (num_ckp_to_load < num_ckp_in_file) ? ckpfile.checkp_list[0] : Checkpoint();
    for (size_t i = 0; i < num_ckp_in_file; ++i)
    {
      if (i < num_ckp_to_load)
        stream >> ckpfile.checkp_list[i];
      else
        stream >> dummy_ckp;
    }
  }

  stream >> ckpfile.best_models;

  stream >> ckpfile.ml_trees;

  stream >> ckpfile.bs_trees;

  return stream;
}


void assign_tree(Checkpoint& ckp, const TreeInfo& treeinfo)
{
  ckp.tree = treeinfo.tree();
}

void assign_model(Checkpoint& ckp, const TreeInfo& treeinfo, size_t index)
{
  assign(ckp.models.at(index), treeinfo, index);
}

void assign_models(Checkpoint& ckp, const TreeInfo& treeinfo)
{
  for (auto p: treeinfo.parts_master())
    assign_model(ckp, treeinfo, p);
}

void assign_models(TreeInfo& treeinfo, const Checkpoint& ckp)
{
  const pllmod_treeinfo_t& pll_treeinfo = treeinfo.pll_treeinfo();
  for (auto& m: ckp.models)
  {
    if (!pll_treeinfo.partitions[m.first])
      continue;

    treeinfo.model(m.first, m.second);
  }
}

void assign(Checkpoint& ckp, const TreeInfo& treeinfo)
{
  assign_tree(ckp, treeinfo);
  assign_models(ckp, treeinfo);
}

void assign(TreeInfo& treeinfo, const Checkpoint& ckp)
{
  // TODO: it is currently not possible to change tree after pll_treeinfo has been created
  // this should be fixed while doing refactoring to change pll_unode_t -> pll_tree_t
  assert(0);
  treeinfo.tree(ckp.tree);

  assign_models(treeinfo, ckp);
}

