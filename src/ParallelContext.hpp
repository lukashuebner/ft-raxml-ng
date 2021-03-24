#ifndef RAXML_PARALLELCONTEXT_HPP_
#define RAXML_PARALLELCONTEXT_HPP_

#include <vector>
#include <set>
#include <unordered_map>
#include <memory>
#include <string>
#include <mutex>
#include <cassert>

#include <functional>

#include "io/BinaryStream.hpp"
#include "common.h"
#include "log.hpp"

#ifdef _RAXML_MPI
#include <mpi.h>
#endif

#ifdef _RAXML_PTHREADS
#include <thread>
#include <mutex>
typedef std::thread ThreadType;
typedef std::thread::id ThreadIDType;
typedef std::mutex MutexType;
typedef std::unique_lock<std::mutex> LockType;
#else
typedef int ThreadType;
typedef size_t ThreadIDType;
typedef int MutexType;
typedef int LockType;
#endif

class Options;

struct ThreadGroup
{
  size_t group_id;           /* global ID */
  size_t local_group_id;     /* ID within MPI rank */
  size_t num_threads;
  std::vector<char> reduction_buf;

  MutexType mtx;
  volatile unsigned int barrier_counter;
  volatile int proceed;

  ThreadGroup(size_t id, size_t local_id, size_t size, size_t bufsize = 0) :
    group_id(id), local_group_id(local_id), num_threads(size), reduction_buf(bufsize),
    mtx(), barrier_counter(0), proceed(0) {}

  ThreadGroup(ThreadGroup&& other):
    group_id(other.group_id), local_group_id(other.local_group_id),
    num_threads(other.num_threads), reduction_buf(std::move(other.reduction_buf)),
    mtx(), barrier_counter(other.barrier_counter), proceed(other.proceed) {}
};

class ParallelContext
{
public:
  static void init_mpi(int argc, char * argv[], void * comm);
  static void init_pthreads(const Options& opts, const std::function<void()>& thread_main);
  static void resize_buffers(size_t reduce_buf_size, size_t worker_buf_size = 0);

  static void saveProfilingData();
  static std::string rankToProcessorName(size_t rank);
  static void finalize(bool force = false);

  static size_t num_procs() { return _num_ranks * _num_threads; }
  static size_t num_threads() { return _num_threads; }
  static size_t num_ranks() { return _num_ranks; }
  static size_t num_nodes() { return _num_nodes; }
  static size_t num_groups() { return _num_groups; }
  static size_t num_local_groups() { return _thread_groups.size(); }
  static size_t ranks_per_node() { return _num_ranks / _num_nodes; }
  static size_t threads_per_group() { return num_procs() / _num_groups; }
  static size_t ranks_per_group() { return _num_ranks / _num_groups; }

  static void parallel_reduce_cb(void * context, double * data, size_t size, int op);
  static void thread_reduce(double * data, size_t size, int op);
  static void thread_broadcast(size_t source_id, void * data, size_t size);
  void thread_send_master(size_t source_id, void * data, size_t size) const;

  static void mpi_broadcast(void * data, size_t size, int root = 0);
  template<typename T> static void mpi_broadcast(T& obj)
  {
    if (_num_ranks > 1)
    {
      size_t size = master() ?
          BinaryStream::serialize(_parallel_buf.data(), _parallel_buf.capacity(), obj) : 0;
      mpi_broadcast((void *) &size, sizeof(size_t));
      mpi_broadcast((void *) _parallel_buf.data(), size);
      if (!master())
      {
        BinaryStream bs(_parallel_buf.data(), size);
        bs >> obj;
      }
    }
  }

  static void log(std::string message);

  static void global_master_broadcast(void * data, size_t size, int root = 0);
  static void global_master_broadcast_custom(std::function<int(void*, int)> prepare_send_cb,
                                             std::function<void(void*,int)> process_recv_cb,
                                             size_t sizfeOfBuffer);

  static void global_broadcast_custom(std::function<int(void*, int)> prepare_send_cb,
                                      std::function<void(void*,int)> process_recv_cb,
                                      size_t sizeOfBuffer, size_t root);

  // If a rank failure occurs during this call, detect_num_nodes() will be called in the failure
  // mitigation routine which itself will call mpi_gather_custom. In this case, we cannot use
  // parallel buf, as we would override its contents.
  static void mpi_gather_custom(std::function<int(void*,int)> prepare_send_cb,
                                std::function<void(void*,int)> process_recv_cb,
                                bool use_parallel_buf = true);

  static std::shared_ptr<std::vector<double>> mpi_allgather(double value);

  static void parallel_reduce(double * data, size_t size, int op);

  #ifdef _RAXML_MPI
  // Simulate failures via MPI_Comm_split instead of signaling ranks with SIGKILL
  #define RAXML_FAILURES_SIMULATE
  // Will throw a RankFailureException and repair the communicator if a rank failed
  static void check_for_rank_failure();
  // For testing purposes; the rank with id <rank> will fail when the fail function is called the n-th time on this rank
  static void fail(size_t rank, int on_nth_call = 1, bool reset = false);
  static void set_failure_prob(float probability);
  #endif

  // Will sleep until a debugger attaches and changes a local variable using for example
  // (gdb) set var i = 1
  static void waitForDebugger();

  static bool master() { return proc_id() == 0; }
  static bool master_rank() { return _rank_id == 0; }
  static bool master_thread() { return _thread_id == 0; }
  static size_t thread_id() { return _thread_id; }
  static size_t rank_id() { return _rank_id; }
  static size_t rank_id(size_t proc) { return proc / _num_threads; }
  static size_t proc_id() { return _rank_id * _num_threads + _thread_id; }
  static size_t group_id() { return _thread_group->group_id; }

  static size_t local_thread_id() { return _local_thread_id; }
  static size_t local_rank_id() { return _local_rank_id; }
  static size_t local_proc_id() { return _local_rank_id * _num_threads + _local_thread_id; }
  static size_t local_group_id() { return _thread_group->local_group_id; }

  static bool group_master() { return local_proc_id() == 0; }
  static bool group_master_rank() { return _local_rank_id == 0; }
  static bool group_master_thread() { return _local_thread_id == 0; }

  static ThreadGroup& thread_group(size_t id);

  static void barrier();
  static void global_barrier();
  static void thread_barrier();
  static void global_thread_barrier();
  static void mpi_barrier();
  static void global_mpi_barrier();

  /* static singleton, no instantiation/copying/moving */
  ParallelContext() = delete;
  ParallelContext(const ParallelContext& other) = delete;
  ParallelContext(ParallelContext&& other) = delete;
  ParallelContext& operator=(const ParallelContext& other) = delete;
  ParallelContext& operator=(ParallelContext&& other) = delete;

  class UniqueLock
  {
  public:
    UniqueLock() : _lock(mtx) {}
    UniqueLock(std::defer_lock_t tag) : _lock(mtx, tag) {}
    void lock() {_lock.lock(); }
  private:
    LockType _lock;
  };

  class GroupLock
  {
  public:
    GroupLock() : _lock(_thread_group->mtx) {}
  private:
    LockType _lock;
  };

  struct RankFailureException : public std::exception
  {
    const char * what() const throw() {
      return "MPI rank failure";
    }
  };

  struct UnrecoverableRankFailureException : public std::exception
  {
    const char * what() const throw() {
      return "Unrecoverable MPI rank failure";
    }
  };

  static long fail_every_nth_call;
  static int max_failures_to_simulate;
private:
#ifdef _RAXML_MPI
  // Has the MPI subsystem been finalized?
  static bool mpi_finalized();
  // Converts a MPI error code and respective error class into human readable format
  static std::string mpi_err_to_string(int errorCode);
  static int failureCounter;
  static bool _simulate_failure;
  static int iFail_at_call;
  static int simulated_failures;
  static float _randomized_failure_prob;
  
  template<class F>
  static void fault_tolerant_mpi_call(F mpi_call)
  {
    assert(_comm != MPI_COMM_NULL);
  
    #ifdef RAXML_FAILURES_SIMULATE
  
    if (_simulate_failure) {
      size_t oldRankId = rank_id();
      MPI_Comm newComm;
      RAXML_UNUSED(mpi_call);
  
      LOG_DEBUG << "Simulating failure" << std::endl;
  
      _simulate_failure = false; // Do this *before* calling detect_num_nodes, which uses ft-MPI calls 
      simulated_failures++;
      MPI_Comm_split(_comm, 0, (rank_id() + 1) % num_ranks(), &newComm);
      MPI_Comm_free(&_comm);
      _comm = newComm;
  
      update_world_parameters();
      assert(num_ranks() == 1 || oldRankId != rank_id());
      assert(num_ranks() != 0 || oldRankId == rank_id());
  
      // log("Rank failure at:");
      // print_stacktrace();
      throw RankFailureException();
    } else {
      mpi_call();
    }
   
    #else
  
    #ifdef RAXML_FAILURES_SIMULATE
    failureCounter++;
    if (failureCounter == iFail_at_call) { // Simulate my failure when my time has come
      raise(SIGKILL);
    }
    #endif
    
    int rc, ec;
    rc = mpi_call();
    MPI_Error_class(rc, &ec);
  
    if (ec == MPI_ERR_PROC_FAILED || ec == MPI_ERR_REVOKED) {
      if (ec == MPI_ERR_PROC_FAILED) {
          MPIX_Comm_revoke(_comm);
      }
  
      // Build a new communicator without the failed ranks
      MPI_Comm newComm;
      if ((rc = MPIX_Comm_shrink(_comm, &newComm)) != MPI_SUCCESS) {
        LOG_ERROR << "A rank failure was detected, but building the new communicator using "
                  << "MPI_Comm_shrink failed with err.: " << mpi_err_to_string(rc) << std::endl; 
        throw UnrecoverableRankFailureException();
      }
      assert(_comm != MPI_COMM_NULL);
      // As for the example in the ULFM documentation, freeing the communicator is necessary.
      // If trying to do so, however, MPI_ERR_COMM (invalid communicator) will be returned. 
      // --mca mpi_show_handle_leaks 1 does not show a leaked handle, so mayble MPIX_Comm_Shrink
      // already frees the communicator.
      // if ((rc = MPI_Comm_free(&_comm)) != MPI_SUCCESS) {
      //   LOG_ERROR << "A rank failure was detected, but freeing the old communicator using "
      //             << "MPI_Comm_free failed with err.: " << mpi_err_to_string(rc) << endl;
      //   throw UnrecoverableRankFailureException();
      // }
      _comm = newComm;
      update_world_parameters();
     throw RankFailureException();
    } else if (rc != MPI_SUCCESS) {
      throw runtime_error("MPI call did non fail because of a faulty rank but still did not return MPI_SUCCESS");
    }
  
    #endif
  }
#endif
  static void update_world_parameters();
  static std::vector<ThreadType> _threads;
  static size_t _num_threads;
  static size_t _num_ranks;
  static size_t _num_nodes;
  static size_t _num_groups;
  static std::vector<char> _parallel_buf;
  static std::unordered_map<ThreadIDType, ParallelContext> _thread_ctx_map;
  static MutexType mtx;

  static std::vector<std::string> _rankToProcessorName;

  static size_t _rank_id;
  static size_t _local_rank_id;
  static thread_local size_t _thread_id;
  static thread_local size_t _local_thread_id;
  static thread_local ThreadGroup * _thread_group;

  static std::vector<ThreadGroup> _thread_groups;

#ifdef _RAXML_MPI
  static bool _owns_comm;
  static MPI_Comm _comm;
#endif

  static void start_thread(size_t thread_id, size_t local_thread_id,
                           ThreadGroup& thread_grp,
                           const std::function<void()>& thread_main);
  static void detect_num_nodes();
};

// Some assertions communicate over the network (for example to compute the loglh)
// They may not handle failures correctly and should be disabled when simulating
// failures.
#ifdef RAXML_FAILURES_SIMULATE
//#define NON_FAILURE_TOLERANT_ASSERTS
#endif

#endif /* RAXML_PARALLELCONTEXT_HPP_ */
