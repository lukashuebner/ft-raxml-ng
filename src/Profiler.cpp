#include <memory>
#include <chrono>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <fstream>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <sstream>

#include "Profiler.hpp"
#include "io/binary_io.hpp"
//#include "io/file_io.hpp"

using namespace std;

FractionalProfiler::FractionalHistogram::FractionalHistogram() {
    bins = std::make_shared<vector<uint64_t>>(numBins, 0);
}

void FractionalProfiler::FractionalHistogram::event(float frac) {
    assert(bins != nullptr);

    int16_t bin = getBin(frac);
    assert(bin >= 0);
    assert(bin < numBins);
    assert(numBins == bins->size());

    (*bins)[bin]++;
}

uint64_t& FractionalProfiler::FractionalHistogram::operator[](size_t idx) {
    return (*bins)[idx];
} 

const uint64_t& FractionalProfiler::FractionalHistogram::operator[](size_t idx) const {
    return (*bins)[idx];
}

int16_t FractionalProfiler::FractionalHistogram::getBin(const float frac) {
    // Memory layout
    // bins [0..99]:    1/frac \in (1.001..1.1), 1/frac \in [1.1..1.2), ..., 1/frac \in [10.0..inf)
    // bin  100:        0.999..1.001
    // bins [101..200]: frac \in (1.001..1.1), frac \in [1.1..1.2), ..., frac \in [10.0..inf)
    assert(frac > 0);
    int16_t bin;
    if (abs(frac - 1) < 0.001) {
        bin = 100;
    } else if (frac > 1) {
        if (frac >= 10) {
            bin = 200;
        } else {
            assert(frac > 1 && frac < 10);
            bin = (int16_t) ((frac - 1) * 10) + 101;
        }
    } else { // frac < 1
        const float oneOverFrac = 1.0 / frac;
        assert(oneOverFrac > 1);
        if (oneOverFrac >= 10) {
            bin = 99;
        } else {
            assert(oneOverFrac < 10);
            bin = (int16_t) ((oneOverFrac - 1) * 10); 
        }
    }
    assert(bin >= 0 && bin < FractionalHistogram::numBins);
    assert(frac > 1 ? bin >= 100 : bin <= 100);
    return bin;
}

FractionalProfiler::FractionalProfiler(string name) : name(name) {}

void FractionalProfiler::startTimer(bool resume) {
    if (!resume && running) {
        invalid = true;
        throw runtime_error("The timer is already running.");
    } else if (invalid) {
        throw runtime_error("The timer is in an invalid state, probably the call sequence was not (start end)^n");
    } else if (!resume && !savedOrDiscarded) {
        invalid = true;
        throw runtime_error("Timer has not been saved yet! Either save or discard it.");
    } else if (resume && savedOrDiscarded) {
        invalid = true;
        throw runtime_error("Can't resume a timer that has already been saved or discarded.");
    }
    
    savedOrDiscarded = false;
    running = true;
    start = chrono::high_resolution_clock::now();
    assert(!(resume && firstStart == chrono::time_point<chrono::high_resolution_clock>::max()));
    if (!resume) {
        nsPassed = 0;
        if (firstStart == chrono::time_point<chrono::high_resolution_clock>::max()) {
            firstStart = start;
        }
    }
}

void FractionalProfiler::endTimer() {
    // We can stop the timing now, as the timer will be left in an invalid state if something goes wrong.
    // The value in stop will then no longer be valid.
    end = chrono::high_resolution_clock::now();
    lastEnd = end;

    if (!running) {
        invalid = true;
        throw runtime_error("The timer is not running.");
    } else if (invalid) {
        throw runtime_error("The timer is in an invalid state, probably the call sequence was not (start end)^n");
    }
    assert(!savedOrDiscarded);
    assert(running);

    running = false;
    int64_t timeDiff = ((chrono::nanoseconds) (end - start)).count();
    assert(timeDiff > 0 && "Timer ended before it started.");
    nsPassed += (uint64_t) timeDiff;
}

uint64_t FractionalProfiler::getTimer() const {
    if (invalid) {
        throw runtime_error("Trying to get an invalid timer.");
    }
    assert(!running);
    return nsPassed;
}

void FractionalProfiler::saveTimer(uint64_t nsMean) {
    if (nsMean == 0) {
        throw runtime_error("An average runtime of 0 ns is a little fast, don't you think?");
    } else if (invalid) {
        throw runtime_error("Trying to save an invalid timer.");
    } else if (running) {
        invalid = true;
        throw runtime_error("Trying to save a running timer.");
    } else if (savedOrDiscarded) {
        invalid = true;
        throw  runtime_error("Timer was already saved (or never ran). Please (re)start if first.");
    }

    savedOrDiscarded = true;
    eventCounter->event((double)nsPassed / nsMean);
}

void FractionalProfiler::discardTimer() {
    if (invalid) {
        throw runtime_error("Invalid timers can't even be discarded, sorry.");
    } if (running) {
        invalid = true;
        throw runtime_error("::discardTimer should be used on ended timers. For running timers, use ::abortTimer");
    } if (savedOrDiscarded && nsPassed == 0) {
        invalid = true;
        throw runtime_error("Are you trying to discard a timer twice in a row?");
    }

    savedOrDiscarded = true;
    nsPassed = 0;
}

void FractionalProfiler::abortTimer() {
    if (invalid) {
        throw runtime_error("Trying to abort invalid timer");
    } else if (!running) {
        invalid = true;
        throw runtime_error("Trying to abort non-running timer. Did you want to call ::discardTimer?");
    }
    running = false;
    assert(!savedOrDiscarded);
    savedOrDiscarded = true;
    // start will be overwritten once the timer starts again
}

FractionalProfiler::~FractionalProfiler() {
    if (!invalid) {
        assert(!running);
        assert(savedOrDiscarded);
    }
}

void FractionalProfiler::resumeTimer() {
    startTimer(true);
}

// shared_ptr<FractionalProfiler::FractionalHistogram> FractionalProfiler::getHistogram() const {
//     if (invalid) {
//         throw runtime_error("Trying to get histogram of invalid timer");
//     }
//     assert(eventCounter != nullptr);
//     return eventCounter;
// }

const string FractionalProfiler::getName() const {
    return this->name;
}

FractionalProfiler::FractionalHistogram::operator std::string () {
    string str = "";
    for (int i = 0; i < numBins; i++) {
        str += to_string((*bins)[i]);
        if (i != numBins - 1) {
            str += ",";
        }
    }
    return str;
}

uint64_t FractionalProfiler::FractionalHistogram::numEvents() const {
    assert(bins != nullptr);
    assert(bins->size() == numBins);
    return accumulate(bins->begin(), bins->end(), 0);
}

const shared_ptr<vector<uint64_t>> FractionalProfiler::FractionalHistogram::data() const {
    assert(bins != nullptr);
    return bins;
}

// unique_ptr<ostream> FractionalProfiler::writeTimingsHeader(unique_ptr<ostream> file) {
//    if (file == nullptr) {
//        throw runtime_error("I will not write to a nullptr.");
//    }

//    *file << "rank,processor,timer,secondsPassed,bin,count" << endl;
//    return file;
// }

// unique_ptr<ostream> FractionalProfiler::writeTimingsStats(shared_ptr<vector<uint64_t>> data, unique_ptr<ostream> file, const string& timerName,
//                                                   string (*rankToProcessorName) (size_t), int secondsPassed) {
//    if (data == nullptr || file == nullptr) {
//        throw runtime_error("nullptr as data or file");
//    }
//    if (secondsPassed < 0) {
//        throw runtime_error("You seem to have invented time travel ;-) (secondsPassed < 0)");
//    } 

//    // Output data
//    for (size_t rank = 0; rank < data->size() / FractionalHistogram::numBins; rank++) {
//        for (uint16_t bin = 0; bin < FractionalHistogram::numBins; bin++) {
//            *file << to_string(rank) << "," << (*rankToProcessorName)(rank) << "," << timerName << ","<< to_string(secondsPassed) << ",";
//            *file << FractionalHistogram::binName(bin) << ",";
//            *file << (*data)[rank * FractionalHistogram::numBins + bin];
//            *file << endl;
//        } 
//    }

//    return file;
// }

uint64_t FractionalProfiler::timesCalled() const {
    if (invalid) {
        throw runtime_error("Trying to get event count on invalid timer.");
    }
    return eventCounter->numEvents();
}

float FractionalProfiler::eventsPerSecond() const {
    if (invalid) {
        throw runtime_error("Trying to get event frequency on invalid timer.");
    }
    return timesCalled() / secondsPassed();
}


bool FractionalProfiler::isRunning() const {
    if (invalid) {
        throw runtime_error("Timer is in invalid state.");
    }
    return running;
}

float FractionalProfiler::secondsPassed() const {
    if (invalid) {
        throw runtime_error("Trying to get passed time on invalid timer.");
    }
    return (lastEnd - firstStart).count() / pow(10, 9);
}

ProfilerRegister::ProfilerRegister(string prefix) {
    createProFile(prefix);
    profilers = make_shared<map<string, shared_ptr<FractionalProfiler>>>();
    stats = make_shared<map<string, Measurement>>();
    //stats->insert({ string("work"), Measurement() });

    //assert(stats->at("work").count == 0 && stats->at("work").nsSum == 0);
    assert(profilers != nullptr);
    assert(stats != nullptr);
}

void ProfilerRegister::createProFile(string prefix) {
    RAXML_UNUSED(prefix);
    //if (proFile != nullptr || callsPerSecondFile != nullptr) {
    if (callsPerSecondFile != nullptr) {
        throw runtime_error("Profiling data logfile has already been created!");
    }

    // ProFile
    //ofstream* file = new ofstream();
    //file->open(prefix + ".proFile.csv");
    //proFile = unique_ptr<ostream>(file);
    //proFile = FractionalProfiler::writeTimingsHeader(move(proFile));
    //assert(proFile != nullptr && *proFile);

    // Calls per second file
    //auto file = new ofstream();
    //file->open(prefix + ".callsPerSecond.csv");
    //callsPerSecondFile = unique_ptr<ostream>(file);
    //callsPerSecondFile = FractionalProfiler::writeCallsPerSecondsHeader(move(callsPerSecondFile));
    //assert(callsPerSecondFile != nullptr && *callsPerSecondFile);

    // Overall statistics file
    // Here, only the filename is computed, the file is reopened (and truncated) on every write.
    overallStatsFilename = "overallStats.csv";

    // Work by rank file
    //file = new ofstream();
    //file->open(prefix + ".workByRank.csv");
    //workByRankFile = unique_ptr<ostream>(file);
    //*workByRankFile << "time,rank,workMs" << endl;
    //assert(workByRankFile != nullptr && *workByRankFile);
}

shared_ptr<FractionalProfiler> ProfilerRegister::registerProfiler(string name) {
    if (profilers->find(name) != profilers->end()) {
        throw runtime_error("A profiler with this name already exists.");
    }
    // This would not throw an exception if a profiler with this name would already exist.
    (*profilers)[name] = make_shared<FractionalProfiler>(name);
    return (*profilers)[name];
}

shared_ptr<FractionalProfiler> ProfilerRegister::getProfiler(string name) const {
    // Will throw an exception if the profiler is non-existent
    return profilers->at(name);
}

shared_ptr<ProfilerRegister> ProfilerRegister::singleton = nullptr;

shared_ptr<ProfilerRegister> ProfilerRegister::getInstance() {
    if (singleton == nullptr) {
        throw runtime_error("No instance has been created yet.");
    }
    return singleton;
}

shared_ptr<ProfilerRegister> ProfilerRegister::createInstance(string logFile) {
    if (singleton != nullptr) {
        throw runtime_error("An instance has already been created.");
    }
    singleton = shared_ptr<ProfilerRegister>(new ProfilerRegister(logFile));
    return singleton;
}

//void ProfilerRegister::startWorkTimer() {
//    assert(!workTimerRunning);
//    workTimerRunning = true;
//
//    workStart = chrono::high_resolution_clock::now();
//}
//
//void ProfilerRegister::endWorkTimer() {
//    assert(stats->find("work") != stats->end());
//    assert(workTimerRunning);
//    workTimerRunning = false;
//
//    auto workEnd = chrono::high_resolution_clock::now();
//    uint64_t workDuration = chrono::duration_cast<chrono::nanoseconds>(workEnd - workStart).count();
//
//    stats->at("work").count++;
//    stats->at("work").nsSum += workDuration;
//}
//
//void ProfilerRegister::discardWorkTimer() {
//    workTimerRunning = false;
//}
//
//void ProfilerRegister::reset_worked_for() {
//    assert(stats->find("work") != stats->end());
//    stats->at("work").count = 0;
//    stats->at("work").nsSum = 0;
//}
//
//double ProfilerRegister::worked_for_ms() {
//    assert(!workTimerRunning);
//    assert(stats->find("work") != stats->end());
//    return (double)stats->at("work").nsSum / (1000*1000);
//}

// void ProfilerRegister::saveProfilingData(bool master, size_t num_ranks, string (*rankToProcessorName) (size_t), MPI_Comm comm) {
//    // Save the time passed once so it will be the same for all measurements collected during one call to this function
// //    int secondsPassed = profilers->begin()->second->secondsPassed();

//    // First, gather the profiling data for all timers
//    // Iterating over a std::map is actually sorted by key -> We will always collect the correct timing data
//    for (auto timerPair: *profilers) {
//        auto timer = timerPair.second;
//        if (timer->isRunning()) {
//            throw runtime_error("Trying to save profiling data while timer " + timer->getName() + " is still running!");
//        }

//        // Get timer histogram
//        //shared_ptr<vector<uint64_t>> timings = timer->getHistogram()->data();
       
//        // Gather data from all ranks and save to file
//        if (master) {
//            //auto allTimingsVec = make_shared<vector<uint64_t>>(num_ranks * timings->size());
//            //MPI_Gather(timings->data(), timings->size(), MPI_UINT64_T,
//            //        allTimingsVec->data(), FractionalProfiler::FractionalHistogram::numBins, MPI_UINT64_T,
//            //        0, comm);

//            //proFile = FractionalProfiler::writeTimingsStats(allTimingsVec, move(proFile), timer->getName(), rankToProcessorName, secondsPassed);
//            //callsPerSecondFile = FractionalProfiler::writeCallsPerSecondsStats(timer->getName(), secondsPassed, timer->eventsPerSecond(), move(callsPerSecondFile));
//        } else {
//            //MPI_Gather(timings->data(), timings->size(), MPI_UINT64_T,
//            //        nullptr, 0, MPI_UINT64_T, 0, comm);
//        }
//    }

// }

//unique_ptr<ostream> FractionalProfiler::writeCallsPerSecondsHeader(unique_ptr<ostream> file) {
//    if (file == nullptr) {
//        throw runtime_error("Stream ptr should not be a nullptr");
//    }
//    *file << "timer,secondsPassed,callsPerSecond" << endl;
//
//    return file;
//}
//
//unique_ptr<ostream> FractionalProfiler::writeCallsPerSecondsStats(string timer, int secondsPassed, float callsPerSecond, unique_ptr<ostream> file) {
//    if (file == nullptr) {
//        throw runtime_error("Ouput stream pointer is nullptr");
//    } else if (secondsPassed < 0) {
//        throw runtime_error("< 0 seconds passed");
//    } else if (callsPerSecond < 0) {
//        throw runtime_error("< 0 calls per second");
//    }
//    *file << timer << "," << to_string(secondsPassed) << "," << to_string(callsPerSecond) << endl;
//
//    return file;
//}

unique_ptr<ostream> ProfilerRegister::writeStatsHeader(unique_ptr<ostream> file) {
    assert(*file);

    *file << "rank,processor,";
    auto stats = getStats();

    size_t field = 0;
    for (auto& stat: *stats) {
        *file << "num" << stat.first << ",nsSum" << stat.first << (field != stats->size() - 1 ? "," : "\n");
        field++;
    }

    return file;
}

void ProfilerRegister::writeStats(string (*rankToProcessorName) (size_t)) {
    unique_ptr<ostream> file;

    // Iterating over a std::map is actually sorted by key -> We will always collect the correct timing data
    if (ParallelContext::master()) {
        assert(overallStatsFilename != "");
        file = unique_ptr<ostream>(new ofstream(overallStatsFilename));

        file = writeStatsHeader(move(file));
        assert(*file);

        auto myRankId = ParallelContext::rank_id();
        // ParallelContext::mpi_gather_custom does not send this ranks data to itself
        *file << myRankId << "," << rankToProcessorName(myRankId) << ",";
        size_t statId = 0;
        auto stats = getStats();
        for (auto& stat: *stats) {
            *file << stat.second.count << "," << stat.second.nsSum << (statId != stats->size() - 1 ? "," : "\n");
            statId++;
        }
    }

    // Gather data from all ranks and save to file
    auto stats = getStats();
    auto prepare_send_cb = [&stats](void * buf, size_t buf_size) -> int
    {
        assert(!ParallelContext::master());
        BinaryStream bs((char*) buf, buf_size);
        bs << ParallelContext::rank_id();
        bs << stats->size();
        for (auto& stat: *stats) {
            bs << stat.second.count;
            bs << stat.second.nsSum;
        }
        return (int) bs.pos();
    };

    auto process_recv_cb = [&stats, &file, &rankToProcessorName](void * buf, size_t buf_size)
    {
        assert(ParallelContext::master());
        BinaryStream bs((char*) buf, buf_size);
        auto rankId = bs.get<size_t>();
        auto numberOfStats = bs.get<size_t>();
        assert(numberOfStats == stats->size()); // The same number of stats have been collected on all ranks

        assert(*file);
        *file << rankId << "," << rankToProcessorName(rankId) << ",";
        for (uint64_t statId = 0; statId < numberOfStats; statId++) {
            auto count = bs.get<size_t>();
            auto nsSum = bs.get<size_t>();
            *file << count << "," << nsSum << (statId != numberOfStats - 1 ? "," : "\n");
        }
    };

    ParallelContext::mpi_gather_custom(prepare_send_cb, process_recv_cb);
}

//const string FractionalProfiler::FractionalHistogram::binName(uint16_t bin) {
//    auto to_string_with_precision = [](const float value, const int n = 3) {
//        std::ostringstream out;
//        out.precision(n);
//        out << std::fixed << value;
//        return out.str();
//    };
//
//    // See getBin() for memory layout
//    if (bin < 100) {
//        float oneOverFrac = (float)bin / 10 + 1;
//        float to = 1 / oneOverFrac;
//        float from = 1 / (oneOverFrac + 0.1);
//        return to_string_with_precision(from) + " to " + to_string_with_precision(to);
//    } else if (bin == 100) {
//        return "0.999 to 1.001";
//    } else {
//        float from = ((float)bin - 101) / 10 + 1;
//        float to = from + 0.1;
//        return to_string_with_precision(from) + " to " + to_string_with_precision(to);
//    } 
//}

shared_ptr<map<string, ProfilerRegister::Measurement>> ProfilerRegister::getStats() {
    assert(stats != nullptr);
    return stats;
}

//void ProfilerRegister::saveWorkByRank(bool reset) {
//    auto work = work_by_rank(); // Will perform an MPI_Allgather
//    if (reset) {
//        reset_worked_for();
//    }
//
//    if (!ParallelContext::master()) {
//        return;
//    }
//    
//    assert(workByRankFile != nullptr && *workByRankFile);
//    
//    for (size_t rank = 0; rank < work->size(); rank++) {
//        *workByRankFile << num_rebalances << "," << rank << "," << work->at(rank) << endl;
//    }
//    assert(num_rebalances < numeric_limits<decltype(num_rebalances)>::max());
//    num_rebalances++;
//}
//
//std::shared_ptr<vector<double>> ProfilerRegister::work_by_rank() {
//    double local_work = worked_for_ms();
//    shared_ptr<doubleVector> work = ParallelContext::mpi_allgather(local_work);
//    assert(work != nullptr);
//    assert(work->size() == ParallelContext::num_ranks());
//    return work;
//}