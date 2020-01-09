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

shared_ptr<FractionalProfiler::FractionalHistogram> FractionalProfiler::getHistogram() const {
    if (invalid) {
        throw runtime_error("Trying to get histogram of invalid timer");
    }
    assert(eventCounter != nullptr);
    return eventCounter;
}

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

unique_ptr<ostream> FractionalProfiler::writeTimingsHeader(unique_ptr<ostream> file) {
    if (file == nullptr) {
        throw runtime_error("I will not write to a nullptr.");
    }

    *file << "rank,processor,timer,secondsPassed,bin,count" << endl;
    return file;
}

unique_ptr<ostream> FractionalProfiler::writeTimingsStats(shared_ptr<vector<uint64_t>> data, unique_ptr<ostream> file, const string& timerName,
                                                   string (*rankToProcessorName) (size_t), int secondsPassed) {
    if (data == nullptr || file == nullptr) {
        throw runtime_error("nullptr as data or file");
    }
    if (secondsPassed < 0) {
        throw runtime_error("You seem to have invented time travel ;-) (secondsPassed < 0)");
    } 

    // Output data
    for (size_t rank = 0; rank < data->size() / FractionalHistogram::numBins; rank++) {
        for (uint16_t bin = 0; bin < FractionalHistogram::numBins; bin++) {
            *file << to_string(rank) << "," << (*rankToProcessorName)(rank) << "," << timerName << ","<< to_string(secondsPassed) << ",";
            *file << FractionalHistogram::binName(bin) << ",";
            *file << (*data)[rank * FractionalHistogram::numBins + bin];
            *file << endl;
        } 
    }

    return file;
}

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
    stats = make_shared<ProfilerStats>();

    assert(profilers != nullptr);
    assert(stats != nullptr);
}

void ProfilerRegister::createProFile(string prefix) {
    if (proFile != nullptr || callsPerSecondFile != nullptr) {
        throw runtime_error("Profiling data logfile has already been created!");
    }

    // ProFile
    ofstream* file = new ofstream();
    file->open(prefix + ".proFile.csv");
    proFile = unique_ptr<ostream>(file);
    proFile = FractionalProfiler::writeTimingsHeader(move(proFile));

    // Calls per second file
    file = new ofstream();
    file->open(prefix + ".callsPerSecond.csv");
    callsPerSecondFile = unique_ptr<ostream>(file);
    callsPerSecondFile = FractionalProfiler::writeCallsPerSecondsHeader(move(callsPerSecondFile));

    // Overall statistics file
    // Here, only the filename is computed, the file is reopened (and truncated) on every write.
    overallStatsFilename = prefix + ".overallStats.csv";
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

void ProfilerRegister::saveProfilingData(bool master, size_t num_ranks, string (*rankToProcessorName) (size_t), MPI_Comm comm) {
    // Save the time passed once so it will be the same for all measurements collected during one call to this function
    int secondsPassed = profilers->begin()->second->secondsPassed();

    // First, gather the profiling data for all timers
    // Iterating over a std::map is actually sorted by key -> We will always collect the correct timing data
    for (auto timerPair: *profilers) {
        auto timer = timerPair.second;
        if (timer->isRunning()) {
            throw runtime_error("Trying to save profiling data while timer " + timer->getName() + " is still running!");
        }

        // Get timer histogram
        shared_ptr<vector<uint64_t>> timings = timer->getHistogram()->data();
        
        // Gather data from all ranks and save to file
        if (master) {
            auto allTimingsVec = make_shared<vector<uint64_t>>(num_ranks * timings->size());
            MPI_Gather(timings->data(), timings->size(), MPI_UINT64_T,
                    allTimingsVec->data(), FractionalProfiler::FractionalHistogram::numBins, MPI_UINT64_T,
                    0, comm);

            proFile = FractionalProfiler::writeTimingsStats(allTimingsVec, move(proFile), timer->getName(), rankToProcessorName, secondsPassed);
            callsPerSecondFile = FractionalProfiler::writeCallsPerSecondsStats(timer->getName(), secondsPassed, timer->eventsPerSecond(), move(callsPerSecondFile));
        } else {
            MPI_Gather(timings->data(), timings->size(), MPI_UINT64_T,
                    nullptr, 0, MPI_UINT64_T, 0, comm);
        }
    }

    // Second, gather all profiler stats
    // Arrange profiler stats data
    auto profilerStats = ProfilerRegister::getInstance()->getStats();
    vector<uint64_t> statsVec = {
        profilerStats->nsSumInsideMPI,
        profilerStats->nsSumOutsideMPI,
        profilerStats->nsSumWait,
        profilerStats->nsSumWork,
        profilerStats->numIterations,
        profilerStats->timesIWasSlowest
    };

    if (master) {
        auto allStatsVec = make_shared<vector<uint64_t>>(num_ranks * statsVec.size());
        MPI_Gather(statsVec.data(), statsVec.size(), MPI_UINT64_T,
            allStatsVec->data(), statsVec.size(), MPI_UINT64_T, 0, comm);

        ofstream* overallStatsFile = new ofstream(overallStatsFilename, ofstream::out | ofstream::trunc);
        FractionalProfiler::writeOverallStats(allStatsVec, rankToProcessorName, move(unique_ptr<ostream>(overallStatsFile)));
    } else {
        MPI_Gather(statsVec.data(), statsVec.size(), MPI_UINT64_T,
                nullptr, 0, MPI_UINT64_T, 0, comm);
    }
}

unique_ptr<ostream> FractionalProfiler::writeCallsPerSecondsHeader(unique_ptr<ostream> file) {
    if (file == nullptr) {
        throw runtime_error("Stream ptr should not be a nullptr");
    }
    *file << "timer,secondsPassed,callsPerSecond" << endl;

    return file;
}

unique_ptr<ostream> FractionalProfiler::writeCallsPerSecondsStats(string timer, int secondsPassed, float callsPerSecond, unique_ptr<ostream> file) {
    if (file == nullptr) {
        throw runtime_error("Ouput stream pointer is nullptr");
    } else if (secondsPassed < 0) {
        throw runtime_error("< 0 seconds passed");
    } else if (callsPerSecond < 0) {
        throw runtime_error("< 0 calls per second");
    }
    *file << timer << "," << to_string(secondsPassed) << "," << to_string(callsPerSecond) << endl;

    return file;
}

unique_ptr<ostream> FractionalProfiler::writeOverallStats(shared_ptr<vector<uint64_t>> allStatsVec, string (*rankToProcessorName) (size_t), unique_ptr<ostream> file) {
    assert(allStatsVec->size() % 6 == 0);
    size_t numRanks = allStatsVec->size() / 6;
    // See saveProfilingData(...) for order of fields
    *file << "rank,processor,nsSumInsideMPI,nsSumOutsideMPI,nsSumWait,nsSumWork,numIterations,timesIWasSlowest" << endl;
    for (size_t rank = 0; rank < numRanks; rank++) {
        size_t base = rank * 6;
        *file << rank << ",";
        *file << (*rankToProcessorName)(rank) << ",";
        for (uint8_t i = 0; i < 6; i++) {
            *file << to_string((*allStatsVec)[base + i]);
            if (i != 5) {
                *file << ",";
            }
        }
        *file << endl;
    }

    return file;
}

const string FractionalProfiler::FractionalHistogram::binName(uint16_t bin) {
    auto to_string_with_precision = [](const float value, const int n = 3) {
        std::ostringstream out;
        out.precision(n);
        out << std::fixed << value;
        return out.str();
    };

    // See getBin() for memory layout
    if (bin < 100) {
        float oneOverFrac = (float)bin / 10 + 1;
        float to = 1 / oneOverFrac;
        float from = 1 / (oneOverFrac + 0.1);
        return to_string_with_precision(from) + " to " + to_string_with_precision(to);
    } else if (bin == 100) {
        return "0.999 to 1.001";
    } else {
        float from = ((float)bin - 101) / 10 + 1;
        float to = from + 0.1;
        return to_string_with_precision(from) + " to " + to_string_with_precision(to);
    } 
}

shared_ptr<ProfilerRegister::ProfilerStats> ProfilerRegister::getStats() {
    assert(stats != nullptr);
    return stats;
}
