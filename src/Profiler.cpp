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

#include "Profiler.hpp"

using namespace std;

LogBinningProfiler::LogarithmicHistogram::LogarithmicHistogram() {
    bins = std::make_shared<vector<uint64_t>>(numBins, 0);
}

void LogBinningProfiler::LogarithmicHistogram::event(uint64_t number) {
    assert(bins != nullptr);

    int16_t bin = log2i(number); // Returns -1 if number == 0
    bin++; // Shift bins by one to be able to save 0
    assert(bin >= 0);
    assert(bin < numBins);
    assert(numBins == bins->size());

    (*bins)[bin]++;
}

uint64_t& LogBinningProfiler::LogarithmicHistogram::operator[](size_t idx) {
    return (*bins)[idx];
} 

const uint64_t& LogBinningProfiler::LogarithmicHistogram::operator[](size_t idx) const {
    return (*bins)[idx];
}

int16_t LogBinningProfiler::LogarithmicHistogram::log2i(uint64_t n) {
    if (n == 0) {
        return -1;
    }

    uint8_t msb = 0;
    while (n >>= 1) {
        msb++;
    }    

    return msb;
}

LogBinningProfiler::LogBinningProfiler(string name) : name(name) {}

void LogBinningProfiler::startTimer(bool resume) {
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

void LogBinningProfiler::endTimer() {
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

uint64_t LogBinningProfiler::getTimer() const {
    if (invalid) {
        throw runtime_error("Trying to get an invalid timer.");
    }
    assert(!running);
    return nsPassed;
}

void LogBinningProfiler::saveTimer(uint64_t min) {
    if (invalid) {
        throw runtime_error("Trying to save an invalid timer.");
    } else if (running) {
        invalid = true;
        throw runtime_error("Trying to save a running timer.");
    } else if (savedOrDiscarded) {
        invalid = true;
        throw  runtime_error("Timer was already saved (or never ran). Please (re)start if first.");
    }

    assert(nsPassed >= min);
    savedOrDiscarded = true;
    eventCounter->event(nsPassed - min);
}

void LogBinningProfiler::discardTimer() {
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

void LogBinningProfiler::abortTimer() {
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

LogBinningProfiler::~LogBinningProfiler() {
    if (!invalid) {
        assert(!running);
        assert(savedOrDiscarded);
    }
}

void LogBinningProfiler::resumeTimer() {
    startTimer(true);
}

shared_ptr<LogBinningProfiler::LogarithmicHistogram> LogBinningProfiler::getHistogram() const {
    if (invalid) {
        throw runtime_error("Trying to get histogram of invalid timer");
    }
    assert(eventCounter != nullptr);
    return eventCounter;
}

const string LogBinningProfiler::getName() const {
    return this->name;
}

LogBinningProfiler::LogarithmicHistogram::operator std::string () {
    string str = "";
    for (int i = 0; i < numBins; i++) {
        str += to_string((*bins)[i]);
        if (i != numBins - 1) {
            str += ",";
        }
    }
    return str;
}

uint64_t LogBinningProfiler::LogarithmicHistogram::numEvents() const {
    assert(bins != nullptr);
    assert(bins->size() == numBins);
    return accumulate(bins->begin(), bins->end(), 0);
}

const shared_ptr<vector<uint64_t>> LogBinningProfiler::LogarithmicHistogram::data() const {
    assert(bins != nullptr);
    return bins;
}

unique_ptr<ostream> LogBinningProfiler::writeStatsHeader(unique_ptr<ostream> file) {
    if (file == nullptr) {
        throw runtime_error("I will not write to a nullptr.");
    }

    *file << "rank,processor,timer,secondsPassed,0 ns,";
    for (int bit = 0; bit < LogarithmicHistogram::numBins - 1; bit++) {
        *file << "\"[2^";
        *file << setfill('0') << setw(2) << bit;
        *file << ",2^";
        *file << setfill('0') << setw(2) << bit + 1;
        *file << ") ns\"";
        if (bit != LogarithmicHistogram::numBins - 2) {
            *file << ",";
        } else {
            *file << endl;
        }
    }
    return file;
}

unique_ptr<ostream> LogBinningProfiler::writeStats(shared_ptr<vector<uint64_t>> data, unique_ptr<ostream> file, const string& timerName,
                                                   string (*rankToProcessorName) (size_t), int secondsPassed) {
    if (data == nullptr || file == nullptr) {
        throw runtime_error("nullptr as data or file");
    }
    if (secondsPassed < 0) {
        throw runtime_error("You seem to have invented time travel ;-) (secondsPassed < 0)");
    } 

    // Output data
    for (size_t rank = 0; rank < data->size() / LogarithmicHistogram::numBins; rank++) {
        *file << to_string(rank) + "," + (*rankToProcessorName)(rank) + "," + timerName + "," + to_string(secondsPassed) + ",";
        for (int timing = 0; timing < LogarithmicHistogram::numBins; timing++) {
            *file << (*data)[rank * LogarithmicHistogram::numBins + timing];
            if (timing != LogarithmicHistogram::numBins - 1) {
                *file << ",";
            } else {
                *file << endl;
            }
        } 
    }

    return file;
}

uint64_t LogBinningProfiler::timesCalled() const {
    if (invalid) {
        throw runtime_error("Trying to get event count on invalid timer.");
    }
    return eventCounter->numEvents();
}

float LogBinningProfiler::eventsPerSecond() const {
    if (invalid) {
        throw runtime_error("Trying to get event frequency on invalid timer.");
    }
    return timesCalled() / secondsPassed();
}


bool LogBinningProfiler::isRunning() const {
    if (invalid) {
        throw runtime_error("Timer is in invalid state.");
    }
    return running;
}

float LogBinningProfiler::secondsPassed() const {
    if (invalid) {
        throw runtime_error("Trying to get passed time on invalid timer.");
    }
    return (lastEnd - firstStart).count() / pow(10, 9);
}

ProfilerRegister::ProfilerRegister(string prefix) {
    createProFile(prefix);
    profilers = make_shared<map<string, shared_ptr<LogBinningProfiler>>>();
}

void ProfilerRegister::createProFile(string prefix) {
    if (proFile != nullptr || callsPerSecondFile != nullptr) {
        throw runtime_error("Profiling data logfile has already been created!");
    }

    // ProFile
    ofstream* file = new ofstream();
    file->open(prefix + "_proFile.csv");
    proFile = unique_ptr<ostream>(file);
    proFile = LogBinningProfiler::writeStatsHeader(move(proFile));

    // Calls per second file
    file = new ofstream();
    file->open(prefix + "_callsPerSecond.csv");
    callsPerSecondFile = unique_ptr<ostream>(file);
    callsPerSecondFile = LogBinningProfiler::writeCallsPerSecondsHeader(move(callsPerSecondFile));
}

shared_ptr<LogBinningProfiler> ProfilerRegister::registerProfiler(string name) {
    if (profilers->find(name) != profilers->end()) {
        throw runtime_error("A profiler with this name already exists.");
    }
    // This would not throw an exception if a profiler with this name would already exist.
    (*profilers)[name] = make_shared<LogBinningProfiler>(name);
    return (*profilers)[name];
}

shared_ptr<LogBinningProfiler> ProfilerRegister::getProfiler(string name) const {
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

    // Iterating over a std::map is actually sorted by key -> We will always collect the correct timing data
    for (auto timerPair: *profilers) {
        auto timer = timerPair.second;
        if (timer->isRunning()) {
            throw runtime_error("Trying to save profiling data while timer " + timer->getName() + " is still running!");
        }
        shared_ptr<vector<uint64_t>> timings = timer->getHistogram()->data();
        if (master) {
            auto allTimingsVec = make_shared<vector<uint64_t>>(num_ranks * timings->size());
            MPI_Gather(timings->data(), timings->size(), MPI_UINT64_T,
                    allTimingsVec->data(), LogBinningProfiler::LogarithmicHistogram::numBins, MPI_UINT64_T,
                    0, comm);

            proFile = LogBinningProfiler::writeStats(allTimingsVec, move(proFile), timer->getName(), rankToProcessorName, secondsPassed);
            callsPerSecondFile = LogBinningProfiler::writeCallsPerSecondsStats(timer->getName(), secondsPassed, timer->eventsPerSecond(), move(callsPerSecondFile));
        } else {
            MPI_Gather(timings->data(), timings->size(), MPI_UINT64_T,
                    nullptr, 0, MPI_UINT64_T, 0, comm);
        }
    }
}

unique_ptr<ostream> LogBinningProfiler::writeCallsPerSecondsHeader(unique_ptr<ostream> file) {
    if (file == nullptr) {
        throw runtime_error("Stream ptr should not be a nullptr");
    }
    *file << "timer,secondsPassed,callsPerSecond" << endl;

    return file;
}

unique_ptr<ostream> LogBinningProfiler::writeCallsPerSecondsStats(string timer, int secondsPassed, float callsPerSecond, unique_ptr<ostream> file) {
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