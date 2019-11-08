#include <memory>
#include <chrono>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>
#include <fstream>
#include <numeric>
#include <cmath>

#include "Profiler.hpp"

using namespace std;

LogBinningProfiler::LogarithmicHistogram::LogarithmicHistogram()
    : numBins(64) {
    bins = std::make_shared<vector<uint64_t>>(numBins, 0);
}

void LogBinningProfiler::LogarithmicHistogram::event(uint64_t number) {
    const uint8_t bin = log2i(number);
    (*bins)[bin]++;
}

uint64_t& LogBinningProfiler::LogarithmicHistogram::operator[](size_t idx) {
    return (*bins)[idx];
} 

const uint64_t& LogBinningProfiler::LogarithmicHistogram::operator[](size_t idx) const {
    return (*bins)[idx];
}

uint8_t LogBinningProfiler::LogarithmicHistogram::log2i(uint64_t n) {
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

void LogBinningProfiler::startTimer() {
    if (running) {
        invalid = true;
        throw runtime_error("The timer is already running.");
    } else if (invalid) {
        throw runtime_error("The timer is in an invalid state, probably the call sequence was not (start end)^n");
    }
    running = true;
    start = chrono::high_resolution_clock::now();
    if (firstStart == chrono::time_point<chrono::high_resolution_clock>::max()) {
        firstStart = start;
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

    running = false;
    int64_t timeDiff = ((chrono::nanoseconds) (end - start)).count();
    assert(timeDiff > 0 && "Timer ended before it started.");
    eventCounter->event((uint64_t) timeDiff);
}

shared_ptr<LogBinningProfiler::LogarithmicHistogram> LogBinningProfiler::getHistogram() const {
    if (invalid) {
        throw runtime_error("Trying to get histogram of invalid timer");
    }
    return eventCounter;
}

const string LogBinningProfiler::getName() const {
    return this->name;
}

LogBinningProfiler::LogarithmicHistogram::operator std::string () {
    string str = "";
    for (int i = 0; i < 63; i++) {
        str += to_string((*bins)[i]) + ",";
    }
    str += to_string((*bins)[63]);
    return str;
}

uint64_t LogBinningProfiler::LogarithmicHistogram::numEvents() const {
    return accumulate(bins->begin(), bins->end(), 0);
}

const shared_ptr<vector<uint64_t>> LogBinningProfiler::LogarithmicHistogram::data() const {
    return bins;
}

void LogBinningProfiler::writeStats(shared_ptr<vector<uint64_t>> data, shared_ptr<ostream> file, const string& timerName, bool printHeader) {
    if (data == nullptr || file == nullptr) {
        throw runtime_error("nullptr as data or file");
    }

    // Output header
    if (printHeader) {
        *file << "rank,timer,";
        for (int bit = 0; bit < 64; bit++) {
            *file << "[2^" << bit << ",2^" << bit + 1 << ") ns";
            if (bit != 63) {
                *file << ",";
            } else {
                *file << endl;
            }
        }
    }

    // Output data
    for (size_t rank = 0; rank < data->size() / 64; rank++) {
        *file << to_string(rank) + "," + timerName + ",";
        for (int timing = 0; timing < 64; timing++) {
            *file << (*data)[rank * 64 + timing];
            if (timing != 63) {
                *file << ",";
            } else {
                *file << endl;
            }
        } 
    }
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
    float passedSeconds = (lastEnd - firstStart).count() / pow(10, 9);
    return timesCalled() / passedSeconds;
}

void LogBinningProfiler::abortTimer() {
    if (invalid) {
        throw runtime_error("Trying to abort invalid timer");
    } else if (!running) {
        throw runtime_error("Trying to abort non-running timer.");
    }
    running = false;
    // start will be overwritten once the timer starts again
}


bool LogBinningProfiler::isRunning() const {
    if (invalid) {
        throw runtime_error("Timer is in invalid state.");
    }
    return running;
}