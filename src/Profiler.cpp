#include <memory>
#include <chrono>
#include <string>
#include <stdexcept>
#include <vector>
#include <cassert>

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
}

void LogBinningProfiler::endTimer() {
    // We can stop the timing now, as the timer will be left in an invalid state if something goes wrong.
    // The value in stop will then no longer be valid.
    end = chrono::high_resolution_clock::now();

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

const shared_ptr<LogBinningProfiler::LogarithmicHistogram> LogBinningProfiler::getHistogram() const {
    return eventCounter;
}

const string LogBinningProfiler::getName() const {
    return this->name;
}