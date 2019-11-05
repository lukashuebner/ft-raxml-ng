#include <chrono>
#include <thread>
#include <iostream> 

#include "RaxmlTest.hpp"
#include "src/Profiler.hpp"

using namespace std;

TEST(ProfilerTest, LogHistogramBasics) {
    // buildup
    auto hist = LogBinningProfiler::LogarithmicHistogram();

    for (int i = 0; i < 64; i++) {
        ASSERT_EQ(hist[i], 0);
    }

    hist.event(1);
    ASSERT_EQ(hist[0], 1);
    ASSERT_EQ(hist[1], 0);
    for (int i = 2; i < 64; i++) {
        ASSERT_EQ(hist[i], 0);
    }


    hist.event(2);
    hist.event(3);
    ASSERT_EQ(hist[1], 2);

    for (int i = 0; i < 1337; i++) {
        hist.event(1 << 10);
    }
    ASSERT_EQ(hist[10], 1337);
}

TEST(ProfilerTest, Timer) {
    // buildup
    auto profiler = LogBinningProfiler("Testing");

    // tests
    profiler.getName() == "Testing";

    auto hist = profiler.getHistogram();
    for (int i = 0; i < 64; i++) {
        ASSERT_EQ((*hist)[i], 0);
    }

    for (int i = 0; i < 100; i++) {
        profiler.startTimer();
        std::this_thread::sleep_for(chrono::nanoseconds((1 << 20) + (1 << 8)));
        profiler.endTimer();
    }
    hist = profiler.getHistogram();

    for (int i = 0; i < 64; i++) {
        if (i == 20) {
            ASSERT_GE((*hist)[i], 70);
        } else {
            ASSERT_LE((*hist)[i], 5);
        }
    }
}