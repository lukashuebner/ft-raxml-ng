#include <chrono>
#include <thread>
#include <iostream> 
#include <cstdio>

#include "RaxmlTest.hpp"
#include "../../src/Profiler.hpp"

using namespace std;

TEST(ProfilerTest, FractionalHistogramBasics) {
    // buildup
    auto hist = FractionalProfiler::FractionalHistogram();

    for (int i = 0; i < FractionalProfiler::FractionalHistogram::numBins; i++) {
        ASSERT_EQ(hist[i], 0);
    }

    hist.event(1);
    hist.event(1.0001);
    for (int i = 0; i < FractionalProfiler::FractionalHistogram::numBins; i++) {
        if (i == 100) {
            ASSERT_EQ(hist[i], 2);
        } else {
            ASSERT_EQ(hist[i], 0);
        }
    }

    hist.event(2);
    hist.event(3);
    ASSERT_EQ(hist[111], 1);
    ASSERT_EQ(hist[121], 1);

    for (int i = 0; i < 1337; i++) {
        hist.event(1.05);
    }
    ASSERT_EQ(hist[101], 1337);

    hist.event(1/4.5);
    ASSERT_EQ(hist[35], 1);
    hist.event(11);
    ASSERT_EQ(hist[200], 1);
    hist.event(0.05);
    ASSERT_EQ(hist[99], 1);

    ASSERT_EQ(hist.numEvents(), 2 + 2 + 1337 + 3);
}

TEST(ProfilerTest, InvalidArguments) {
    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return to_string(*(rank == 0 ? "Schneewittchen" : "Zwerg"));
    };
    unique_ptr<ostream> output = unique_ptr<ostream>(new ostringstream());

    auto profiler = FractionalProfiler("Testing");
    ASSERT_ANY_THROW(output = profiler.writeTimingsStats(nullptr, move(output), "hello_world()", dummyLabeller, 0));
    ASSERT_ANY_THROW(profiler.writeTimingsStats(profiler.getHistogram()->data(), nullptr, "bye_world()", dummyLabeller, 0));
    ASSERT_ANY_THROW(output = profiler.writeTimingsStats(profiler.getHistogram()->data(), move(output), "bye_world()", dummyLabeller, -1));

    ASSERT_ANY_THROW(FractionalProfiler::writeCallsPerSecondsHeader(nullptr));
    ASSERT_ANY_THROW(FractionalProfiler::writeCallsPerSecondsStats("hello_world", -1, 1, move(output)));
    output = unique_ptr<ostream>(new ostringstream());
    ASSERT_ANY_THROW(FractionalProfiler::writeCallsPerSecondsStats("hello_world", 1, -1, move(output)));
    output = unique_ptr<ostream>(new ostringstream());
    ASSERT_ANY_THROW(FractionalProfiler::writeCallsPerSecondsStats("hello_world", -1, 1, nullptr));

    ASSERT_NO_THROW(output = FractionalProfiler::writeCallsPerSecondsStats("", 0, 0, move(output)));
}

TEST(ProfilerTest, TimerBasics) {
    // buildup
    auto profiler = FractionalProfiler("Testing");

    // tests
    profiler.getName() == "Testing";

    auto hist = profiler.getHistogram();
    for (int i = 0; i < hist->numBins; i++) {
        ASSERT_EQ((*hist)[i], 0);
    }

    chrono::time_point<chrono::high_resolution_clock> start, end;
    start = chrono::high_resolution_clock::now();
    const unsigned int NUM_EVENTS = 100;
    for (unsigned int i = 0; i < NUM_EVENTS; i++) {
        profiler.startTimer();
        std::this_thread::sleep_for(chrono::nanoseconds((1 << 20) + (1 << 8)));
        profiler.endTimer();
        profiler.saveTimer(profiler.getTimer());
    }
    end = chrono::high_resolution_clock::now();
    float eventsPerSecond = NUM_EVENTS / ((end - start).count() / pow(10, 9));
    hist = profiler.getHistogram();

    for (int i = 0; i < 64; i++) {
        if (i == 100) {
            ASSERT_EQ((*hist)[i], 100);
        } else {
            ASSERT_EQ((*hist)[i], 0);
        }
    }

    ASSERT_EQ(hist->numEvents(), 100);
    ASSERT_GT(profiler.eventsPerSecond(), 100);
    ASSERT_NEAR(profiler.eventsPerSecond(), eventsPerSecond, eventsPerSecond * 0.01);
    ASSERT_EQ(profiler.eventsPerSecond(), profiler.eventsPerSecond());

    profiler.startTimer();
    profiler.endTimer();
    ASSERT_ANY_THROW(profiler.saveTimer(0));
    profiler.saveTimer(1); // Leaves the profiler in a valid state so the destructor won't throw
}

TEST(ProfilerTest, TimerAdvanced) {
    auto profiler = FractionalProfiler("Testing");

    const size_t ONEHALF_EVENTS = 2;
    const size_t ONETHIRD_EVENTS = 4;
    const size_t ONESIXTH_EVENTS = 8;
    const size_t ONEPOINTFIVE_EVENTS = 16;
    const size_t DOUBLE_EVENTS = 32;
    const size_t THREEPOINTSIX_EVENTS = 64; 

    auto registerWait = [&profiler](float frac) {
        profiler.startTimer();
        std::this_thread::sleep_for(chrono::nanoseconds((1 << 20)));
        profiler.endTimer();
        profiler.saveTimer(profiler.getTimer() / frac);
    };

    for (size_t i = 0; i < ONEHALF_EVENTS; i++) {
        registerWait(0.5);
    }
    for (size_t i = 0; i < ONETHIRD_EVENTS; i++) {
        registerWait(0.333); // 1/3
    }
    for (size_t i = 0; i < ONESIXTH_EVENTS; i++) {
        registerWait(0.166); // 1/6
    }
    for (size_t i = 0; i < ONEPOINTFIVE_EVENTS; i++) {
        registerWait(1.5);
    }
    for (size_t i = 0; i < DOUBLE_EVENTS; i++) {
        registerWait(2.0);
    }
    for (size_t i = 0; i < THREEPOINTSIX_EVENTS; i++) {
        registerWait(3.6);
    }

    auto hist = *(profiler.getHistogram());
    for (uint16_t i = 0; i < FractionalProfiler::FractionalHistogram::numBins; i++) {
        if (i == 10) { // 1/2
            ASSERT_EQ(hist[i], ONEHALF_EVENTS);
        } else if (i == 20) {
            ASSERT_EQ(hist[i], ONETHIRD_EVENTS);
        } else if (i == 50) {
            ASSERT_EQ(hist[i], ONESIXTH_EVENTS);
        } else if (i == 106) {
            ASSERT_EQ(hist[i], ONEPOINTFIVE_EVENTS);
        } else if (i == 111) {
            ASSERT_EQ(hist[i], DOUBLE_EVENTS);
        } else if (i == 127) {
            ASSERT_EQ(hist[i], THREEPOINTSIX_EVENTS);
        } else {
            ASSERT_EQ(hist[i], 0);
        }
    }
}

TEST(ProfilerTest, InvalidStates) {
    auto profiler = FractionalProfiler("Testing");
    // Try to stop the timer before it started
    ASSERT_ANY_THROW(profiler.endTimer());
    // After that, most other methods should throw an exception because of invalid state, too.
    ASSERT_ANY_THROW(profiler.startTimer());
    ASSERT_ANY_THROW(profiler.abortTimer());
    ASSERT_ANY_THROW(profiler.eventsPerSecond());
    ASSERT_ANY_THROW(profiler.getHistogram());

    // Try to start the timer twice in a row
    profiler = FractionalProfiler("Testing");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.startTimer());

    // Try to abort a non-running timer
    profiler = FractionalProfiler("Testing");
    ASSERT_ANY_THROW(profiler.abortTimer());

    // Try to end an aborted timer
    profiler = FractionalProfiler("Testing");
    profiler.startTimer();
    profiler.abortTimer();
    ASSERT_ANY_THROW(profiler.endTimer());

    // Try to save a running timer and an invalid timer
    profiler = FractionalProfiler("Testing");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.saveTimer(0));

    profiler = FractionalProfiler("Testing");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.startTimer());
    ASSERT_ANY_THROW(profiler.saveTimer(0));

    // Try to discard an invalid timer
    ASSERT_ANY_THROW(profiler.discardTimer());
}

TEST(ProfilerTest, TestIsRunning) {
    auto profiler = FractionalProfiler("Testing");
    
    ASSERT_FALSE(profiler.isRunning());
    profiler.startTimer();
    ASSERT_TRUE(profiler.isRunning());
    profiler.endTimer();
    ASSERT_FALSE(profiler.isRunning());
    profiler.saveTimer(1);
    ASSERT_FALSE(profiler.isRunning());

    profiler.startTimer();
    profiler.abortTimer();
    ASSERT_FALSE(profiler.isRunning());

    ASSERT_ANY_THROW(profiler.endTimer());
    ASSERT_ANY_THROW(profiler.isRunning());
}

TEST(ProfilerTest, TimerCanBeAborted) {
    auto profiler = FractionalProfiler("Testing");

    profiler.startTimer();
    ASSERT_NO_THROW(profiler.endTimer());
    profiler.saveTimer(1);
    profiler.startTimer();
    ASSERT_NO_THROW(profiler.abortTimer());
    ASSERT_EQ(profiler.getHistogram()->numEvents(), 1);
    ASSERT_NO_THROW(profiler.startTimer());
    profiler.abortTimer();
    ASSERT_ANY_THROW(profiler.endTimer());

    profiler = FractionalProfiler("Testing");
    profiler.startTimer();
    profiler.endTimer();
    ASSERT_ANY_THROW(profiler.abortTimer());
}

const string statsHeader = "rank,processor,timer,secondsPassed,bin,count\n";
const string callsPerSecondsHeader = "timer,secondsPassed,callsPerSecond\n";
TEST(ProfilerTest, writeStats_Empty) {
    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return to_string(*(rank == 0 ? "Schneewittchen" : "Zwerg"));
    };
    auto stream = new ostringstream();
    unique_ptr<ostream> output = unique_ptr<ostream>(stream);
    auto input = make_shared<vector<uint64_t>>();

    output = FractionalProfiler::writeTimingsStats(input, move(output), "hello_world()", dummyLabeller, 0);
    ASSERT_TRUE(output);
    ASSERT_EQ(stream->str(), "");

    output = FractionalProfiler::writeTimingsHeader(move(output));
    ASSERT_TRUE(output);
    ASSERT_EQ(stream->str(), statsHeader);

    stream = new ostringstream();
    output = unique_ptr<ostream>(stream);
    output = FractionalProfiler::writeCallsPerSecondsHeader(move(output));
    ASSERT_TRUE(output);
    ASSERT_EQ(stream->str(), callsPerSecondsHeader);
}

TEST(ProfilerTest, writeStats_Basic) {
    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return (rank == 0 ? string("Schneewittchen") : string("Zwerg"));
    };

    auto stream = new ostringstream();
    unique_ptr<ostream> output = unique_ptr<ostream>(stream);
    auto input = make_shared<vector<uint64_t>>(FractionalProfiler::FractionalHistogram::numBins * 4);

    for (int node = 0; node < 4; node++) {
        for (int timing = 0; timing < FractionalProfiler::FractionalHistogram::numBins; timing++) {
            (*input)[FractionalProfiler::FractionalHistogram::numBins * node + timing] = 100 * node + timing;   
        }
    } 

    // timing stats
    ifstream referenceFile("../../../test/writeStatsReference.proFile.csv");
    referenceFile.seekg(0, ios::end);
    size_t size = referenceFile.tellg();
    string referenceString(size, ' ');
    referenceFile.seekg(0);
    referenceFile.read(&referenceString[0], size); // In C++11 and above, string is required to have continuous storage

    output = FractionalProfiler::writeTimingsHeader(move(output));
    output = FractionalProfiler::writeTimingsStats(input, move(output), "hello_world()", dummyLabeller, 42);
    ASSERT_EQ(stream->str(), referenceString);

    // Calls per second
    stream = new ostringstream();
    output = unique_ptr<ostream>(stream);
    output = FractionalProfiler::writeCallsPerSecondsStats("hello_world()", 42, 1.337, move(output));
    ASSERT_EQ(stream->str(), "hello_world(),42,1.337000\n");
}

TEST(ProfilerTest, ProfilerRegister) {
    // An instance is already created in ParallelContext.cpp
    ProfilerRegister::createInstance("dummyFile");
    // TODO: Create a properly named file (filename determined by prefix passed on the command line)
    auto profilerRegister = ProfilerRegister::getInstance();
    ASSERT_ANY_THROW(profilerRegister->getProfiler("non-existent"));
    ASSERT_ANY_THROW(profilerRegister->getProfiler(""));

    profilerRegister->registerProfiler("hello_world()");
    ASSERT_NO_THROW(profilerRegister->getProfiler("hello_world()"));
    ASSERT_EQ(profilerRegister->getProfiler("hello_world()")->getName(), "hello_world()");
    ASSERT_ANY_THROW(profilerRegister->getProfiler("Hello_World()"));
    
    auto hwProfiler = profilerRegister->getProfiler("hello_world()");
    hwProfiler->startTimer();
    hwProfiler->endTimer();
    hwProfiler->saveTimer(111970);
    ASSERT_EQ(hwProfiler->getHistogram()->numEvents(), 1);

    // ASSERT_EQ(0, remove("profile_2.callsPerSecond.csv"));
    // ASSERT_EQ(0, remove("profile_2.proFile.csv"));
    // ASSERT_EQ(0, remove("dummyFile"));
}

TEST(ProfilerTest, PausingTimers) {
    auto profiler = FractionalProfiler("hello_world()");

    chrono::time_point<chrono::high_resolution_clock> start, end;
    start = chrono::high_resolution_clock::now();
    profiler.startTimer();
    for (int i = 0; i < 15; i++) {
        std::this_thread::sleep_for(chrono::nanoseconds((1 << 20)));
        profiler.endTimer();
        std::this_thread::sleep_for(chrono::nanoseconds((1 << 24) + (1 << 15)));
        profiler.resumeTimer();
    }
    profiler.endTimer();
    end = chrono::high_resolution_clock::now();

    profiler.saveTimer(1);
    ASSERT_FALSE(profiler.isRunning());

    auto hist = profiler.getHistogram();
    ASSERT_EQ(hist->numEvents(), 1);

    float eventsPerSecond = 1 / ((end - start).count() / pow(10, 9));
    ASSERT_NEAR(profiler.eventsPerSecond(), eventsPerSecond, eventsPerSecond * 0.01);
}

TEST(ProfilerTest, SavingTimersBasic) {
    auto profiler = FractionalProfiler("hello_world()");
    ASSERT_ANY_THROW(profiler.saveTimer(0));

    profiler = FractionalProfiler("hello_world()");
    ASSERT_ANY_THROW(profiler.discardTimer());

    profiler = FractionalProfiler("hello_world()");
    ASSERT_ANY_THROW(profiler.abortTimer());

    profiler = FractionalProfiler("hello_world()");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.discardTimer());
}

TEST(ProfilerTest, SavingTimersAdvanced) {
    auto profiler = FractionalProfiler("hello_world()");
    profiler.startTimer();
    ASSERT_TRUE(profiler.isRunning());
    ASSERT_NO_THROW(profiler.abortTimer());
    ASSERT_FALSE(profiler.isRunning());
    
    profiler = FractionalProfiler("hello_world()");
    profiler.startTimer();
    profiler.endTimer();
    profiler.resumeTimer();
    ASSERT_ANY_THROW(profiler.discardTimer());
    
    profiler = FractionalProfiler("hello_world()");
    profiler.startTimer();
    profiler.endTimer();
    ASSERT_NO_THROW(profiler.saveTimer(1));
    ASSERT_EQ(profiler.getHistogram()->numEvents(), 1);
    
    profiler.startTimer();
    ASSERT_TRUE(profiler.isRunning());
    profiler.endTimer();
    profiler.discardTimer();
    ASSERT_FALSE(profiler.isRunning());
    ASSERT_EQ(profiler.getHistogram()->numEvents(), 1);
    ASSERT_ANY_THROW(profiler.discardTimer());
}

TEST(ProfilerTest, SavingTimersOften) {
    auto profiler = FractionalProfiler("hello_word()");
    for (int i = 0; i < 10; i++) {
        profiler.startTimer();
        profiler.endTimer();
        ASSERT_FALSE(profiler.isRunning());
        ASSERT_NO_THROW(profiler.discardTimer());
        
        profiler.startTimer();
        profiler.endTimer();
        profiler.resumeTimer();
        profiler.endTimer();
        profiler.saveTimer(1);
        ASSERT_FALSE(profiler.isRunning());
    }
    ASSERT_EQ(profiler.getHistogram()->numEvents(), 10);
}

TEST(ProfilerTest, OverallStatistics) {
    //ProfilerRegister::createInstance("tmplog.log");
    auto profilerRegister = ProfilerRegister::getInstance();
    auto stats = profilerRegister->getStats();

    // In the beginning, are all values set to 0?
    ASSERT_EQ(stats->nsSumInsideMPI, 0);
    ASSERT_EQ(stats->nsSumOutsideMPI, 0);
    ASSERT_EQ(stats->nsSumWait, 0);
    ASSERT_EQ(stats->nsSumWork, 0);
    ASSERT_EQ(stats->numIterations, 0);
    ASSERT_EQ(stats->timesIWasSlowest, 0);

    // Are the values unchanged when running timers?
    // The values have to be manually changed and are not to be modified by timers.
    auto zuse = profilerRegister->registerProfiler("Zuse");
    zuse->startTimer();
    zuse->endTimer();
    zuse->saveTimer(zuse->getTimer());

    ASSERT_EQ(stats->nsSumInsideMPI, 0);
    ASSERT_EQ(stats->nsSumOutsideMPI, 0);
    ASSERT_EQ(stats->nsSumWait, 0);
    ASSERT_EQ(stats->nsSumWork, 0);
    ASSERT_EQ(stats->numIterations, 0);
    ASSERT_EQ(stats->timesIWasSlowest, 0);

    // Is the information printed correctly?
    vector<uint64_t> allStatsVec = {
        1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16
    };
    auto allStatsVecPtr = make_shared<vector<uint64_t>>(allStatsVec); // Makes a ptr to a copy!

    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return (rank == 0 ? string("Schneewittchen") : string("Zwerg"));
    };

    auto stream = new ostringstream();
    unique_ptr<ostream> output = unique_ptr<ostream>(stream);
    output = FractionalProfiler::writeOverallStats(allStatsVecPtr, dummyLabeller, move(output));

    ASSERT_EQ(stream->str(),
        "rank,processor,nsSumInsideMPI,nsSumOutsideMPI,nsSumWait,nsSumWork,numIterations,timesIWasSlowest\n0,Schneewittchen,1,2,3,4,5,6\n1,Zwerg,11,12,13,14,15,16\n"
    );
    //ASSERT_EQ(0, remove("tmplog.log"));
}