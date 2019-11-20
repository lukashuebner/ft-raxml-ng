#include <chrono>
#include <thread>
#include <iostream> 
#include <cstdio>

#include "RaxmlTest.hpp"
#include "../../src/Profiler.hpp"

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

    ASSERT_EQ((string)hist, "1,2,0,0,0,0,0,0,0,0,1337,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");

    auto data = hist.data();
    string res = "";
    for (auto& bin: *data) {
        res += to_string(bin) + ",";
    }
    res.pop_back();
    ASSERT_EQ((string)hist, res);

    ASSERT_EQ(hist.numEvents(), 1337 + 2 +1);
}

TEST(ProfilerTest, InvalidArguments) {
    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return to_string(*(rank == 0 ? "Schneewittchen" : "Zwerg"));
    };
    unique_ptr<ostream> output = unique_ptr<ostream>(new ostringstream());

    auto profiler = LogBinningProfiler("Testing");
    ASSERT_ANY_THROW(output = profiler.writeStats(nullptr, move(output), "hello_world()", dummyLabeller, 0));
    ASSERT_ANY_THROW(profiler.writeStats(profiler.getHistogram()->data(), nullptr, "bye_world()", dummyLabeller, 0));
    ASSERT_ANY_THROW(output = profiler.writeStats(profiler.getHistogram()->data(), move(output), "bye_world()", dummyLabeller, -1));

    ASSERT_ANY_THROW(LogBinningProfiler::writeCallsPerSecondsHeader(nullptr));
    ASSERT_ANY_THROW(LogBinningProfiler::writeCallsPerSecondsStats("hello_world", -1, 1, move(output)));
    output = unique_ptr<ostream>(new ostringstream());
    ASSERT_ANY_THROW(LogBinningProfiler::writeCallsPerSecondsStats("hello_world", 1, -1, move(output)));
    output = unique_ptr<ostream>(new ostringstream());
    ASSERT_ANY_THROW(LogBinningProfiler::writeCallsPerSecondsStats("hello_world", -1, 1, nullptr));

    ASSERT_NO_THROW(output = LogBinningProfiler::writeCallsPerSecondsStats("", 0, 0, move(output)));
}

TEST(ProfilerTest, TimerBasics) {
    // buildup
    auto profiler = LogBinningProfiler("Testing");

    // tests
    profiler.getName() == "Testing";

    auto hist = profiler.getHistogram();
    for (int i = 0; i < 64; i++) {
        ASSERT_EQ((*hist)[i], 0);
    }

    chrono::time_point<chrono::high_resolution_clock> start, end;
    start = chrono::high_resolution_clock::now();
    const unsigned int NUM_EVENTS = 100;
    for (unsigned int i = 0; i < NUM_EVENTS; i++) {
        profiler.startTimer();
        std::this_thread::sleep_for(chrono::nanoseconds((1 << 20) + (1 << 8)));
        profiler.endTimer();
        profiler.saveTimer(0);
    }
    end = chrono::high_resolution_clock::now();
    float eventsPerSecond = NUM_EVENTS / ((end - start).count() / pow(10, 9));
    hist = profiler.getHistogram();

    for (int i = 0; i < 64; i++) {
        if (i == 20) {
            ASSERT_GE((*hist)[i], 70);
        } else {
            ASSERT_LE((*hist)[i], 7);
        }
    }

    ASSERT_EQ(hist->numEvents(), 100);
    ASSERT_GT(profiler.eventsPerSecond(), 100);
    ASSERT_NEAR(profiler.eventsPerSecond(), eventsPerSecond, eventsPerSecond * 0.01);
    ASSERT_EQ(profiler.eventsPerSecond(), profiler.eventsPerSecond());
}

TEST(ProfilerTest, InvalidStates) {
    auto profiler = LogBinningProfiler("Testing");
    // Try to stop the timer before it started
    ASSERT_ANY_THROW(profiler.endTimer());
    // After that, most other methods should throw an exception because of invalid state, too.
    ASSERT_ANY_THROW(profiler.startTimer());
    ASSERT_ANY_THROW(profiler.abortTimer());
    ASSERT_ANY_THROW(profiler.eventsPerSecond());
    ASSERT_ANY_THROW(profiler.getHistogram());

    // Try to start the timer twice in a row
    profiler = LogBinningProfiler("Testing");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.startTimer());

    // Try to abort a non-running timer
    profiler = LogBinningProfiler("Testing");
    ASSERT_ANY_THROW(profiler.abortTimer());

    // Try to end an aborted timer
    profiler = LogBinningProfiler("Testing");
    profiler.startTimer();
    profiler.abortTimer();
    ASSERT_ANY_THROW(profiler.endTimer());

    // Try to save a running timer and an invalid timer
    profiler = LogBinningProfiler("Testing");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.saveTimer(0));

    profiler = LogBinningProfiler("Testing");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.startTimer());
    ASSERT_ANY_THROW(profiler.saveTimer(0));

    // Try to discard an invalid timer
    ASSERT_ANY_THROW(profiler.discardTimer());
}

TEST(ProfilerTest, TestIsRunning) {
    auto profiler = LogBinningProfiler("Testing");
    
    ASSERT_FALSE(profiler.isRunning());
    profiler.startTimer();
    ASSERT_TRUE(profiler.isRunning());
    profiler.endTimer();
    ASSERT_FALSE(profiler.isRunning());
    profiler.saveTimer(0);
    ASSERT_FALSE(profiler.isRunning());

    profiler.startTimer();
    profiler.abortTimer();
    ASSERT_FALSE(profiler.isRunning());

    ASSERT_ANY_THROW(profiler.endTimer());
    ASSERT_ANY_THROW(profiler.isRunning());
}

TEST(ProfilerTest, TimerCanBeAborted) {
    auto profiler = LogBinningProfiler("Testing");

    profiler.startTimer();
    ASSERT_NO_THROW(profiler.endTimer());
    profiler.saveTimer(0);
    profiler.startTimer();
    ASSERT_NO_THROW(profiler.abortTimer());
    ASSERT_EQ(profiler.getHistogram()->numEvents(), 1);
    ASSERT_NO_THROW(profiler.startTimer());
    profiler.abortTimer();
    ASSERT_ANY_THROW(profiler.endTimer());

    profiler = LogBinningProfiler("Testing");
    profiler.startTimer();
    profiler.endTimer();
    ASSERT_ANY_THROW(profiler.abortTimer());
}

const string statsHeader = "rank,processor,timer,secondsPassed,\"[2^00,2^01) ns\",\"[2^01,2^02) ns\",\"[2^02,2^03) ns\",\"[2^03,2^04) ns\",\"[2^04,2^05) ns\",\"[2^05,2^06) ns\",\"[2^06,2^07) ns\",\"[2^07,2^08) ns\",\"[2^08,2^09) ns\",\"[2^09,2^10) ns\",\"[2^10,2^11) ns\",\"[2^11,2^12) ns\",\"[2^12,2^13) ns\",\"[2^13,2^14) ns\",\"[2^14,2^15) ns\",\"[2^15,2^16) ns\",\"[2^16,2^17) ns\",\"[2^17,2^18) ns\",\"[2^18,2^19) ns\",\"[2^19,2^20) ns\",\"[2^20,2^21) ns\",\"[2^21,2^22) ns\",\"[2^22,2^23) ns\",\"[2^23,2^24) ns\",\"[2^24,2^25) ns\",\"[2^25,2^26) ns\",\"[2^26,2^27) ns\",\"[2^27,2^28) ns\",\"[2^28,2^29) ns\",\"[2^29,2^30) ns\",\"[2^30,2^31) ns\",\"[2^31,2^32) ns\",\"[2^32,2^33) ns\",\"[2^33,2^34) ns\",\"[2^34,2^35) ns\",\"[2^35,2^36) ns\",\"[2^36,2^37) ns\",\"[2^37,2^38) ns\",\"[2^38,2^39) ns\",\"[2^39,2^40) ns\",\"[2^40,2^41) ns\",\"[2^41,2^42) ns\",\"[2^42,2^43) ns\",\"[2^43,2^44) ns\",\"[2^44,2^45) ns\",\"[2^45,2^46) ns\",\"[2^46,2^47) ns\",\"[2^47,2^48) ns\",\"[2^48,2^49) ns\",\"[2^49,2^50) ns\",\"[2^50,2^51) ns\",\"[2^51,2^52) ns\",\"[2^52,2^53) ns\",\"[2^53,2^54) ns\",\"[2^54,2^55) ns\",\"[2^55,2^56) ns\",\"[2^56,2^57) ns\",\"[2^57,2^58) ns\",\"[2^58,2^59) ns\",\"[2^59,2^60) ns\",\"[2^60,2^61) ns\",\"[2^61,2^62) ns\",\"[2^62,2^63) ns\",\"[2^63,2^64) ns\"\n";
const string callsPerSecondsHeader = "timer,secondsPassed,callsPerSecond\n";
TEST(ProfilerTest, writeStats_Empty) {
    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return to_string(*(rank == 0 ? "Schneewittchen" : "Zwerg"));
    };
    auto stream = new ostringstream();
    unique_ptr<ostream> output = unique_ptr<ostream>(stream);
    auto input = make_shared<vector<uint64_t>>();

    output = LogBinningProfiler::writeStats(input, move(output), "hello_world()", dummyLabeller, 0);
    ASSERT_TRUE(output);
    ASSERT_EQ(stream->str(), "");

    output = LogBinningProfiler::writeStatsHeader(move(output));
    ASSERT_TRUE(output);
    ASSERT_EQ(stream->str(), statsHeader);

    stream = new ostringstream();
    output = unique_ptr<ostream>(stream);
    output = LogBinningProfiler::writeCallsPerSecondsHeader(move(output));
    ASSERT_TRUE(output);
    ASSERT_EQ(stream->str(), callsPerSecondsHeader);
}

TEST(ProfilerTest, writeStats_Basic) {
    string (*dummyLabeller) (size_t) = [](size_t rank) {
        return (rank == 0 ? string("Schneewittchen") : string("Zwerg"));
    };

    auto stream = new ostringstream();
    unique_ptr<ostream> output = unique_ptr<ostream>(stream);
    auto input = make_shared<vector<uint64_t>>(64 * 4);

    for (int node = 0; node < 4; node++) {
        for (int timing = 0; timing < 64; timing++) {
            (*input)[64 * node + timing] = 100 * node + timing;   
        }
    } 

    // timing stats
    output = LogBinningProfiler::writeStatsHeader(move(output));
    output = LogBinningProfiler::writeStats(input, move(output), "hello_world()", dummyLabeller, 42);
    ASSERT_EQ(stream->str(),
        statsHeader + 
        "0,Schneewittchen,hello_world(),42,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63\n"
        "1,Zwerg,hello_world(),42,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163\n" +
        "2,Zwerg,hello_world(),42,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263\n" +
        "3,Zwerg,hello_world(),42,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363\n"
    );

    // Calls per second
    stream = new ostringstream();
    output = unique_ptr<ostream>(stream);
    output = LogBinningProfiler::writeCallsPerSecondsStats("hello_world()", 42, 1.337, move(output));
    ASSERT_EQ(stream->str(), "hello_world(),42,1.337000\n");
}

TEST(ProfilerTest, ProfilerRegister) {
    auto profilerRegister = ProfilerRegister::createInstance("dummyFile");
    ASSERT_ANY_THROW(profilerRegister->getProfiler("non-existent"));
    ASSERT_ANY_THROW(profilerRegister->getProfiler(""));

    profilerRegister->registerProfiler("hello_world()");
    ASSERT_NO_THROW(profilerRegister->getProfiler("hello_world()"));
    ASSERT_EQ(profilerRegister->getProfiler("hello_world()")->getName(), "hello_world()");
    ASSERT_ANY_THROW(profilerRegister->getProfiler("Hello_World()"));
    
    auto hwProfiler = profilerRegister->getProfiler("hello_world()");
    hwProfiler->startTimer();
    hwProfiler->endTimer();
    hwProfiler->saveTimer(0);
    ASSERT_EQ(hwProfiler->getHistogram()->numEvents(), 1);

    ASSERT_EQ(0, remove("dummyFile_callsPerSecond.csv"));
    ASSERT_EQ(0, remove("dummyFile_proFile.csv"));
}

TEST(ProfilerTest, PausingTimers) {
    auto profiler = LogBinningProfiler("hello_world()");

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

    profiler.saveTimer(0);
    ASSERT_FALSE(profiler.isRunning());

    float eventsPerSecond = 1 / ((end - start).count() / pow(10, 9));
    auto hist = profiler.getHistogram();

    for (int i = 0; i < 64; i++) {
        if (i == 24) {
            ASSERT_EQ((*hist)[i], 1);
        } else {
            ASSERT_EQ((*hist)[i], 0);
        }
    }

    ASSERT_EQ(hist->numEvents(), 1);
    ASSERT_NEAR(profiler.eventsPerSecond(), eventsPerSecond, eventsPerSecond * 0.01);
}

TEST(ProfilerTest, SavingTimersBasic) {
    auto profiler = LogBinningProfiler("hello_world()");
    ASSERT_ANY_THROW(profiler.saveTimer(0));

    profiler = LogBinningProfiler("hello_world()");
    ASSERT_ANY_THROW(profiler.discardTimer());

    profiler = LogBinningProfiler("hello_world()");
    ASSERT_ANY_THROW(profiler.abortTimer());

    profiler = LogBinningProfiler("hello_world()");
    profiler.startTimer();
    ASSERT_ANY_THROW(profiler.discardTimer());
}

TEST(ProfilerTest, SavingTimersAdvanced) {
    auto profiler = LogBinningProfiler("hello_world()");
    profiler.startTimer();
    ASSERT_TRUE(profiler.isRunning());
    ASSERT_NO_THROW(profiler.abortTimer());
    ASSERT_FALSE(profiler.isRunning());
    
    profiler = LogBinningProfiler("hello_world()");
    profiler.startTimer();
    profiler.endTimer();
    profiler.resumeTimer();
    ASSERT_ANY_THROW(profiler.discardTimer());
    
    profiler = LogBinningProfiler("hello_world()");
    profiler.startTimer();
    profiler.endTimer();
    ASSERT_NO_THROW(profiler.saveTimer(100));
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
    auto profiler = LogBinningProfiler("hello_word()");
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