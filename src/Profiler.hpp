#include <cstdint>
#include <map>

using namespace std;

class LogBinningProfiler {
    public:
        class LogarithmicHistogram {
            private:
                shared_ptr<vector<uint64_t>> bins;
                int16_t log2i(uint64_t n);

            public:
                explicit LogarithmicHistogram();
                void event(uint64_t number);
                uint64_t& operator[](size_t idx);
                const uint64_t& operator[](size_t idx) const;
                const shared_ptr<vector<uint64_t>> data() const;
                operator std::string();
                uint64_t numEvents() const;
                static const uint16_t numBins = 65;
        };

    private:
        string name;
        bool running = false;
        bool savedOrDiscarded = true;
        bool invalid = false;
        chrono::time_point<chrono::high_resolution_clock> start, end;
        uint64_t nsPassed = 0;
        chrono::time_point<chrono::high_resolution_clock> firstStart = chrono::time_point<chrono::high_resolution_clock>::max();
        chrono::time_point<chrono::high_resolution_clock> lastEnd;
        uint64_t timesCalled() const;

    public:
        shared_ptr<LogarithmicHistogram> eventCounter = make_shared<LogarithmicHistogram>();
        explicit LogBinningProfiler(string name);
        void startTimer(bool resume = false);
        bool isRunning() const;
        void endTimer();
        void abortTimer();
        void resumeTimer();
        void discardTimer();
        uint64_t getTimer() const;
        void saveTimer(uint64_t min);
        shared_ptr<LogarithmicHistogram> getHistogram() const;
        const string getName() const; 
        static unique_ptr<ostream> writeStatsHeader(unique_ptr<ostream> file);
        static unique_ptr<ostream> writeStats(shared_ptr<vector<uint64_t>> data, unique_ptr<ostream> file, const string& timerName,
                                              string (*rankToProcessorName) (size_t), int secondsPassed);
        static unique_ptr<ostream> writeCallsPerSecondsHeader(unique_ptr<ostream> file);
        static unique_ptr<ostream> writeCallsPerSecondsStats(string timer, int secondsPassed, float callsPerSecond, unique_ptr<ostream> file);
        float eventsPerSecond() const;
        float secondsPassed() const;
        ~LogBinningProfiler();
};

class ProfilerRegister {
    private:
        unique_ptr<ostream> proFile = nullptr;
        unique_ptr<ostream> callsPerSecondFile = nullptr;
        shared_ptr<map<string, shared_ptr<LogBinningProfiler>>> profilers = nullptr;
        void createProFile(string path);
        
        ProfilerRegister(string logFile);
        static shared_ptr<ProfilerRegister> singleton;

    public:
        static shared_ptr<ProfilerRegister> getInstance();
        static shared_ptr<ProfilerRegister> createInstance(string logFile);
        ProfilerRegister() = delete;

        shared_ptr<LogBinningProfiler> registerProfiler(string name);
        shared_ptr<LogBinningProfiler> getProfiler(string name) const;
        void saveProfilingData(bool master, size_t num_ranks, string (*rankToProcessorName) (size_t), MPI_Comm comm);
};