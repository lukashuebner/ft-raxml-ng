#include <cstdint>
#include <map>

using namespace std;

class FractionalProfiler {
    public:
        class FractionalHistogram {
            private:
                shared_ptr<vector<uint64_t>> bins;
                int16_t getBin(float n);

            public:
                explicit FractionalHistogram();
                void event(float number);
                uint64_t& operator[](size_t idx);
                const uint64_t& operator[](size_t idx) const;
                const shared_ptr<vector<uint64_t>> data() const;
                operator std::string();
                uint64_t numEvents() const;
                static const uint16_t numBins = 201;
                static const string binName(uint16_t bin);
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
        shared_ptr<FractionalHistogram> eventCounter = make_shared<FractionalHistogram>();
        explicit FractionalProfiler(string name);
        void startTimer(bool resume = false);
        bool isRunning() const;
        void endTimer();
        void abortTimer();
        void resumeTimer();
        void discardTimer();
        uint64_t getTimer() const;
        void saveTimer(uint64_t min);
        shared_ptr<FractionalHistogram> getHistogram() const;
        const string getName() const; 
        
        static unique_ptr<ostream> writeTimingsHeader(unique_ptr<ostream> file);
        static unique_ptr<ostream> writeTimingsStats(shared_ptr<vector<uint64_t>> data, unique_ptr<ostream> file, const string& timerName,
                                                     string (*rankToProcessorName) (size_t), int secondsPassed);
        static unique_ptr<ostream> writeCallsPerSecondsHeader(unique_ptr<ostream> file);
        static unique_ptr<ostream> writeCallsPerSecondsStats(string timer, int secondsPassed, float callsPerSecond, unique_ptr<ostream> file);
        static unique_ptr<ostream> writeOverallStats(shared_ptr<vector<uint64_t>> allStatsVec, string (*rankToProcessorName) (size_t), unique_ptr<ostream> file);
        
        float eventsPerSecond() const;
        float secondsPassed() const;
        ~FractionalProfiler();
};

class ProfilerRegister {
    public:
        struct ProfilerStats {
            uint64_t nsSumWork = 0;
            uint64_t nsSumWait = 0;
            uint64_t nsSumInsideMPI = 0;
            uint64_t nsSumOutsideMPI = 0;
            uint64_t timesIWasSlowest = 0;
            uint64_t numIterations = 0;
        };

    private:
        unique_ptr<ostream> proFile = nullptr;
        unique_ptr<ostream> callsPerSecondFile = nullptr;
        shared_ptr<map<string, shared_ptr<FractionalProfiler>>> profilers = nullptr;
        shared_ptr<ProfilerStats> stats = nullptr;
        void createProFile(string path);
        
        ProfilerRegister(string logFile);
        static shared_ptr<ProfilerRegister> singleton;
        string overallStatsFilename = "";

    public:
        static shared_ptr<ProfilerRegister> getInstance();
        static shared_ptr<ProfilerRegister> createInstance(string logFile);
        ProfilerRegister() = delete;

        shared_ptr<ProfilerStats> getStats();

        shared_ptr<FractionalProfiler> registerProfiler(string name);
        shared_ptr<FractionalProfiler> getProfiler(string name) const;
        void saveProfilingData(bool master, size_t num_ranks, string (*rankToProcessorName) (size_t), MPI_Comm comm);
};