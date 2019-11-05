#include <cstdint>

using namespace std;

class LogBinningProfiler {
    public:
        class LogarithmicHistogram {
            private:
                const uint8_t numBins;
                shared_ptr<vector<uint64_t>> bins;
                uint8_t log2i(uint64_t n);

            public:
                explicit LogarithmicHistogram();
                void event(uint64_t number);
                uint64_t& operator[](size_t idx);
                const uint64_t& operator[](size_t idx) const;
        };

    private:
        string name;
        bool running = false;
        bool invalid = false;
        chrono::time_point<chrono::high_resolution_clock> start, end;

    public:
        shared_ptr<LogarithmicHistogram> eventCounter = make_shared<LogarithmicHistogram>();
        explicit LogBinningProfiler(string name);
        void startTimer();
        void endTimer();
        const shared_ptr<LogarithmicHistogram> getHistogram() const;
        const string getName() const;
};