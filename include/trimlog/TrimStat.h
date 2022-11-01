#ifndef TRIM_STAT_H
#define TRIM_STAT_H
#include <vector>
#include <string>
#include <fstream>
#include "Logger.h"

namespace rabbit
{
    namespace trim
    {
        typedef unsigned long long uint64;
        class TrimStat
        {
            public:
                uint64 readsInput;
                uint64 readsSurvivingBoth;
                uint64 readsSurvivingForward;
                uint64 readsSurvivingReverse;

                TrimStat();
                TrimStat(rabbit::Logger& logger);
                TrimStat(const TrimStat &trimStat_);
                ~TrimStat();
                void merge(std::vector<TrimStat>& trimStatArr);
                void printSE(std::string filename);
                void printPE(std::string filename);

            private:
                rabbit::Logger logger;
        };
    } // namespace trim
    
} // namespace rabbit

#endif
