#ifndef TRIM_STAT_H
#define TRIM_STAT_H
#include <vector>
#include <string>
#include <fstream>
#include "Logger.h"

namespace rabbit
{
    namespace log
    {
        typedef unsigned long long uint64;
        class TrimStat
        {
            public:
                uint64 readsInput;
                uint64 readsSurvivingBoth;
                uint64 readsSurvivingForward;
                uint64 readsSurvivingReverse;

                uint64 readsDropped;
                uint64 realHit;
                uint64 tailHit;
                uint64 dimer;
                
                TrimStat();
                TrimStat(rabbit::Logger& logger);
                TrimStat(const TrimStat &trimStat_);
                ~TrimStat();
                void merge(std::vector<TrimStat>& trimStatArr);
                void printSE(std::string filename);
                void printPE(std::string filename);
                void print(std::string filename);

            private:
                rabbit::Logger logger;
        };
    } // namespace log
    
} // namespace rabbit

#endif
