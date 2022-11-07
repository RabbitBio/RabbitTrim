#ifndef TRIM_STAT_2_H
#define TRIM_STAT_2_H
namespace rabbit
{
    namespace log
    {
        typedef unsigned long long uint64; 
        class TrimStat2
        {
        public:
            uint64 readsInput;
            uint64 dropped;
            uint64 realAdapter;
            uint64 tailHit;
            
            TrimStat2();
            ~TrimStat2();
        };
    } // namespace log
} // namespace rabbit

#endif
