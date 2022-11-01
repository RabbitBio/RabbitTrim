#ifndef TRIMMER_H
#define TRIMMER_H
#include "Reference.h"
#include <vector>

namespace rabbit{
    namespace trim
    {
        class Trimmer{
            public:
                virtual void processOneRecord(neoReference& rec) = 0;
                virtual void processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse) = 0;
                virtual ~Trimmer() = default;
        };
    } // namespace trim
}
#endif
