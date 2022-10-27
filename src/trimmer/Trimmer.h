#ifndef TRIMMER_H
#define TRIMMER_H
#include "../io/Reference.h"
#include <vector>

namespace rabbit{
    class Trimmer{
        public:
            virtual void processOneRecord(Reference& rec) = 0;
            virtual void processRecords(std::vector<Reference&> recs, bool isPair, bool isReverse);
            virtual ~Trimmer() = default; 
    };
}
#endif
