#ifndef TRIMMER_H
#define TRIMMER_H
#include "../io/Reference.h"

namespace rabbit{
    class Trimmer{
        public:
            virtual void processOneRecord(Reference& rec) = 0;
            virtual ~Trimmer() = default;
    };
}
#endif
