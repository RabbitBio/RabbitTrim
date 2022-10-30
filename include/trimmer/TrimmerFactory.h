#ifndef TRIMMER_FACTORY_H
#define TRIMMER_FACTORY_H

#include <limits.h>
#include <string>
#include "Logger.h"
#include "param.h"
#include "trimmer/Trimmer.h"
#include "trimmer/IlluminaClippingTrimmer.h"
#include "trimmer/LeadingTrimmer.h"
#include "trimmer/TrailingTrimmer.h"
#include "trimmer/HeadCropTrimmer.h"
#include "trimmer/TailCropTrimmer.h"
#include "trimmer/CropTrimmer.h"
#include "trimmer/SlidingWindowTrimmer.h"
#include "trimmer/MaximumInformationTrimmer.h"
#include "trimmer/MinLenTrimmer.h"
#include "trimmer/MaxLenTrimmer.h"
#include "trimmer/AvgQualTrimmer.h"
#include "trimmer/BaseCountTrimmer.h"
#include "trimmer/ToPhred33Trimmer.h"
#include "trimmer/ToPhred64Trimmer.h"
#include <sstream>
#include <vector>

namespace rabbit
{
    namespace trim
    {
        class TrimmerFactory{
            public:
                TrimmerFactory();
                TrimmerFactory(rabbit::Logger& logger_);
                ~TrimmerFactory();
                
                Trimmer* makeOneTrimmer(std::string step, int phred);
                void makeTrimmers(std::string steps, int phred, std::vector<Trimmer*>& trimmers);
                void makeTrimmers(std::vector<std::string> steps, int phred, std::vector<Trimmer*>& trimmers);
                
            private:
                rabbit::Logger logger;
        };
    } // namespace trim
    
} // namespace rabbit


#endif

