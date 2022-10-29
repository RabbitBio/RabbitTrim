#ifndef TRIMMER_FACTORY_H
#define TRIMMER_FACTORY_H

#include <limits.h>
#include <string>
#include "../Logger.h"
#include "param.h"
#include "Trimmer.h"
#include "IlluminaClippingTrimmer.h"
#include "LeadingTrimmer.h"
#include "TrailingTrimmer.h"
#include "HeadCropTrimmer.h"
#include "TailCropTrimmer.h"
#include "CropTrimmer.h"
#include "SlidingWindowTrimmer.h"
#include "MaximumInformationTrimmer.h"
#include "MinLenTrimmer.h"
#include "MaxLenTrimmer.h"
#include "AvgQualTrimmer.h"
#include "BaseCountTrimmer.h"
#include "ToPhred33Trimmer.h"
#include "ToPhred64Trimmer.h"
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

