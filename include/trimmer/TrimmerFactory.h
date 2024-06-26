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
#include "trimmer/SeedClippingTrimmer.h"
#include <sstream>
#include <vector>

namespace rabbit
{
    namespace trim
    {
        class TrimmerFactory{
            public:
                TrimmerFactory();
                TrimmerFactory(rabbit::Logger& logger_, int threadId);
                ~TrimmerFactory();
                
                Trimmer* makeOneTrimmer(std::string step, rabbit::trim::RabbitTrimParam& rp);
                // void makeTrimmers(std::string steps, int phred, std::vector<Trimmer*>& trimmers);
                void makeTrimmers(rabbit::trim::RabbitTrimParam& rp, std::vector<std::string> steps, std::vector<Trimmer*>& trimmers);
                
            private:
                rabbit::Logger logger;
                int threadId;
        };
    } // namespace trim
    
} // namespace rabbit


#endif

