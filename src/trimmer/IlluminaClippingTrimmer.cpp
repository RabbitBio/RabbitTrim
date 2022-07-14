#include "IlluminaClippingTrimmer.h"

using namespace rabbit;


IlluminaClippingTrimmer::IlluminaClippingTrimmer(rabbit::Logger& logger_, std::string fastaAdapterFile, int seedMaxMiss_, int minPalindromeLikelihood_, int minSequenceLikelihood_, int minPrefix_ = 1, int palindromeKeepBoth_ = false) : logger(logger_){
    // adapter file 
    std::ifstream fastaAdapter(fastaAdapterFile.to_str(), std::ifstream::in);
    if(!fastaAdapter.is_open()){
        logger.errorln("Can not open fastaAdapterFile : " + fastaAdapterFile);
    }
    
    seedMaxMiss = seedMaxMiss_;
    minPalindromeLikelihood = minPalindromeLikelihood_;
    minSequenceLikelihood = minSequenceLikelihood_;
    minPrefix = minPrefix_;
    palindromeKeepBoth = palindromeKeepBoth_;
    
    minSequenceOverlap = (int)(minSequenceLikelihood / LOG10_4);
    minSequenceOverlap = minSequenceOverlap > 15 ? 15 : minSequenceOverlap; // TODO 15 是依据什么呢？
    

}
IlluminaClippingTrimmer::processOneRecord(Reference& rec){

}