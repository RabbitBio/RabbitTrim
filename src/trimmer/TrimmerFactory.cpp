#include "trimmer/TrimmerFactory.h"

using namespace rabbit::trim;

TrimmerFactory::TrimmerFactory() = default;
TrimmerFactory::TrimmerFactory(rabbit::Logger& logger_):logger(logger_){}
TrimmerFactory::~TrimmerFactory() = default;
                
Trimmer* TrimmerFactory::makeOneTrimmer(std::string step, int phred){
    std::string trimmerName;
    std::string trimmerArgs;
    std::size_t idx = step.find(":");
    trimmerName = step.substr(0,idx);
    if(idx != std::string::npos)
        trimmerArgs = step.substr(idx+1);

    if(trimmerName.compare("ILLUMINACLIP") == 0) {
        std::string fastaAdapterFile;
        int seedMaxMiss;
        int minPalindromeLikelihood;
        int minSequenceLikelihood;
        int minPrefix = 1;
        bool palindromeKeepBoth = false; 
        
        std::size_t pos = trimmerArgs.find(":");
        fastaAdapterFile = trimmerArgs.substr(0, pos);
        std::string::size_type sz;
        std::string::size_type cur_pos = pos + 1;

        seedMaxMiss = std::stoi(trimmerArgs.substr(cur_pos), &sz);
        cur_pos += sz + 1;
        minPalindromeLikelihood = std::stoi(trimmerArgs.substr(cur_pos), &sz);
        cur_pos += sz + 1;
        minSequenceLikelihood = std::stoi(trimmerArgs.substr(cur_pos), &sz);
        cur_pos += sz + 1;
        
        if(trimmerArgs.size() > cur_pos){
            minPrefix = std::stoi(trimmerArgs.substr(cur_pos), &sz);
            cur_pos += sz + 1;
            if(trimmerArgs.size() > cur_pos){
                std::string keepBothStr = trimmerArgs.substr(cur_pos);
                if(keepBothStr.compare("true") == 0) palindromeKeepBoth = true;
            }
        }
        
        return new IlluminaClippingTrimmer(logger, phred, fastaAdapterFile, seedMaxMiss, minPalindromeLikelihood, minSequenceLikelihood, minPrefix, palindromeKeepBoth);
    }
    if(trimmerName.compare("LEADING") == 0) {
        int qual = std::stoi(trimmerArgs, nullptr);
        return new LeadingTrimmer(qual, phred);
    }
    if(trimmerName.compare("TRAILING") == 0) {
        int qual = std::stoi(trimmerArgs, nullptr);
        return new TrailingTrimmer(qual, phred);
    }
    if(trimmerName.compare("HEADCROP") == 0) {
        std::size_t pos = trimmerArgs.find(":");
        int bases = std::stoi(trimmerArgs, nullptr);
        int maxLength = 1 << 30;
        if(pos != std::string::npos){
            std::string maxLength_str = trimmerArgs.substr(pos+1);
            maxLength = std::stoi(maxLength_str, nullptr);
        }
        return new HeadCropTrimmer(bases, maxLength);
    }
    if(trimmerName.compare("TAILCROP") == 0) {
        std::size_t pos = trimmerArgs.find(":");
        int bases = std::stoi(trimmerArgs, nullptr);
        int maxLength = 1 << 30;
        if(pos != std::string::npos){
            std::string maxLength_str = trimmerArgs.substr(pos+1);
            maxLength = std::stoi(maxLength_str, nullptr);
        }
        return new TailCropTrimmer(bases, maxLength);

    }
    if(trimmerName.compare("CROP") == 0) {
        int len = std::stoi(trimmerArgs,nullptr);
        return new CropTrimmer(len);
    }
    if(trimmerName.compare("SLIDINGWINDOW") == 0) {
        std::size_t pos = trimmerArgs.find(":");
        int windowLength = std::stoi(trimmerArgs, nullptr);
        std::string requiredQuality_str = trimmerArgs.substr(pos+1);
        int requiredQuality = std::stoi(requiredQuality_str, nullptr);
        return new SlidingWindowTrimmer(windowLength, requiredQuality, phred);
    }
    if(trimmerName.compare("MAXINFO") == 0) {
        std::size_t pos = trimmerArgs.find(":");
        int parLength = std::stoi(trimmerArgs, nullptr);
        std::string strictness_str = trimmerArgs.substr(pos+1);
        float strictness = std::stof(strictness_str, nullptr);
        return new MaximumInformationTrimmer(parLength, strictness);
    }
    if(trimmerName.compare("MINLEN") == 0) {
        int minLen = std::stoi(trimmerArgs, nullptr) ;
        return new MinLenTrimmer(minLen);
    }
    if(trimmerName.compare("MAXLEN") == 0) {
        int maxLen = std::stoi(trimmerArgs, nullptr) ;
        return new MaxLenTrimmer(maxLen);
    }
    if(trimmerName.compare("AVGQUAL") == 0) {
        int qual = std::stoi(trimmerArgs, nullptr);
        return new AvgQualTrimmer(qual, phred);
    }
    if(trimmerName.compare("BASECOUNT") == 0) {
        std::size_t pos = trimmerArgs.find(":");
        std::string bases = trimmerArgs.substr(0, pos);
        if(pos != std::string::npos){
            std::string::size_type sz;
            std::string minCount_str = trimmerArgs.substr(pos + 1);
            int minCount  = std::stoi(minCount_str, &sz);
            int maxCount = 1 << 30;
            if(sz < minCount_str.size()) maxCount = std::stoi(minCount_str.substr(1), nullptr);
            return new BaseCountTrimmer(bases, minCount, maxCount);
        }
    }
    if(trimmerName.compare("TOPHRED33") == 0) {
        return new ToPhred33Trimmer(phred) ;
    }
    if(trimmerName.compare("TOPHRED64") == 0) {
        return new ToPhred64Trimmer(phred) ;
    }
    
    logger.errorln("Unknown trimmer: " + trimmerName);
    exit(1);
}

void TrimmerFactory::makeTrimmers(std::string steps, int phred, std::vector<Trimmer*>& trimmers){
    std::istringstream is(steps);
    while(is){
        std::string oneStep;
        is >> oneStep;
        if(oneStep.size()){
            Trimmer* trimmer = makeOneTrimmer(oneStep, phred);
            trimmers.push_back(trimmer);
        }
    }
}

void TrimmerFactory::makeTrimmers(std::vector<std::string> steps, int phred, std::vector<Trimmer*>& trimmers){
    for(auto oneStep : steps){
        Trimmer* trimmer = makeOneTrimmer(oneStep, phred);
        trimmers.push_back(trimmer);
    }
}
