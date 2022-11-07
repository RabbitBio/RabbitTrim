#include "trimmer/LowQualityTrimmer.h"

using namespace rabbit::trim;

void LowQualityTrimmer::processOneRecord(Reference& rec){}
void LowQualityTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse)
{
    for(auto& rec : recs){
        processOneRecord(rec);
    }
}

void LowQualityTrimmer::processOneRecord(neoReference& rec)
{

}
void LowQualityTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse)
{
    for(auto& rec : recs){
        processOneRecord(rec);
    }
}