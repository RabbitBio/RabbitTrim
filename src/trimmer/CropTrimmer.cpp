#include "CropTrimmer.h"

using namespace rabbit;

void CropTrimmer::processOneRecord(Reference& rec){
    int cur_len = rec.length;
    cur_len = cur_len > len ? len : cur_len;
    rec.length = cur_len;
}