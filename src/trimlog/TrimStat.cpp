#include "trimlog/TrimStat.h"
using namespace rabbit::trim;

TrimStat::TrimStat(){
    readsInput = 0ULL;
    readsSurvivingBoth = 0ULL;
    readsSurvivingForward = 0ULL;
    readsSurvivingReverse = 0ULL;
} 

TrimStat::TrimStat(rabbit::Logger& logger) : logger(logger){
    readsInput = 0ULL;
    readsSurvivingBoth = 0ULL;
    readsSurvivingForward = 0ULL;
    readsSurvivingReverse = 0ULL;
    
} 

TrimStat::TrimStat(const TrimStat& trimStat_){
  readsInput = trimStat_.readsInput;
  readsSurvivingBoth = trimStat_.readsSurvivingBoth;
  readsSurvivingForward = trimStat_.readsSurvivingForward;
  readsSurvivingReverse = trimStat_.readsSurvivingReverse;
}

TrimStat::~TrimStat() = default;

// void TrimStat::operator=(const TrimStat& ts){
//   readsInput = ts.readsInput;
//   readsSurvivingBoth = ts.readsSurvivingBoth;
//   readsSurvivingForward = ts.readsSurvivingForward;
//   readsSurvivingReverse = ts.readsSurvivingReverse;
// }

void TrimStat::merge(std::vector<TrimStat>& trimStatArr){
    for(auto t : trimStatArr){
        readsInput += t.readsInput;
        readsSurvivingBoth += t.readsSurvivingBoth;
        readsSurvivingForward += t.readsSurvivingForward;
        readsSurvivingReverse += t.readsSurvivingReverse;
    }

}

void TrimStat::printSE(std::string filename){
    const char* path = filename.c_str();
    uint64 dropped = readsInput - readsSurvivingForward;
    double survivingPercent = (100.0 * readsSurvivingForward) / readsInput;
    double droppedPercent = (100.0 * dropped) / readsInput;
    
    std::ofstream fout(path);
    if(fout.fail()) { 
          logger.errorln("\033[1;34mError: cannot open file " + filename +"!\033[0m\n");
          return;
    }
    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(2);
    fout << "Input Reads: " << readsInput << '\n'
         << "Surviving Reads: " << readsSurvivingForward << '\n'
         << "Surviving Read Percent: "  << survivingPercent << "%" << '\n'
         << "Dropped Reads: " << dropped << '\n'
         << "Dropped Read Percent: "  << droppedPercent << "%" << '\n';
    fout.close();
    
}

void TrimStat::printPE(std::string filename){
    const char* path = filename.c_str();
    uint64 dropped = readsInput - readsSurvivingBoth - readsSurvivingForward - readsSurvivingReverse;
    double survivingBothPercent  = (100.0 * readsSurvivingBoth) / readsInput;
    double survivingForwardPercent = (100.0 * readsSurvivingForward) / readsInput;
    double survivingReversePercent = (100.0 * readsSurvivingReverse) / readsInput;
    double droppedPercent = (100.0 * dropped) / readsInput;
    
    std::ofstream fout(path);
    if(fout.fail()) { 
          logger.errorln("\033[1;34mError: cannot open file " + filename +"!\033[0m\n");
          return;
    }
    fout.setf(std::ios::fixed, std::ios::floatfield);
    fout.precision(2);
    fout << "Input Read Pairs: " << readsInput << '\n'
         << "Both Surviving Reads: " << readsSurvivingBoth << '\n'
         << "Both Surviving Read Percent: " << survivingBothPercent << "%" << '\n'
         << "Forward Only Surviving Reads: " << readsSurvivingForward << '\n'
         << "Forward Only Surviving Read Percent: " << survivingForwardPercent << "%" << '\n'
         << "Reverse Only Surviving Reads: " << readsSurvivingReverse << '\n'
         << "Reverse Only Surviving Read Percent: " << survivingReversePercent << "%" << '\n'
         << "Dropped Reads: " << dropped << '\n'
         << "Dropped Read Percent: "  << droppedPercent << "%" << '\n';
    fout.close();
}
