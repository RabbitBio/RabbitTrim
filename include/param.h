#ifndef PARAM_H
#define PARAM_H
#include <string>
#include <vector>
#include "Logger.h"
#include "util.h"

namespace rabbit
{
    namespace trim
    {
        const int MEM_PER_CHUNK = 1 << 22; // size of chunk
        const int MEM_PER_TRIMLOG_BUFFER = 1 << 20; // size of chunk
        const int MAX_READ_LENGTH = 512;
        const int MAX_ADAPTER_LENGTH = 128;
        const int PREREAD_COUNT = 10000;

        struct RabbitTrimParam{
            RabbitTrimParam() = default;
            void prepare(){
                // set steps && stats
                steps.push_back("SEED");
                stats = output.substr(0, output.find(".gz")) + ".trim.log";
              
                // set seqKit of Ktrim
                if(seqKit.size() || (seqA.size() == 0 && seqB.size() == 0)){
                    // use built-in adapter
                    std::pair<std::string, std::string> seq = rabbit::trim::util::getBuiltInAdapter(seqKit);
                    seqA = seq.first;
                    seqB = seq.second;
                }else{
                    if(seqA.size() == 0) seqA = seqB;
                    if(seqB.size() == 0) seqB = seqA;
                    int lenA = seqA.size();
                    int lenB = seqB.size();
                    int len = std::min(lenA, lenB);
                    seqA.substr(0,len);
                    seqB.substr(0,len);
                }

                std::cout << "Using seqA: " << seqA << '\n'
                          << "      seqB: " << seqB << '\n';
                
                
            }
            std::vector<std::string> forwardFiles;
            std::vector<std::string> reverseFiles;
            std::string output;
            std::string trimLog;
            std::string stats;
            // std::string steps;
            std::vector<std::string> steps;
            int threads;
            int phred;
            bool validatePairing;
            bool useIgzip;
            bool usePigz;
            int pigzThreadsNum;
            int compressLevel;
            
            std::string seqA;
            std::string seqB;
            std::string seqKit;
            int minLen;
            int minQual;
            int window;
            double mismatch;
            bool use_default_mismatch;
        };

    } // namespace trim
} // namespace name



#endif
