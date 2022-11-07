#ifndef PARAM_H
#define PARAM_H
#include <string>
#include <vector>
#include "Logger.h"

namespace rabbit
{
    namespace trim
    {
        const int MEM_PER_CHUNK = 1 << 22; // size of chunk
        const int MEM_PER_TRIMLOG_BUFFER = 1 << 20; // size of chunk
        const std::string Illumina_adapter_r1 = "AGATCGGAAGAGC";
        const std::string Illumina_adapter_r2 = "AGATCGGAAGAGC";
        const std::string Nextera_adapter_r1 = "CTGTCTCTTATACACATCT";
        const std::string Nextera_adapter_r2 = "CTGTCTCTTATACACATCT";
        const std::string Transposase_adapter_r1 = "TCGTCGGCAGCGTC";
        const std::string Transposase_adapter_r2 = "GTCTCGTGGGCTCG";
        const std::string BGI_adapter_r1 = "AAGTCGGAGGCCAAGCGGTC";
        const std::string BGI_adapter_r2 = "AAGTCGGATCGTAGCCATGT";

        std::pair<std::string, std::string> getBuiltInAdapter(std::string seqKit){
            if(seqKit.compare("Nextera") == 0) return std::make_pair(Nextera_adapter_r1, Nextera_adapter_r2);
            if(seqKit.compare("Transposase") == 0) return std::make_pair(Transposase_adapter_r1, Transposase_adapter_r2);
            if(seqKit.compare("BGI") == 0) return std::make_pair(BGI_adapter_r1, BGI_adapter_r2);
            return std::make_pair(Illumina_adapter_r1, Illumina_adapter_r2);
        }

        struct RabbitTrimParam{
            RabbitTrimParam() = default;
            void prepare(){
                if(seqKit.size() || (seqA.size() == 0 && seqB.size() == 0)){
                    // use built-in adapter
                    std::pair<std::string, std::string> seq = getBuiltInAdapter(seqKit);
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
            
            std::string seqA;
            std::string seqB;
            std::string seqKit;
            int minLen;
            int minQual;
            int window;
        };

    } // namespace trim
} // namespace name



#endif
