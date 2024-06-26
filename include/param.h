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

        // pragzip coefficients
        const double PRAGZIP_COEF_0 = 157.8302;
        const double PRAGZIP_COEF_1 = -3.949971;
        const double PRAGZIP_COEF_2 = 0.413682;
        const double PRAGZIP_COEF_3 = -0.0147582;

        const double PRAGZIP_COEF_4 = 248.915;
        const double PRAGZIP_COEF_5 = -7.99773;
        const double PRAGZIP_COEF_6 = 0.075041;

        // pigz coefficients
        const double PIGZ_COEF_0 = 59.2099;
        const double PIGZ_COEF_1 = -0.5408;
        const double PIGZ_COEF_2 = 0.0128;

        const double PIGZ_COEF_3 = -18.3622;
        const double PIGZ_COEF_4 = 7.06845;
        const double PIGZ_COEF_5 = -0.21403;
        const double PIGZ_COEF_6 = 0.00184882;
        
        // preset of speed
        const double INPUT_READ_SPEED = 2000.0;
        const double IGZIP_DECOMPRESS_SPEED = 800.0;
        const double IGZIP_COMPRESS_SPEED = 490.0;
        const double RTRIM_PE_SPEED_PER_SEC = 500.0;
        const double RTRIM_SE_SPEED_PER_SEC = 260.0;
        const double PIGZ_SPEED_PER_SEC = 50.0;
        const double PRAGZIP_SPEED_PER_SEC = 180.0;
        const int MIN_PIGZ_THREAD_NUM = 10;
        const int MIN_PRAGZIP_THREAD_NUM = 4;

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
            bool usePragzip;
            bool usePigz;
            int pigzThreadsNum;
            int pragzipThreadsNum;
            int compressLevel;
            
            std::string seqA;
            std::string seqB;
            std::string seqKit;
            int minLen;
            int minQual;
            int window;
            double mismatch;
            bool use_default_mismatch;

            int workerThreadNum;
        };

    } // namespace trim
} // namespace name



#endif
