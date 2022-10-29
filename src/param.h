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
        struct RabbitTrimParam{
            RabbitTrimParam() = default;
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
        };

        
        void exec_subcommand_trimmomatic(std::string mode, RabbitTrimParam& rp, Logger& logger){
            if(mode.compare("PE") == 0){
                
            }
            else{
                if(mode.compare("SE") == 0){

                }
                else{
                    logger.errorln("The value of mode int Trimmomatic must be either PE or SE !");
                }
            }

        }
        void exec_subcommand_ktrim(RabbitTrimParam& rp, Logger& logger){

        }
    } // namespace trim
} // namespace name



#endif
