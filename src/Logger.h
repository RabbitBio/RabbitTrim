#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>

namespace rabbit{
    class Logger{
        private:
            bool showError, showWarning, showInfo;

        public:
            Logger(bool showError_, bool showWarning_, bool showInfo_){
                showError = showError_;
                showWarning = showWarning_;
                showInfo = showInfo_;
            }

            Logger(Logger& Logger_){
                showError = logger_.showError;
                showWarning = logger_.showWarning;
                showInfo = logger_.showInfo;
            }

            void errorln(std::string message){
                if(showError){
                    std::cerr << "ERROR: " << message << std::endl;
                }
            }

            void error(std::string message){
                if(showError){
                    std::cerr << "ERROR: " << message;
                }
            }

            void warningln(std::string message){
                if(showWarning){
                    std::cerr << "WARNING: " << message << std::endl;
                }
            }

            void warning(std::string message){
                if(showWarning){
                    std::cerr << "WARNING: " << message;
                }
            }

            void infoln(std::string message){
                if(showInfo){
                    std::cerr << "INFO: " << message << std::endl;
                }
            }

            void info(std::string message){
                if(showError){
                    std::cerr << "INFO: " << message;
                }
            }
        
    };
}

#endif