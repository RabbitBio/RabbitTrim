#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>

namespace rabbit{
    class Logger{
        private:
            bool showError, showWarning, showInfo;

        public:
            Logger(){} 
            Logger(bool showError_, bool showWarning_, bool queit){
                showError = showError_;
                showWarning = showWarning_;
                showInfo = !queit;
            }

            Logger(const Logger& logger_){
                showError = logger_.showError;
                showWarning = logger_.showWarning;
                showInfo = logger_.showInfo;
            }

            void errorln(std::string message){
                if(showError){
                    std::cerr << "\033[1;31mERROR: " << message << "\033[0m"<< std::endl;
                }
            }

            void error(std::string message){
                if(showError){
                    std::cerr << "\033[1;31mERROR: " << message << "\033[0m";
                }
            }

            void warningln(std::string message){
                if(showWarning){
                    std::cerr << "\033[1;33mWARNING: " << message << "\033[0m" << std::endl;
                }
            }

            void warning(std::string message){
                if(showWarning){
                    std::cerr << "\033[1;33mWARNING: " << message << "\033[0m";
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
