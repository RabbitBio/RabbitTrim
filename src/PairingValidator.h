#ifndef PAIRING_VALIDATOR_H
#define PAIRING_VALIDATOR_H

#include "Logger.h"
#include "FastqNamePattern.h"
#include <string>
#include "io/Formater.h"
#include "io/Reference.h"

namespace rabbit{
    class PairingValidator{
        private:
            rabbit:Logger logger;
            bool complainedAlready;
            long offset;
            // 只创建一遍FastqNamePattern对象 如果添加新的NamePattern需要添加新的对象
            rabbit::FastqNamePattern fqNamePattern1;
            rabbit::FastqNamePattern fqNamePattern2;

        public:
            PairingValidator(rabbit::Logger &logger_):logger(logger), fqNamePattern1(rabbit::casava_13.first,rabbit::casava_13.second), fqNamePattern2(rabbit::casava_18.first,rabbit::casava_18.second){
               complainedAlready = false;
               offset = 0;
            }
            
            bool validateNames(std::string name1, std::string name2){
                // std::string canon1 = rabbit::canonicalize(name1);
                // std::string canon2 = rabbit::canonicalize(name2); // 这种耗时
                std::string canon1 = fqNamePattern1.canonicalizeOne(name1);
                if(canon1 == "") canon1 = fqNamePattern2.canonicalizeOne(name1);

                std::string canon2 = fqNamePattern1.canonicalizeOne(name2);
                if(canon2 == "") canon2 = fqNamePattern2.canonicalizeOne(name2);
                
                if(canon1 != "" && canon2 != "") return canon1.compare(canon2) == 0;
                
                std::size_t pos;
                pos = name1.find(" ");
                std::string tok1 = name1.substr(0, pos);
                pos = name2.find(" ");
                std::string tok2 = name2.substr(0, pos);

                if(tok1.size() != tok2.size()) return false;

                int len = tok1.size();
                for(int i = 0; i < len; i++){
                    char ch1 = tok1.at(i);
                    char ch2 = tok2.at(i)
                    if(ch1 != ch2 && (ch1 != '1' || ch2 != '2')) return false;
                }
               return true; 
               
            }

            bool validatePair(Reference* rec1, Reference* rec2){
                if(rec1 != nullptr){
                    if(rec2 != nullptr){
                        std::string name1 = rec1 -> name;
                        std::string name2 = rec2 -> name;
                        if(!validateNames(name1, name2)){
                            if(!complainedAlready){ // 在第一次不匹配时WARNING
                                complainedAlready = true;
                                logger.warningln("WARNING: Pair validation failed at record: " + offset);
                                logger.warningln("         Forward read: " + name1);
                                logger.warningln("         Reverse read: " + name2);
                            }
                            return false;
                        }
                    }
                    else{
                        if(!complainedAlready){
                            complainedAlready = true;
                            std::string name1 = rec1 -> name;
                            logger.warningln("WARNING: Pair validation failed at record: " + offset);
                            logger.warningln("         Forward read: " + name1);
                            logger.warningln("         No more reverse reads");

                        }
                        return false;
                    }
                    
                }
                else{
                    if(rec2 != nullptr){
                        if(!complainedAlready){
                            complainedAlready = true;
                            std::string name2 = rec2 -> name;
                            logger.warningln("WARNING: Pair validation failed at record: " + offset);
                            logger.warningln("         No more reverse reads");
                            logger.warningln("         Reverse read: " + name2);

                        }
                        return false;

                    }
                }
                
                offset++;
                return true;
                
            }
             
            
            // TODO validatePairs

    };
}

#endif