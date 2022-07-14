#ifndef FASTQ_NAME_PATTERN_H
#define FASTQ_NAME_PATTERN_H

#include <string>
#include <regex> 
#include <vector>

namespace rabbit{
    const int fastqNamePatternNums = 2;
    constexpr std::pair<std::string, std::string> casava_13{"(.* )?([^:]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+#[A-Z0-9]+).*", "$2"}; 
    constexpr std::pair<std::string, std::string> casava_18 = std::make_pair("(.* )?([^:]+:[0-9]+:[A-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+).*","$2");
    std::vector<std::pair<std::string, std::string>> fastqNamePatternArray{
        casava_13, casava_18
    };
    
    class FastqNamePattern{
        public:
            FastqNamePattern(std::string patternStr, std::string replacement_):base_regex(patternStr){
                replacement = replacement_;
            }
            
            bool match(std::string str){
                return std::regex_match(str, base_match, base_regex);
            }
            
            std::string canonicalizeOne(std::string str){
                if(!match(str)) 
                    return "";
                return std::regex_replace(str, base_regex, replacement);
            }
            
        
        private:
            std::regex base_regex;
            std::smatch base_match;
            std::string replacement;
            
    };
    
    // 匹配所有的fastqNamePattern
    std::string canonicalize(std::string str){
        for(int i = 0; i < fastqNamePatternNums; i++){
            FastqNamePattern fqNamePattern(FastqNamePatternArray[i].first, FastqNamePatternArray[i].second); // 构造FastqNamePattern对象 
            // TODO : 对于每一组PE Read 都需要验证，就会反复的构造FastqNamePattern对象，有点耗时
            std::string canon = fqNamePattern.canonicalizeOne(str);
            if(canon != "") return canon;
        }
        return ""; // casava_13 and casava_18 都不能匹配

    }
}

#endif