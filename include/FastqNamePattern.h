#ifndef FASTQ_NAME_PATTERN_H
#define FASTQ_NAME_PATTERN_H

#include <string>
#include <regex> 
#include <vector>

namespace rabbit{
    
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
            // $2 表示第二个分组匹配结果
            // std::vector<std::pair<std::string, std::string>> fastqNamePatternArray;
            
    };
    
//    // 匹配所有的fastqNamePattern
//    std::string canonicalize(std::string str){
//      
//            FastqNamePattern fqNamePattern_13(casava_13, "$2"); // 构造FastqNamePattern对象 
//            FastqNamePattern fqNamePattern_18(casava_18, "$2"); 
//            // TODO : 对于每一组PE Read 都需要验证，就会反复的构造FastqNamePattern对象，有点耗时
//            std::string canon = fqNamePattern_13.canonicalizeOne(str);
//            if(canon.size()) return canon;
//            canon = fqNamePattern_18.canonicalizeOne(str);
//            if(canon.size()) return canon;
//        }
//        return ""; // casava_13 and casava_18 都不能匹配
//
//    }
} // namespace rabbit

#endif
