#ifndef TRIM_LOG_RECORD_H
#define TRIM_LOG_RECORD_H
#include <string>

namespace rabbit
{
    class TrimLogRecord
    {
    private:
        /* data */
        std::string readName;
        int length;
        int startPos;
        int endPos;
        int trimTail;
    public:
        TrimLogRecord(/* args */);
        ~TrimLogRecord();
    };
    
    
    
} // namespace rabbit


#endif
