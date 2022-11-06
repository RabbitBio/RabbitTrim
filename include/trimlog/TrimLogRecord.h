#ifndef TRIM_LOG_RECORD_H
#define TRIM_LOG_RECORD_H
#include <string>
#include "DataQueue.h"

namespace rabbit
{
    namespace log
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

      struct TrimLogBuffer{
        char* data;
        TrimLogBuffer() = default;
        TrimLogBuffer(int size_){
          data = new char[size_];
        }
        ~TrimLogBuffer(){
          if(data != nullptr)
            delete[] data;
        }
      };

      // typedef rabbit::core::TDataQueue<TrimLogBuffer> WriterDataQueue;

    } // namespace trim
    
} // namespace rabbit


#endif
