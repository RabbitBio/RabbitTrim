#ifndef RT_WRITER_BUFFER_H
#define RT_WRITER_BUFFER_H

#include "DataPool.h"
#include "DataQueue.h"

namespace rabbit
{
    namespace trim
    {
        typedef rabbit::core::TDataPool<WriterBuffer> WriterDataPool;
        typedef rabbit::core::TDataQueue<WriterBuffer> WriterDataQueue;
        struct WriterBuffer
        {
            char* data;
            unsigned int size;
            WriterBuffer():data(nullptr), size(0){
            }
            WriterBuffer(unsigned int size_){
                size = 0;
                data = new char[size_];
            }
            void Reset(){
                size = 0;
            }
        };
        
        
    } // namespace trim
    
} // namespace rabbit

#endif
