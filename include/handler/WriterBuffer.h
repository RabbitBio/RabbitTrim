#ifndef RT_WRITER_BUFFER_H
#define RT_WRITER_BUFFER_H

#include "DataPool.h"
#include "DataQueue.h"

namespace rabbit
{
  namespace trim
  {
    struct WriterBuffer
    {
      char* data;
      unsigned int size;
      WriterBuffer():data(nullptr), size(0){}
      WriterBuffer(unsigned int size_):size(0){
        data = new char[size_];
      }
      ~WriterBuffer(){
        if(data != nullptr){
          delete[] data;
        }
      }
      void Reset(){
        size = 0;
      }
    };

    struct PEWriterBuffer
    {
      char* d1_p;
      char* d1_u;
      char* d2_p;
      char* d2_u;
      unsigned int size;
      PEWriterBuffer():d1_p(nullptr), d2_p(nullptr), d1_u(nullptr), d2_u(nullptr), size(0){}
      PEWriterBuffer(unsigned int size_):size(0){
        d1_p = new char[size_];
        d1_u = new char[size_];
        d2_p = new char[size_];
        d2_u = new char[size_];
      }
      ~PEWriterBuffer(){
        if(d1_p != nullptr){
          delete [] d1_p;
          delete [] d1_u;
          delete [] d2_p;
          delete [] d2_u;
        }
      }
    };

    // typedef rabbit::core::TDataPool<WriterBuffer> WriterDataPool;
    // typedef rabbit::core::TDataQueue<WriterBuffer> WriterDataQueue;
    // typedef rabbit::core::TDataQueue<PEWriterBuffer> PEWriterDataQueue;


  } // namespace trim

} // namespace rabbit

#endif
