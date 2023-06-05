
#ifndef _RABBIT_TRIM_PRAGZIP_H
#define _RABBIT_TRIM_PRAGZIP_H


#include "DataQueue.h"

typedef rabbit::core::TDataQueue<std::pair<char*, int> > PragzipQueue;

int main_pragzip(int argc, char *argv[], PragzipQueue& pragzipQueue);

#endif
