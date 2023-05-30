#ifndef _THREAD_ASSIGN_H
#define _THREAD_ASSIGN_H
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <cmath>
#include "util.h"
#include "param.h"

namespace rabbit
{
  namespace trim
  {
    namespace util
    {
      void countAdapter(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp, int& prefixPairCnt, int& forwardCnt, int& reverseCnt, int& commonCnt, bool& palindromeKeepBoth_);
      void threadAssignForSE(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp);
      void threadAssignForPE(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp);
    } //namespace util
  } // namespace trim
} // namespace rabbit

#endif

