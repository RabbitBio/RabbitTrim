#ifndef _THREAD_ASSIGN_H
#define _THREAD_ASSIGN_H
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <cmath>
#include <vector>
#include "util.h"
#include "param.h"

namespace rabbit
{
  namespace trim
  {
    namespace util
    {
      template <typename T>
      void genCoefVec(vector<T>& coefVec, const T& coef);
      template <typename T, typename... Args>
      void genCoefVec(vector<T>& coefVec, const T& coef, const Args&... args);
      
      double cal_trim_speed(int forwardCnt, int reverseCnt, int commonCnt, int prefixPairCnt, bool palindromeKeepBoth, bool isPE);
      double cal_speed(vector<double>& coef_1, vector<double>& coef_2, int x);
      vector<int> cal_thread_num_1(vector<double>& coef_1, vector<double>& coef_2, double trim_speed_per_thread, int T);
      vector<int> cal_thread_num_2(vector<double>& pigz_coef_1, vector<double>& pigz_coef_2, vector<double>& pragzip_coef_1, vector<double>& pragzip_coef_2, double trim_speed_per_thread, int T);

      void countAdapter(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp, int& prefixPairCnt, int& forwardCnt, int& reverseCnt, int& commonCnt, bool& palindromeKeepBoth_);
      void threadAssign(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp, bool isPE);
    } //namespace util
  } // namespace trim
} // namespace rabbit

#endif

