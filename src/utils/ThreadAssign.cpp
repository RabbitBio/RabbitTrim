#include "ThreadAssign.h"
// using namespace rabbit::trim::util;

template <typename T>
void rabbit::trim::util::genCoefVec(vector<T>& coefVec, const T& coef)
{
    coefVec.push_back(coef);
}

template <typename T, typename... Args>
void rabbit::trim::util::genCoefVec(vector<T>& coefVec, const T& coef, const Args&... args)
{
    coefVec.push_back(coef);
      genCoefVec(coefVec, args...);
}

double rabbit::trim::util::cal_trim_speed(int forwardCnt, int reverseCnt, int commonCnt, int prefixPairCnt, bool palindromeKeepBoth, bool isPE)
{
  double trim_speed;
  if(isPE) 
  {
    int cnt = forwardCnt + commonCnt;
    if(palindromeKeepBoth)
    {
      cnt += reverseCnt + commonCnt;
    }
    if(prefixPairCnt == 0)
    {
      if(cnt == 0)
      {
        trim_speed = rabbit::trim::RTRIM_PE_SPEED_PER_SEC;
      }
      else
      {
        trim_speed = rabbit::trim::RTRIM_SE_SPEED_PER_SEC / (cnt);

      }
    }
    else
    {
      if(cnt == 0)
      {
        trim_speed = rabbit::trim::RTRIM_PE_SPEED_PER_SEC / prefixPairCnt;
      }
      else
      {
        double tmp_pair_speed = rabbit::trim::RTRIM_PE_SPEED_PER_SEC / prefixPairCnt;
        double tmp_single_speed = rabbit::trim::RTRIM_SE_SPEED_PER_SEC / (cnt);
        trim_speed = 1 / ((1 / tmp_pair_speed) + (1 / tmp_single_speed));
      }
    }
  }
  else // SE
  {
    int cnt = forwardCnt + commonCnt;
    if(cnt == 0) cnt = 1;
    trim_speed = rabbit::trim::RTRIM_SE_SPEED_PER_SEC / (cnt);
  }
  return trim_speed;

}

double rabbit::trim::util::cal_speed(vector<double>& coef_1, vector<double>& coef_2, int x)
{
  if(x <= 0) return 0.0;
  double speed = 0.0;
  double x_d = 1.0;
  if(x <= 20)
  {
    for(auto coef : coef_1)
    {
      speed += coef * x_d;
      x_d *= x;
    }
  }
  else
  {
    for(auto coef : coef_2)
    {
      speed += coef * x_d;
      x_d *= x;
    }
  }
  return speed;
}

vector<int> rabbit::trim::util::cal_thread_num_1(vector<double>& coef_1, vector<double>& coef_2, double trim_speed_per_thread, int T)
{
  if(T <= 2) return {1, 1};
  int x = 1;
  int y = T - x;
  while(x > 0 && y > 1)
  {
    double compression_speed = cal_speed(coef_1, coef_2, x) * x;
    double trim_speed = trim_speed_per_thread * y;
    if(compression_speed < trim_speed)
    {
      x++;
      y--;
    }
    else break;
  }
  // judge which case is better 
  double tmp_compression_speed = cal_speed(coef_1, coef_2, x) * x;
  double tmp_trim_speed = trim_speed_per_thread * y;
  if(tmp_compression_speed >= tmp_trim_speed && x > 1)
  {
    double diff_1 = tmp_compression_speed - tmp_trim_speed;
    double diff_2 = (trim_speed_per_thread * (y + 1)) - (cal_speed(coef_1, coef_2, x - 1) * (x - 1));
    if(diff_2 <= diff_1)
    {
      x--;
      y++;
    }
  }
  return {x, y};
}

vector<int> rabbit::trim::util::cal_thread_num_2(vector<double>& pigz_coef_1, vector<double>& pigz_coef_2, vector<double>& pragzip_coef_1, vector<double>& pragzip_coef_2, double trim_speed_per_thread, int T)
{
  if(T <= 3) return {1, 1, 1};
  // x: pigz, y: pragzip, z: trim
  int x = 1, y = 1;
  int z = T - x - y;
  while(x > 0 & y > 0 & z > 1)
  {
    double pigz_speed = cal_speed(pigz_coef_1, pigz_coef_2, x) * x;
    double pragzip_speed = cal_speed(pragzip_coef_1, pragzip_coef_2, y) * y;
    double trim_speed = trim_speed_per_thread * z;
    if(pigz_speed < pragzip_speed)
    {
      if(trim_speed > pragzip_speed)
      {
        x++;
        z--;
      }
      else break;
    }
    else
    {
      if(trim_speed > pigz_speed)
      {
        y++;
        z--;
      }
      else break;
    }
  }
  // adjust
  double pigz_speed =  cal_speed(pigz_coef_1, pigz_coef_2, x) * x;
  double pragzip_speed = cal_speed(pragzip_coef_1, pragzip_coef_2, y) * y;
  double trim_speed = trim_speed_per_thread * z;
  if(pigz_speed >= pragzip_speed)
  {
    if(pigz_speed >= trim_speed && x > 1)
    {
      double diff_1 = pigz_speed - trim_speed;
      double diff_2 = (trim_speed_per_thread * (z + 1)) - (cal_speed(pigz_coef_1, pigz_coef_2, x - 1) * (x - 1));
      if(diff_2 <= diff_1)
      {
        x--;
        z++;
      }
    }
  }
  else
  {
    if(pragzip_speed >= trim_speed && y > 1)
    {
      double diff_1 = pragzip_speed - trim_speed;
      double diff_2 = (trim_speed_per_thread * (z + 1)) - (cal_speed(pragzip_coef_1, pragzip_coef_2, y - 1) * (y - 1));
      if(diff_2 <= diff_1)
      {
        y--;
        z++;
      }
    }
  }
  return {x, y, z};

}

void rabbit::trim::util::countAdapter(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp, int& prefixPairCnt, int& forwardCnt, int& reverseCnt, int& commonCnt, bool& palindromeKeepBoth_)
{
  const std::string PREFIX = "Prefix"; 
  const std::string SUFFIX_F = "/1";
  const std::string SUFFIX_R = "/2";
  for(auto step : rp.steps)
  {
    std::string trimmerName;
    std::string trimmerArgs;
    std::size_t idx = step.find(":");
    trimmerName = step.substr(0,idx);
    if(idx != std::string::npos)
      trimmerArgs = step.substr(idx+1);
    if(trimmerName.compare("ILLUMINACLIP") == 0) 
    {
      std::string fastaAdapterFile;
      int seedMaxMiss;
      int minPalindromeLikelihood;
      int minSequenceLikelihood;
      int minPrefix = 1;
      bool palindromeKeepBoth = false; 

      std::size_t pos = trimmerArgs.find(":");
      fastaAdapterFile = trimmerArgs.substr(0, pos);
      std::string::size_type sz;
      std::string::size_type cur_pos = pos + 1;

      seedMaxMiss = std::stoi(trimmerArgs.substr(cur_pos), &sz);
      cur_pos += sz + 1;
      minPalindromeLikelihood = std::stoi(trimmerArgs.substr(cur_pos), &sz);
      cur_pos += sz + 1;
      minSequenceLikelihood = std::stoi(trimmerArgs.substr(cur_pos), &sz);
      cur_pos += sz + 1;

      if(trimmerArgs.size() > cur_pos){
        minPrefix = std::stoi(trimmerArgs.substr(cur_pos), &sz);
        cur_pos += sz + 1;
        if(trimmerArgs.size() > cur_pos){
          std::string keepBothStr = trimmerArgs.substr(cur_pos);
          if(keepBothStr.compare("true") == 0) palindromeKeepBoth = true;
        }
      }
      palindromeKeepBoth_ = palindromeKeepBoth;

      std::ifstream fastaAdapter(fastaAdapterFile.c_str(), std::ifstream::in);
      if(!fastaAdapter.is_open()){
        logger.errorln("Can not open fastaAdapterFile : " + fastaAdapterFile);
        exit(0);
      }

      // read fasta file
      std::map<std::string, std::string> forwardSeqMap;
      std::map<std::string, std::string> reverseSeqMap;
      std::map<std::string, std::string> commonSeqMap;
      std::set<std::string> forwardPrefix;
      std::set<std::string> reversePrefix;
      std::set<std::string> prefixSet;
      std::string curLine;
      std::string name;
      std::string sequence;
      std::getline(fastaAdapter, curLine);
      // get all adapter
      while(!fastaAdapter.eof()){
        rabbit::trim::util::ClearHeadTailSpace(curLine);
        while(curLine.size() && curLine.at(0) != '>'){
          if(!std::getline(fastaAdapter, curLine)) break;
          rabbit::trim::util::ClearHeadTailSpace(curLine);
        }
        if(curLine.size() && curLine.at(0) == '>'){
          std::string fullName = curLine.substr(1);
          std::vector<std::string> tokens = rabbit::trim::util::split(fullName, "[\\| ]");
          name = tokens[0];

          sequence = "";
          while(std::getline(fastaAdapter, curLine)){
            rabbit::trim::util::ClearHeadTailSpace(curLine);
            if(curLine.size() == 0) continue;
            if(curLine.at(0) == '>') break;
            if(curLine.at(0) == ';') continue;
            sequence += curLine;
          }
        }
        if(sequence.size() == 0) continue;

        if(rabbit::trim::util::endsWith(name, SUFFIX_F)){
          forwardSeqMap.insert(std::make_pair(name, sequence));
          if(rabbit::trim::util::startsWith(name, PREFIX)){
            forwardPrefix.insert(name.substr(0, name.size() - SUFFIX_F.size()));
          }
        }
        else{
          if(rabbit::trim::util::endsWith(name, SUFFIX_R)){
            reverseSeqMap.insert(std::make_pair(name, sequence));
            if(rabbit::trim::util::startsWith(name, PREFIX)){
              reversePrefix.insert(name.substr(0, name.size() - SUFFIX_R.size()));
            }
          }
          else{
            commonSeqMap.insert(std::make_pair(name, sequence));
          }
        }

      }


      // load sequence
      // find same elements in forwardPrefix and reversePrefix
      for(auto iter = forwardPrefix.begin(); iter != forwardPrefix.end(); iter++){
        if(reversePrefix.count(*iter))
        {
          prefixSet.insert(*iter);
          prefixPairCnt++;
        }
      }
      for(auto iter = prefixSet.begin(); iter != prefixSet.end(); iter++){
        std::string forwardName = *iter + SUFFIX_F;
        std::string reverseName = *iter + SUFFIX_R;

        std::string forwardRec = forwardSeqMap[forwardName];
        std::string reverseRec = reverseSeqMap[reverseName];
        // delete same elements in forwardSeqMap and reverseSeqMap
        forwardSeqMap.erase(forwardName);
        reverseSeqMap.erase(reverseName);
      }

      // mapClippingSet
      // 01 forwardSeqMap 
      for(auto iter = forwardSeqMap.begin(); iter != forwardSeqMap.end(); iter++){
        std::set<std::string> uniqueSeq;
        std::string sequence =  iter -> second;
        if(uniqueSeq.count(sequence)) continue;
        else{
          uniqueSeq.insert(sequence);
          forwardCnt++;
        }

      }
      // 02 reverseSeqMap 
      for(auto iter = reverseSeqMap.begin(); iter != reverseSeqMap.end(); iter++){
        std::set<std::string> uniqueSeq;
        std::string sequence =  iter -> second;
        if(uniqueSeq.count(sequence)) continue;
        else{
          uniqueSeq.insert(sequence);
          reverseCnt++;
        }
      }
      // 03 commonSeqMap 
      for(auto iter = commonSeqMap.begin(); iter != commonSeqMap.end(); iter++){
        std::set<std::string> uniqueSeq;
        std::string sequence =  iter -> second;
        if(uniqueSeq.count(sequence)) continue;
        else{
          uniqueSeq.insert(sequence);
          commonCnt++;
        }
      }

    }
  }
}
void rabbit::trim::util::threadAssign(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp, bool isPE)
{
  int prefixPairCnt = 0; 
  int forwardCnt = 0;
  int reverseCnt = 0;
  int commonCnt = 0;
  bool palindromeKeepBoth;
  countAdapter(logger, rp, prefixPairCnt, forwardCnt, reverseCnt, commonCnt, palindromeKeepBoth);
    
  double trim_speed_per_thread = cal_trim_speed(forwardCnt, reverseCnt, commonCnt, prefixPairCnt, palindromeKeepBoth, isPE);

  // get coefficients vector
  vector<double> pragzip_coef_vec_1;
  vector<double> pragzip_coef_vec_2;
  vector<double> pigz_coef_vec_1;
  vector<double> pigz_coef_vec_2;
  pragzip_coef_vec_1.reserve(4);
  pragzip_coef_vec_2.reserve(3);
  pigz_coef_vec_1.reserve(3);
  pigz_coef_vec_2.reserve(4);
  genCoefVec(pragzip_coef_vec_1, PRAGZIP_COEF_0, PRAGZIP_COEF_1, PRAGZIP_COEF_2, PRAGZIP_COEF_3);
  genCoefVec(pragzip_coef_vec_2, PRAGZIP_COEF_4, PRAGZIP_COEF_5, PRAGZIP_COEF_6);
  genCoefVec(pigz_coef_vec_1, PIGZ_COEF_0, PIGZ_COEF_1, PIGZ_COEF_2);
  genCoefVec(pigz_coef_vec_2, PIGZ_COEF_3, PIGZ_COEF_4, PIGZ_COEF_5, PIGZ_COEF_6);


  if((rabbit::trim::util::endsWith(rp.output, ".gz") && rp.usePigz))
  {
    if(rabbit::trim::util::endsWith(rp.forwardFiles[0], ".gz"))
    {
      double max_speed;
      // using pragzip and pigz
      vector<int> res_1 = cal_thread_num_2(pigz_coef_vec_1, pigz_coef_vec_2, pragzip_coef_vec_1, pragzip_coef_vec_2, trim_speed_per_thread, rp.threads);
      double pigz_speed_1 = cal_speed(pigz_coef_vec_1, pigz_coef_vec_2, res_1[0]) * res_1[0];
      double pragzip_speed_1 = cal_speed(pragzip_coef_vec_1, pragzip_coef_vec_2, res_1[1]) * res_1[1];
      double trim_speed_1 = trim_speed_per_thread * res_1[2];
      double min_speed_1 = std::min(std::min(pigz_speed_1, pragzip_speed_1), trim_speed_1);
      max_speed = min_speed_1;
      rp.usePragzip = true;
      rp.usePigz = true;
      rp.pigzThreadsNum = res_1[0];
      rp.pragzipThreadsNum = res_1[1];
      rp.workerThreadNum = res_1[2];
      // using pragzip and igzip
      vector<int> res_2 = cal_thread_num_1(pragzip_coef_vec_1, pragzip_coef_vec_2, trim_speed_per_thread, rp.threads - 1);
      double pigz_speed_2 = rabbit::trim::IGZIP_COMPRESS_SPEED;
      double pragzip_speed_2 = cal_speed(pragzip_coef_vec_1, pragzip_coef_vec_2, res_2[0]) * res_2[0];
      double trim_speed_2 = trim_speed_per_thread * res_2[1];
      double min_speed_2 = std::min(std::min(pigz_speed_2, pragzip_speed_2), trim_speed_2);
      if(min_speed_2 > max_speed)
      {
        max_speed = min_speed_2;
        rp.usePragzip = true;
        rp.usePigz = false;
        rp.pigzThreadsNum = 0;
        rp.pragzipThreadsNum = res_2[0];
        rp.workerThreadNum = res_2[1];
      }
      // using igzip and pigz
      vector<int> res_3 = cal_thread_num_1(pigz_coef_vec_1, pigz_coef_vec_2, trim_speed_per_thread, rp.threads - 1);
      double pigz_speed_3 = cal_speed(pigz_coef_vec_1, pigz_coef_vec_2, res_3[0]) * res_3[0];
      double pragzip_speed_3 = rabbit::trim::IGZIP_DECOMPRESS_SPEED;
      double trim_speed_3 = trim_speed_per_thread * res_3[1];
      double min_speed_3 = std::min(std::min(pigz_speed_3, pragzip_speed_3), trim_speed_3);
      if(min_speed_3 > max_speed)
      {
        max_speed = min_speed_3;
        rp.usePragzip = false;
        rp.usePigz = true;
        rp.pigzThreadsNum = res_3[0];
        rp.pragzipThreadsNum = 0;
        rp.workerThreadNum = res_3[1];
      }
      // using igzip and igzip
      double pigz_speed_4 = rabbit::trim::IGZIP_COMPRESS_SPEED;
      double pragzip_speed_4 = rabbit::trim::IGZIP_DECOMPRESS_SPEED;
      double trim_speed_4 = trim_speed_per_thread * (rp.threads > 2 ?  rp.threads - 2 : 1); 
      double min_speed_4 = std::min(std::min(pigz_speed_4, pragzip_speed_4), trim_speed_4);
      if(min_speed_4 >= max_speed)
      {
        max_speed = min_speed_4;
        rp.usePragzip = false;
        rp.usePigz = false;
        rp.pigzThreadsNum = 0;
        rp.pragzipThreadsNum = 0;
        rp.workerThreadNum = (rp.threads > 2 ? rp.threads - 2 : 1);
      }
      
    }
    else
    {
      // assign for trimmer and pigz
      rp.usePragzip = false;
      rp.pragzipThreadsNum = 0;
      vector<int> res = cal_thread_num_1(pigz_coef_vec_1, pigz_coef_vec_2, trim_speed_per_thread, rp.threads);
      double pigz_speed = cal_speed(pigz_coef_vec_1, pigz_coef_vec_2, res[0]) * res[0];
      double trim_speed = trim_speed_per_thread * res[1];
      rp.usePigz = true;
      rp.pigzThreadsNum  = res[0];
      rp.workerThreadNum = res[1];
      if(pigz_speed <= rabbit::trim::IGZIP_COMPRESS_SPEED)
      {
        rp.usePigz = false;
        rp.pigzThreadsNum = 0;
        rp.workerThreadNum = (rp.threads > 1 ? rp.threads - 1 : 1); 
      }
    }
  }
  else
  { // output is not gzip
    if(rabbit::trim::util::endsWith(rp.forwardFiles[0], ".gz"))
    {
      // assign for trimmer and pragzip
      rp.usePigz = false;
      rp.pigzThreadsNum = 0;
      vector<int> res = cal_thread_num_1(pragzip_coef_vec_1, pragzip_coef_vec_2, trim_speed_per_thread, rp.threads);
      double pragzip_speed = cal_speed(pragzip_coef_vec_1, pragzip_coef_vec_2 , res[0]) * res[0];
      double trim_speed = trim_speed_per_thread * res[1];
      rp.usePragzip = true;
      rp.pragzipThreadsNum = res[0];
      rp.workerThreadNum   = res[1];
      if(pragzip_speed <= rabbit::trim::IGZIP_DECOMPRESS_SPEED)
      {
        rp.usePragzip = false;
        rp.pragzipThreadsNum = 0;
        rp.workerThreadNum = (rp.threads > 1 ? rp.threads - 1 : 1); 
      }
    }
    else
    {
      rp.usePigz = false;
      rp.usePragzip = false;
      rp.workerThreadNum = rp.threads;
      rp.pigzThreadsNum = 0;
      rp.pragzipThreadsNum = 0;
    }
  }
  if(rp.usePragzip && rp.pragzipThreadsNum == 1) rp.usePragzip = false;

#ifdef TRIM_DEBUG
  logger.debugln("======== Thread Assign Result ======== ");
  logger.debugln("total thread num : " + std::to_string(rp.threads));
  logger.debugln("use pragzip: " + std::to_string(rp.usePragzip));
  logger.debugln("pragzip thread num : " + std::to_string(rp.pragzipThreadsNum));
  logger.debugln("worker thread num : " + std::to_string(rp.workerThreadNum));
  logger.debugln("use pigz: " + std::to_string(rp.usePigz));
  logger.debugln("pigz thread num : " + std::to_string(rp.pigzThreadsNum));
  logger.debugln("======== Thread Assign Result ======== ");
#endif
}

