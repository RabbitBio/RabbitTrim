#include "ThreadAssign.h"
// using namespace rabbit::trim::util;

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
void rabbit::trim::util::threadAssignForSE(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp)
{
  int prefixPairCnt = 0; 
  int forwardCnt = 0;
  int reverseCnt = 0;
  int commonCnt = 0;
  bool palindromeKeepBoth;
  countAdapter(logger, rp, prefixPairCnt, forwardCnt, reverseCnt, commonCnt, palindromeKeepBoth);

  int tmp_worker_thread_num;
  int tmp_pigz_thread_num;
  int tmp_pragzip_thread_num;
  int max_thread_num;
  int min_thread_num;

  int cnt = forwardCnt + commonCnt;
  if(cnt == 0) cnt = 1;
  double tmp_rtrim_se_speed = rabbit::trim::RTRIM_SE_SPEED_PER_SEC / (cnt);

  if((rabbit::trim::util::endsWith(rp.output, ".gz") && rp.usePigz))
  {
    if(rabbit::trim::util::endsWith(rp.forwardFiles[0], ".gz"))
    {
      rp.usePragzip = true;
      // assign for worker, pigz and pragzip
      tmp_worker_thread_num = std::ceil(rp.threads / (1.0 + tmp_rtrim_se_speed / rabbit::trim::PIGZ_SPEED_PER_SEC + tmp_rtrim_se_speed / rabbit::trim::PRAGZIP_SPEED_PER_SEC));
      tmp_pragzip_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PRAGZIP_SPEED_PER_SEC / tmp_rtrim_se_speed + rabbit::trim::PRAGZIP_SPEED_PER_SEC / rabbit::trim::PIGZ_SPEED_PER_SEC));
      tmp_pigz_thread_num = rp.threads - tmp_worker_thread_num - tmp_pragzip_thread_num;
      double min_speed = std::min(tmp_worker_thread_num * tmp_rtrim_se_speed, tmp_pragzip_thread_num * rabbit::trim::PRAGZIP_SPEED_PER_SEC );
      min_speed = std::min(min_speed, tmp_pigz_thread_num * rabbit::trim::PIGZ_SPEED_PER_SEC );
      if(min_speed < rabbit::trim::IGZIP_DECOMPRESS_SPEED)
      {
        // assign for worker and pigzer
        tmp_pigz_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PIGZ_SPEED_PER_SEC / tmp_rtrim_se_speed ));
        tmp_worker_thread_num = rp.threads - tmp_pigz_thread_num;
        max_thread_num = std::ceil(rabbit::trim::IGZIP_DECOMPRESS_SPEED / tmp_rtrim_se_speed);
        if(tmp_worker_thread_num > max_thread_num) 
        {
          tmp_worker_thread_num = max_thread_num;
          tmp_pigz_thread_num = rp.threads - tmp_worker_thread_num;
        }
        if(tmp_pigz_thread_num * rabbit::trim::PIGZ_SPEED_PER_SEC >= rabbit::trim::IGZIP_COMPRESS_SPEED)
        {
          rp.usePragzip = false;
          tmp_pragzip_thread_num = 0;
        }
        else
        {
          rp.usePragzip = false;
          rp.usePigz = false;
          tmp_pragzip_thread_num = 0;
          tmp_pigz_thread_num = 0;
          tmp_worker_thread_num = rp.threads;
        }

      }

    }
    else
    {
      rp.usePragzip = false;
      tmp_pragzip_thread_num = 0;
      // assign for worker and pigz
      min_thread_num = rabbit::trim::MIN_PIGZ_THREAD_NUM + std::ceil(rabbit::trim::IGZIP_COMPRESS_SPEED / tmp_rtrim_se_speed);
      if(rp.threads < min_thread_num)
      {
        rp.usePigz = false;
        tmp_worker_thread_num = rp.threads;
        tmp_pigz_thread_num = 0;
        tmp_pragzip_thread_num = 0;
      }
      else
      {
        // double tmp_ = (tmp_rtrim_se_speed + rabbit::trim::PIGZ_SPEED_PER_SEC) / 2;
        // tmp_pigz_thread_num = std::floor(tmp_ - (std::sqrt(std::pow(tmp_, 2) - tmp_rtrim_se_speed * rp.threads)));
        tmp_pigz_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PIGZ_SPEED_PER_SEC / tmp_rtrim_se_speed ));
        tmp_worker_thread_num = rp.threads - tmp_pigz_thread_num;
        max_thread_num = std::ceil(rabbit::trim::INPUT_READ_SPEED / tmp_rtrim_se_speed);
        if(tmp_worker_thread_num > max_thread_num) 
        {
          tmp_worker_thread_num = max_thread_num;
          tmp_pigz_thread_num = rp.threads - tmp_worker_thread_num;
        }

        if(tmp_pigz_thread_num < rabbit::trim::MIN_PIGZ_THREAD_NUM)
        {
          rp.usePigz = false;
          tmp_pigz_thread_num = 0;
          tmp_worker_thread_num = rp.threads;
        }
      }

    }
  }
  else
  { // output is not gzip
    if(rabbit::trim::util::endsWith(rp.forwardFiles[0], ".gz"))
    {
      rp.usePragzip = true;
      rp.usePigz = false;
      tmp_pigz_thread_num = 0;
      // assign for worker and pragzip
      tmp_pragzip_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PRAGZIP_SPEED_PER_SEC / tmp_rtrim_se_speed));
      tmp_worker_thread_num = rp.threads - tmp_pragzip_thread_num;
      max_thread_num = std::ceil(rabbit::trim::INPUT_READ_SPEED / tmp_rtrim_se_speed);
      if(tmp_worker_thread_num > max_thread_num) 
      {
        tmp_worker_thread_num = max_thread_num;
        tmp_pragzip_thread_num = rp.threads - tmp_worker_thread_num;
      }

      if(tmp_pragzip_thread_num < rabbit::trim::MIN_PRAGZIP_THREAD_NUM)
      {
        rp.usePragzip = false;
        tmp_pragzip_thread_num = 0;
        tmp_worker_thread_num = rp.threads;
      }

    }
    else
    {
      rp.usePigz = false;
      rp.usePragzip = false;
      tmp_worker_thread_num = rp.threads;
      tmp_pigz_thread_num = 0;
      tmp_pragzip_thread_num = 0;
    }
  }

  rp.workerThreadNum = tmp_worker_thread_num;
  rp.pigzThreadsNum = tmp_pigz_thread_num;
  rp.pragzipThreadsNum = tmp_pragzip_thread_num;
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

void rabbit::trim::util::threadAssignForPE(rabbit::Logger& logger, rabbit::trim::RabbitTrimParam& rp)
{

  int prefixPairCnt = 0; 
  int forwardCnt = 0;
  int reverseCnt = 0;
  int commonCnt = 0;
  bool palindromeKeepBoth;
  countAdapter(logger, rp, prefixPairCnt, forwardCnt, reverseCnt, commonCnt, palindromeKeepBoth);

  int tmp_worker_thread_num;
  int tmp_pigz_thread_num;
  int tmp_pragzip_thread_num;
  int max_thread_num;
  int min_thread_num;

  // PE
  double tmp_rtrim_pe_speed;
  int cnt = forwardCnt + commonCnt;
  if(palindromeKeepBoth)
  {
    cnt += reverseCnt + commonCnt;
  }
  if(prefixPairCnt == 0)
  {
    if(cnt == 0)
    {
      tmp_rtrim_pe_speed = rabbit::trim::RTRIM_PE_SPEED_PER_SEC;
    }
    else
    {
      tmp_rtrim_pe_speed = rabbit::trim::RTRIM_SE_SPEED_PER_SEC / (cnt);

    }
  }
  else
  {
    if(cnt == 0)
    {
      tmp_rtrim_pe_speed = rabbit::trim::RTRIM_PE_SPEED_PER_SEC / prefixPairCnt;
    }
    else
    {
      double tmp_pair_speed = rabbit::trim::RTRIM_PE_SPEED_PER_SEC / prefixPairCnt;
      double tmp_single_speed = rabbit::trim::RTRIM_SE_SPEED_PER_SEC / (cnt);
      tmp_rtrim_pe_speed = 1 / ((1 / tmp_pair_speed) + (1 / tmp_single_speed));
    }
  }

  if((rabbit::trim::util::endsWith(rp.output, ".gz") && rp.usePigz))
  {
    if(rabbit::trim::util::endsWith(rp.forwardFiles[0], ".gz"))
    {
      rp.usePragzip = true;
      // assign for worker, pigz and pragzip
      tmp_worker_thread_num = std::ceil(rp.threads / (1.0 + tmp_rtrim_pe_speed / rabbit::trim::PIGZ_SPEED_PER_SEC + tmp_rtrim_pe_speed / rabbit::trim::PRAGZIP_SPEED_PER_SEC));
      tmp_pragzip_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PRAGZIP_SPEED_PER_SEC / tmp_rtrim_pe_speed + rabbit::trim::PRAGZIP_SPEED_PER_SEC / rabbit::trim::PIGZ_SPEED_PER_SEC));
      tmp_pigz_thread_num = rp.threads - tmp_worker_thread_num - tmp_pragzip_thread_num;
      double min_speed = std::min(tmp_worker_thread_num * tmp_rtrim_pe_speed, tmp_pragzip_thread_num * rabbit::trim::PRAGZIP_SPEED_PER_SEC );
      min_speed = std::min(min_speed, tmp_pigz_thread_num * rabbit::trim::PIGZ_SPEED_PER_SEC );
      if(min_speed < rabbit::trim::IGZIP_DECOMPRESS_SPEED)
      {
        // assign for worker and pigzer
        tmp_pigz_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PIGZ_SPEED_PER_SEC / tmp_rtrim_pe_speed ));
        tmp_worker_thread_num = rp.threads - tmp_pigz_thread_num;
        max_thread_num = std::ceil(rabbit::trim::IGZIP_DECOMPRESS_SPEED / tmp_rtrim_pe_speed);
        if(tmp_worker_thread_num > max_thread_num) 
        {
          tmp_worker_thread_num = max_thread_num;
          tmp_pigz_thread_num = rp.threads - tmp_worker_thread_num;
        }
        if(tmp_pigz_thread_num * rabbit::trim::PIGZ_SPEED_PER_SEC >= rabbit::trim::IGZIP_COMPRESS_SPEED)
        {
          rp.usePragzip = false;
          tmp_pragzip_thread_num = 0;
        }
        else
        {
          rp.usePragzip = false;
          rp.usePigz = false;
          tmp_pragzip_thread_num = 0;
          tmp_pigz_thread_num = 0;
          tmp_worker_thread_num = rp.threads;
        }

      }

    }
    else
    {
      // intput is not gzip
      rp.usePragzip = false;
      tmp_pragzip_thread_num = 0;
      // assign for worker and pigz
      min_thread_num = rabbit::trim::MIN_PIGZ_THREAD_NUM + std::ceil(rabbit::trim::IGZIP_COMPRESS_SPEED / tmp_rtrim_pe_speed);
      if(rp.threads < min_thread_num)
      {
        rp.usePigz = false;
        tmp_worker_thread_num = rp.threads;
        tmp_pigz_thread_num = 0;
        tmp_pragzip_thread_num = 0;
      }
      else
      {
        // double tmp_ = (tmp_rtrim_pe_speed + rabbit::trim::PIGZ_SPEED_PER_SEC) / 2;
        // tmp_pigz_thread_num = std::floor(tmp_ - (std::sqrt(std::pow(tmp_, 2) - tmp_rtrim_pe_speed * rp.threads)));
        tmp_pigz_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PIGZ_SPEED_PER_SEC / tmp_rtrim_pe_speed ));
        tmp_worker_thread_num = rp.threads - tmp_pigz_thread_num;
        max_thread_num = std::ceil(rabbit::trim::INPUT_READ_SPEED / tmp_rtrim_pe_speed);
        if(tmp_worker_thread_num > max_thread_num) 
        {
          tmp_worker_thread_num = max_thread_num;
          tmp_pigz_thread_num = rp.threads - tmp_worker_thread_num;
        }

        if(tmp_pigz_thread_num < rabbit::trim::MIN_PIGZ_THREAD_NUM)
        {
          rp.usePigz = false;
          tmp_pigz_thread_num = 0;
          tmp_worker_thread_num = rp.threads;
        }
      }
    }

  }
  else
  { // output is not gzip
    if(rabbit::trim::util::endsWith(rp.forwardFiles[0], ".gz"))
    {
      rp.usePragzip = true;
      rp.usePigz = false;
      tmp_pigz_thread_num = 0;
      // assign for worker and pragzip
      tmp_pragzip_thread_num = std::ceil(rp.threads / (1.0 + rabbit::trim::PRAGZIP_SPEED_PER_SEC / tmp_rtrim_pe_speed));
      tmp_worker_thread_num = rp.threads - tmp_pragzip_thread_num;
      max_thread_num = std::ceil(rabbit::trim::INPUT_READ_SPEED / tmp_rtrim_pe_speed);
      if(tmp_worker_thread_num > max_thread_num) 
      {
        tmp_worker_thread_num = max_thread_num;
        tmp_pragzip_thread_num = rp.threads - tmp_worker_thread_num;
      }

      if(tmp_pragzip_thread_num < rabbit::trim::MIN_PRAGZIP_THREAD_NUM)
      {
        rp.usePragzip = false;
        tmp_pragzip_thread_num = 0;
        tmp_worker_thread_num = rp.threads;
      }

    }
    else
    {
      rp.usePigz = false;
      rp.usePragzip = false;
      tmp_worker_thread_num = rp.threads;
      tmp_pigz_thread_num = 0;
      tmp_pragzip_thread_num = 0;
    }
  }

  rp.workerThreadNum = tmp_worker_thread_num;
  rp.pigzThreadsNum = tmp_pigz_thread_num;
  rp.pragzipThreadsNum = tmp_pragzip_thread_num;
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

