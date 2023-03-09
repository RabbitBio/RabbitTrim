#include "handler/HandlerSE.h"
#include "pigz.h"

using namespace rabbit::trim;

int rabbit::trim::process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
  // the number of consumer 
  int consumer_num = rp.threads;

  rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
  rabbit::trim::FastqDataChunkQueue queue1(256,1);
  rabbit::trim::WriterDataQueue queue2(512, consumer_num);
  TrimLogDataQueue queue3(256, consumer_num);
  // output DataPool
  WriterBufferDataPool* wbDataPool = new WriterBufferDataPool(512, MEM_PER_CHUNK);

  // Trim Stats for each thread
  std::vector<rabbit::log::TrimStat> statsArr(consumer_num);

  // // PairingValidator
  // rabbit::PairingValidator* pairingValidator;
  // if(rp.validatePairing)
  //   pairingValidator = new rabbit::PairingValidator(logger);

  // trimmers
  TrimmerFactory *trimmerFactory = new TrimmerFactory(logger);
  std::vector<Trimmer*> trimmers;
  trimmerFactory -> makeTrimmers(rp, rp.steps, trimmers);


  // producer
  std::thread producer(rabbit::trim::producer_se_task, std::ref(rp), std::ref(logger), fastqPool, std::ref(queue1));


  // consumer
  std::thread **consumer_threads = new std::thread* [consumer_num]; 
  for(int tn = 0; tn < consumer_num; tn++){
    // asm ("":::"memory");
    consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_se_task, std::ref(rp), fastqPool, std::ref(queue1), wbDataPool, std::ref(queue2), std::ref(queue3), std::ref(statsArr[tn]), std::ref(trimmers), tn));
  }

  // writer
  std::thread* writer = NULL;

  // pigzer
  std::thread* pigzer = NULL;
  std::pair<char*, int> pigzLast;
  if(rabbit::trim::util::endsWith(rp.output, ".gz") && rp.usePigz){
    std::string tmp_output = rp.output.substr(0, rp.output.find(".gz"));
    // logger.errorln(tmp_output);
    std::ofstream fout(tmp_output.c_str());
    if(fout.fail()){
      logger.errorln("Can not open file " + rp.output);
      exit(1);
    }
    fout.close();
    pigzLast.first = new char [ MEM_PER_CHUNK ];
    pigzLast.second = 0;
    pigzer = new std::thread(std::bind(rabbit::trim::pigzer_se_task, std::ref(rp),  wbDataPool, std::ref(queue2), std::ref(pigzLast)));
  }
  else{
    // writer 
    writer = new std::thread(std::bind(rabbit::trim::writer_se_task, std::ref(rp), wbDataPool, std::ref(queue2), std::ref(logger)));
  }

  // trim log
  std::thread* trimlog_thread;
  if(rp.trimLog.size()){
    trimlog_thread = new std::thread(std::bind(trimlog_se_task, std::ref(rp), std::ref(queue3), std::ref(logger)));
  }


  producer.join();
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn]->join();
  }
  if(writer) writer -> join();
  else pigzer -> join();
  if(rp.trimLog.size())
    trimlog_thread -> join();

  // write trim stats
  rabbit::log::TrimStat * totalStat = new rabbit::log::TrimStat(logger);
  totalStat -> merge(statsArr);
  if(rp.seqA.size()) totalStat -> print(rp.stats);
  else totalStat -> printSE(rp.stats);
  // auto aa = dynamic_cast<rabbit::trim::IlluminaClippingTrimmer*>(trimmers[0]);
  // aa -> printCnt();
  return 0;
}

// producer
int rabbit::trim::producer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataChunkQueue& dq){
  bool isReverse = (bool)rp.reverseFiles.size();
  int totalFileCount = isReverse ? rp.reverseFiles.size() : rp.forwardFiles.size();
  std::vector<std::string> inputFiles = isReverse ? rp.reverseFiles : rp.forwardFiles;
  rabbit::int64 n_chunks = 0;
  for(int fileCnt = 0; fileCnt < totalFileCount; fileCnt++){
    rabbit::fq::FastqFileReader * fqFileReader;
    bool file_is_gz = false;
    // unsigned int i = inputFiles[fileCnt].size() - 3;
    // const char * p = inputFiles[fileCnt].c_str();
    // if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
    if(rabbit::trim::util::endsWith(inputFiles[fileCnt], ".gz")){
      file_is_gz = true;
    }
    if(file_is_gz){
      fqFileReader = new rabbit::fq::FastqFileReader(inputFiles[fileCnt],*fastqPool, true, "");
    }else{
      fqFileReader = new rabbit::fq::FastqFileReader(inputFiles[fileCnt],*fastqPool, false, "");
    }

    while (true){
      rabbit::fq::FastqChunk *fqChunk = new rabbit::fq::FastqChunk; 
      fqChunk->chunk = fqFileReader->readNextChunk();
      if(fqChunk->chunk == NULL) break;
      dq.Push(n_chunks, fqChunk->chunk);
      n_chunks++;
    }
    delete fqFileReader;
  }
  dq.SetCompleted();
  return 0;
  }


  // comusmer task
  void rabbit::trim::consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterDataQueue& dq2, TrimLogDataQueue& dq3, rabbit::log::TrimStat& rstats, const std::vector<rabbit::trim::Trimmer*>& trimmers, int threadId)
  {
    bool isReverse = rp.reverseFiles.size();
    bool isLog = rp.trimLog.size();
    rabbit::int64 chunk_id;
    rabbit::fq::FastqChunk* fqChunk = new rabbit::fq::FastqChunk;
    // IlluminaClippingTrimmer* trimmer = dynamic_cast<IlluminaClippingTrimmer*>(trimmers[0]);
    while(dq.Pop(chunk_id,fqChunk->chunk))
    {
      // std::vector<Reference> data;
      std::vector<neoReference> data;
      data.reserve(1e4);
      int loaded  = rabbit::fq::chunkFormat(fqChunk->chunk, data, true);
      rstats.readsInput += loaded;

      for(auto trimmer : trimmers){
        trimmer -> processRecords(data, threadId, false, isReverse);
      }

      // copy data to WriterBuffer 
      // WriterBuffer* wb = new WriterBuffer((unsigned int)MEM_PER_CHUNK);
      WriterBuffer* wb = NULL;
      wbDataPool -> Acquire(wb);
      rabbit::log::TrimLogBuffer* tb;
      if(isLog) tb = new rabbit::log::TrimLogBuffer(MEM_PER_TRIMLOG_BUFFER);
      uint64_t pos = 0;
      for(auto rec : data)
      {
        if(rec.lorigin == 1) rstats.realHit++;
        if(rec.lorigin == 2) rstats.tailHit++;
        if(rec.lorigin == 3)
        {
          rstats.realHit++;
          rstats.dimer++;
        }
        if(rec.lseq == 0)
        {
          rstats.readsDropped++;
          continue;
        }
        rstats.readsSurvivingForward++;

        std::memcpy(wb->data + pos, (rec.base + rec.pname), rec.lname);
        pos += rec.lname;
        wb->data[pos++] = '\n';
        std::memcpy(wb->data + pos, (rec.base + rec.pseq), rec.lseq);
        pos += rec.lseq;
        wb->data[pos++] = '\n';
        std::memcpy(wb->data + pos, (rec.base + rec.pstrand), rec.lstrand);
        pos += rec.lstrand;
        wb->data[pos++] = '\n';
        std::memcpy(wb->data + pos, (rec.base + rec.pqual), rec.lqual);
        pos += rec.lqual;
        wb->data[pos++] = '\n';
      }
      wb->data[pos] = '\0';
      wb->size = pos;
      dq2.Push(chunk_id, wb);
      if(isLog)
      {
        uint64_t pos_ = 0;
        for(auto rec : data){
          std::memcpy(tb->data + pos_, (rec.base + (rec.pname + 1)), rec.lname - 1);
          pos_ += rec.lname - 1;
          tb->data[pos_++] = ' ';
          int length = 0;
          int startPos = 0;
          int endPos = 0;
          int trimTail = 0;
          if(rec.lseq != 0){
            length = rec.lseq;
            startPos = rec.pseq - rec.porigin;
            endPos = startPos + length;
            trimTail = rec.lorigin - endPos;
          }
          std::string tmp = std::to_string(length) + " " + std::to_string(startPos) + " " + std::to_string(endPos) + " " + std::to_string(trimTail) + "\n";
          std::memcpy(tb->data + pos_, tmp.c_str(), tmp.size());
          pos_ += tmp.size();
        }
        tb->data[pos_++] = '\0';
        dq3.Push(chunk_id, tb);
      }
      fastqPool->Release(fqChunk->chunk);
    }
    dq2.SetCompleted();
    dq3.SetCompleted();
  }

  void rabbit::trim::writer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterDataQueue& dq2, rabbit::Logger& logger){
    std::ofstream fout(rp.output.c_str());
    if(fout.fail()){
      logger.errorln("Can not open file " + rp.output);
      exit(1);
    }
    rabbit::int64 chunk_id;
    WriterBuffer* wb;
    while(dq2.Pop(chunk_id, wb)){
      // write to file
      fout << wb->data;
      wbDataPool -> Release(wb);
    }
    fout.close();
  } 

  void rabbit::trim::trimlog_se_task(rabbit::trim::RabbitTrimParam& rp, TrimLogDataQueue& dq, rabbit::Logger& logger){
    std::ofstream fout(rp.trimLog.c_str());
    if(fout.fail()){
      logger.errorln("Can not open file " + rp.trimLog);
      exit(1);
    }
    rabbit::int64 chunk_id;
    rabbit::log::TrimLogBuffer* tb = new rabbit::log::TrimLogBuffer;
    while(dq.Pop(chunk_id, tb)){
      // write to file
      fout << tb->data;
      // delelte ?
    }
    fout.close();

  }

  void rabbit::trim::pigzer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterDataQueue& dq2, std::pair<char*, int>& pigzLast)
  {
    /*
       argc 9
       argv ./pigz
       argv -p
       argv 16
       argv -k
       argv -4
       argv -f
       argv -b
       argv -4096
       argv p.fq
       */
    int cnt = 9;

    char **infos = new char *[9];
    infos[0] = "./pigz";
    infos[1] = "-p";
    int th_num = rp.pigzThreadsNum;
    std::string th_num_s = std::to_string(th_num);

    infos[2] = new char[th_num_s.length() + 1];
    std::memcpy(infos[2], th_num_s.c_str(), th_num_s.length());
    infos[2][th_num_s.length()] = '\0';
    infos[3] = "-k";


    std::string tmp_level = std::to_string(rp.compressLevel);
    tmp_level = "-" + tmp_level;
    infos[4] = new char[tmp_level.length() + 1];
    std::memcpy(infos[4], tmp_level.c_str(), tmp_level.length());
    infos[4][tmp_level.length()] = '\0';

    infos[5] = "-f";
    infos[6] = "-b";
    infos[7] = "4096";
    std::string out_name1 = rp.output;
    std::string out_file = out_name1.substr(0, out_name1.find(".gz"));

    infos[8] = new char[out_file.length() + 1];
    std::memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';
    main_pigz(cnt, infos, wbDataPool, dq2, pigzLast);
  }
