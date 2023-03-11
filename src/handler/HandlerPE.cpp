#include "handler/HandlerPE.h"
#include "pigz.h"

using namespace rabbit::trim;

int rabbit::trim::process_pe(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
  int consumer_num = rp.threads;
  rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(128, MEM_PER_CHUNK);
  rabbit::trim::FastqDataPairChunkQueue queue1(128,1);

  // output data pool
  rabbit::trim::WriterBufferDataPool* wbDataPool = new rabbit::trim::WriterBufferDataPool(128 << 2, MEM_PER_CHUNK); 
  rabbit::trim::WriterBufferDataQueue wbQueue1(128, consumer_num);
  rabbit::trim::WriterBufferDataQueue wbQueue2(128, consumer_num);
  rabbit::trim::WriterBufferDataQueue wbQueue3(128, consumer_num);
  rabbit::trim::WriterBufferDataQueue wbQueue4(128, consumer_num);
  rabbit::trim::WriterBufferDataQueue logQueue(128, consumer_num);

  std::vector<rabbit::log::TrimStat> statsArr(consumer_num);
  // PairingValidator
  rabbit::PairingValidator* pairingValidator;
  if(rp.validatePairing)
    pairingValidator = new rabbit::PairingValidator(logger);

  // trimmers
  TrimmerFactory *trimmerFactory = new TrimmerFactory(logger);
  std::vector<Trimmer*> trimmers;
  trimmerFactory -> makeTrimmers(rp, rp.steps, trimmers);

  // trim log
  std::thread* trimlog_thread = NULL;
  std::thread* trimlog_pigz = NULL;
  if(rp.trimLog.size()){
    // use pigz
    if(rabbit::trim::util::endsWith(rp.trimLog, ".gz") && rp.usePigz)
    {
      // the file must exist before write while using pigz 
      std::string out_file = rp.trimLog.substr(0, rp.trimLog.size() - 3);
      std::ofstream fout(out_file.c_str());
      if(fout.fail()){
        logger.errorln("Can not open file " + out_file);
        exit(1);
      }
      fout.close();
      pair<char*, int> pigzLast;
      pigzLast.first = new char[MEM_PER_CHUNK];
      pigzLast.second = 0;
      trimlog_pigz = new std::thread(std::bind(pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(logQueue), std::ref(pigzLast), out_file));
    }
    else{
      std::string out_file = rp.trimLog;
      trimlog_thread = new std::thread(std::bind(writer_pe_task, std::ref(rp), wbDataPool, std::ref(logQueue), out_file, std::ref(logger)));
    }
  }

  // producer
  std::thread producer(rabbit::trim::producer_pe_task, std::ref(rp), std::ref(logger), fastqPool, std::ref(queue1));

  // consumer
  std::thread **consumer_threads = new std::thread* [consumer_num]; 
  if(rp.seqA.size()) 
  {
    // Ktrim
    for(int tn = 0; tn < consumer_num; tn++){
      consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_pe_task2, std::ref(rp), fastqPool, std::ref(queue1), wbDataPool, std::ref(wbQueue1), std::ref(wbQueue2), std::ref(statsArr[tn]), std::ref(trimmers)));
    }
  }
  else
  {
    // Trimmomatic
    for(int tn = 0; tn < consumer_num; tn++){
      consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_pe_task, std::ref(rp), fastqPool, std::ref(queue1), wbDataPool, std::ref(wbQueue1), std::ref(wbQueue2), std::ref(wbQueue3), std::ref(wbQueue4), std::ref(logQueue), std::ref(statsArr[tn]), std::ref(trimmers), tn));
    }
  }

  // writer 
  std::thread** writers = NULL;
  // pigzer
  std::thread** pigzers = NULL;

  if(rabbit::trim::util::endsWith(rp.output, ".gz") && rp.usePigz)
  {
    std::string out_file = rp.output.substr(0, rp.output.find(".gz"));
    
    if(rp.seqA.size()) // Ktrim
    {
      pigzers = new std::thread* [2];
      std::string out_file_1 = out_file + ".read1";
      std::string out_file_2 = out_file + ".read2";
      std::ofstream fout1(out_file_1.c_str());
      if(fout1.fail()){
        logger.errorln("Can not open file " + out_file_1); 
        exit(1);
      }
      std::ofstream fout2(out_file_2.c_str());
      if(fout2.fail()){
        logger.errorln("Can not open file " + out_file_2);
        exit(1);
      }

      std::pair<char*, int> pigzLast1;
      std::pair<char*, int> pigzLast2;
      pigzLast1.first = new char [ MEM_PER_CHUNK ];
      pigzLast1.second = 0;
      pigzLast2.first = new char [ MEM_PER_CHUNK ];
      pigzLast2.second = 0;


      pigzers[0] = new std::thread(rabbit::trim::pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue1), std::ref(pigzLast1), out_file_1);
      pigzers[1] = new std::thread(rabbit::trim::pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue2), std::ref(pigzLast2), out_file_2);
    }
    else // Trimmomatic
    {
      pigzers = new std::thread* [4];
      std::string read1_p_name = out_file + ".read1_p";
      std::string read2_p_name = out_file + ".read2_p";
      std::string read1_u_name = out_file + ".read1_u";
      std::string read2_u_name = out_file + ".read2_u";

      std::ofstream fout1(read1_p_name.c_str());
      if(fout1.fail()){
        logger.errorln("Can not open file " + read1_p_name);
        exit(1);
      }
      std::ofstream fout2(read2_p_name.c_str());
      if(fout2.fail()){
        logger.errorln("Can not open file " + read2_p_name);
        fout1.close();
        exit(1);
      }
      std::ofstream fout3(read1_u_name.c_str());
      if(fout3.fail()){
        logger.errorln("Can not open file " + read1_u_name);
        fout1.close(); fout2.close();
        exit(1);
      }
      std::ofstream fout4(read2_u_name.c_str());
      if(fout4.fail()){
        logger.errorln("Can not open file " + read2_u_name);
        fout1.close(); fout2.close(); fout3.close();
        exit(1);
      }

      std::pair<char*, int> pigzLast1;
      std::pair<char*, int> pigzLast2;
      std::pair<char*, int> pigzLast3;
      std::pair<char*, int> pigzLast4;
      pigzLast1.first = new char [ MEM_PER_CHUNK ];
      pigzLast1.second = 0;
      pigzLast2.first = new char [ MEM_PER_CHUNK ];
      pigzLast2.second = 0;
      pigzLast3.first = new char [ MEM_PER_CHUNK ];
      pigzLast3.second = 0;
      pigzLast4.first = new char [ MEM_PER_CHUNK ];
      pigzLast4.second = 0;


      pigzers[0] = new std::thread(rabbit::trim::pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue1), std::ref(pigzLast1), read1_p_name);
      pigzers[1] = new std::thread(rabbit::trim::pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue2), std::ref(pigzLast2), read2_p_name);
      pigzers[2] = new std::thread(rabbit::trim::pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue3), std::ref(pigzLast3), read1_u_name);
      pigzers[3] = new std::thread(rabbit::trim::pigzer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue4), std::ref(pigzLast4), read2_u_name);
    }
    
  }
  else
  {
    // don't use pigz 
    
    std::string out_file = rp.output; 
    std::string surfix = "";
    if(rabbit::trim::util::endsWith(rp.output, ".gz")){
      out_file = out_file.substr(0, out_file.size() - 3);
      surfix = ".gz";
    }
  
    if(rp.seqA.size()) // Ktrim
    {
      writers = new std::thread* [2];
      std::string out_file_1 = out_file + ".read1" + surfix;
      std::string out_file_2 = out_file + ".read2" + surfix;
      writers[0] = new std::thread(rabbit::trim::writer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue1), out_file_1, std::ref(logger));
      writers[1] = new std::thread(rabbit::trim::writer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue2), out_file_2, std::ref(logger));
    }
    else // Trimmomatic
    {
      writers = new std::thread* [4];
      std::string out_file_1 = out_file + ".read1_p" + surfix;
      std::string out_file_2 = out_file + ".read2_p" + surfix;
      std::string out_file_3 = out_file + ".read1_u" + surfix;
      std::string out_file_4 = out_file + ".read2_u" + surfix;
      writers[0] = new std::thread(rabbit::trim::writer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue1), out_file_1, std::ref(logger));
      writers[1] = new std::thread(rabbit::trim::writer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue2), out_file_2, std::ref(logger));
      writers[2] = new std::thread(rabbit::trim::writer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue3), out_file_3, std::ref(logger));
      writers[3] = new std::thread(rabbit::trim::writer_pe_task, std::ref(rp), wbDataPool, std::ref(wbQueue4), out_file_4, std::ref(logger));
    }
  }

  // release thread resource
  producer.join();
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn]->join();
  }

  if(pigzers)
  {
    if(rp.seqA.size())
    {
      pigzers[0]->join();
      pigzers[1]->join();
    }
    else
    {
      pigzers[0]->join();
      pigzers[1]->join();
      pigzers[2]->join();
      pigzers[3]->join();
    }
  }
  else
  {
    if(rp.seqA.size())
    {
      writers[0]->join();
      writers[1]->join();
    }
    else
    {
      writers[0]->join();
      writers[1]->join();
      writers[2]->join();
      writers[3]->join();
    }
  }

  if(trimlog_thread) trimlog_thread -> join();
  if(trimlog_pigz)   trimlog_pigz -> join();
  
  // write trim stats
  rabbit::log::TrimStat * totalStat = new rabbit::log::TrimStat(logger);
  totalStat -> merge(statsArr);
  if(rp.seqA.size()) totalStat -> print(rp.stats);
  else totalStat -> printPE(rp.stats);
  return 0;
}

// producer
int rabbit::trim::producer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataPairChunkQueue& dq){
  if(rp.forwardFiles.size() != rp.reverseFiles.size()) {
    logger.errorln("\033[1;31mError: Read1 and Read2 do not contain equal sized files!\033[0m\n");
    exit(0);
  }
  int totalFileCount = rp.forwardFiles.size();
  rabbit::int64 n_chunks = 0;
  for(int fileCnt = 0; fileCnt < totalFileCount; fileCnt++){
    rabbit::fq::FastqFileReader * fqFileReader;
    // 判断是否是gz文件 如果不一致需要退出终止
    bool file_is_gz = false;
    unsigned int i = rp.forwardFiles[fileCnt].size() - 3;
    const char * p = rp.forwardFiles[fileCnt].c_str();
    if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
      file_is_gz = true;
    }
    i = rp.reverseFiles[fileCnt].size() - 3;
    p = rp.reverseFiles[fileCnt].c_str();
    if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ){
      if(file_is_gz){
        fqFileReader = new rabbit::fq::FastqFileReader(rp.forwardFiles[fileCnt], *fastqPool, true, rp.reverseFiles[fileCnt]);
      }else{
        fprintf(stderr,"\033[1;31mError: file %s and %s must have same file type!\033[0m\n",rp.forwardFiles[fileCnt].c_str(),rp.reverseFiles[fileCnt].c_str());
      }
    }else{
      if(!file_is_gz){
        fqFileReader = new rabbit::fq::FastqFileReader(rp.forwardFiles[fileCnt], *fastqPool, false, rp.reverseFiles[fileCnt]);
      }else{
        fprintf(stderr,"\033[1;31mError: file %s and %s must have same file type!\033[0m\n",rp.forwardFiles[fileCnt].c_str(),rp.reverseFiles[fileCnt].c_str());
      }
    }

    while (true){
      rabbit::fq::FastqPairChunk *fqPairChunk = new rabbit::fq::FastqPairChunk; 
      fqPairChunk->chunk = fqFileReader->readNextPairChunk();
      if(fqPairChunk->chunk == NULL) break;
      dq.Push(n_chunks,fqPairChunk->chunk);
      n_chunks++;
    }
    delete fqFileReader;
  }
  dq.SetCompleted();
  std::cout<<"Input files have "<<n_chunks<<" chunks "<<std::endl;
  return 0;
}


// trimmomatic comusmer task
void rabbit::trim::consumer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue1, WriterBufferDataQueue& wbQueue2,  WriterBufferDataQueue& wbQueue3, WriterBufferDataQueue& wbQueue4, WriterBufferDataQueue& logQueue, rabbit::log::TrimStat& rstats, std::vector<rabbit::trim::Trimmer*>& trimmers, int threadId)
{
  rabbit::int64 chunk_id;
  rabbit::fq::FastqPairChunk* fqPairChunk = new rabbit::fq::FastqPairChunk;
  while(dq.Pop(chunk_id,fqPairChunk->chunk)){
    std::vector<neoReference> data;
    data.reserve(2e4);
    int loaded = rabbit::fq::chunkFormat(fqPairChunk->chunk->left_part , data, true);
    int loaded2  = rabbit::fq::chunkFormat(fqPairChunk->chunk->right_part, data, true);
    ASSERT(loaded == loaded2);
    rstats.readsInput += loaded;
    for(auto trimmer : trimmers){
      trimmer -> processRecords(data, threadId, true, false);
    }
    // prepare trim log data
    bool isLog = (bool)rp.trimLog.size();

    // copy data to WriterBuffer which is gotten from WriterBufferDataPool
    rabbit::trim::WriterBuffer* wb1 = NULL;
    rabbit::trim::WriterBuffer* wb2 = NULL;
    rabbit::trim::WriterBuffer* wb3 = NULL;
    rabbit::trim::WriterBuffer* wb4 = NULL;
    rabbit::trim::WriterBuffer* tb = NULL;

    wbDataPool -> Acquire(wb1);
    wbDataPool -> Acquire(wb2);
    wbDataPool -> Acquire(wb3);
    wbDataPool -> Acquire(wb4);
    
    uint32_t pos1 = 0; // d1_p
    uint32_t pos2 = 0; // d2_p
    uint32_t pos3 = 0; // d1_u
    uint32_t pos4 = 0; // d2_u
    
    for(int i = 0; i < loaded; i++)
    {
      auto& rec1 = data[i];
      auto& rec2 = data[i + loaded];
      if(rec1.lseq == 0 && rec2.lseq == 0) continue;
      if(rec1.lseq && rec2.lseq)
      {
          rstats.readsSurvivingBoth++;
          std::memcpy(wb1->data + pos1, (rec1.base + rec1.pname), rec1.lname);
          pos1 += rec1.lname;
          wb1->data[pos1++] = '\n';
          std::memcpy(wb1->data + pos1, (rec1.base + rec1.pseq), rec1.lseq);
          pos1 += rec1.lseq;
          wb1->data[pos1++] = '\n';
          std::memcpy(wb1->data + pos1, (rec1.base + rec1.pstrand), rec1.lstrand);
          pos1 += rec1.lstrand;
          wb1->data[pos1++] = '\n';
          std::memcpy(wb1->data + pos1, (rec1.base + rec1.pqual), rec1.lqual);
          pos1 += rec1.lqual;
          wb1->data[pos1++] = '\n';
          // reverse
          std::memcpy(wb2->data + pos2, (rec2.base + rec2.pname), rec2.lname);
          pos2 += rec2.lname;
          wb2->data[pos2++] = '\n';
          std::memcpy(wb2->data + pos2, (rec2.base + rec2.pseq), rec2.lseq);
          pos2 += rec2.lseq;
          wb2->data[pos2++] = '\n';
          std::memcpy(wb2->data + pos2, (rec2.base + rec2.pstrand), rec2.lstrand);
          pos2 += rec2.lstrand;
          wb2->data[pos2++] = '\n';
          std::memcpy(wb2->data + pos2, (rec2.base + rec2.pqual), rec2.lqual);
          pos2 += rec2.lqual;
          wb2->data[pos2++] = '\n';
      }
      else
      {
        if(rec2.lseq)
        {
          rstats.readsSurvivingReverse++;
          std::memcpy(wb4->data + pos4, (rec2.base + rec2.pname), rec2.lname);
          pos4 += rec2.lname;
          wb4->data[pos4++] = '\n';
          std::memcpy(wb4->data + pos4, (rec2.base + rec2.pseq), rec2.lseq);
          pos4 += rec2.lseq;
          wb4->data[pos4++] = '\n';
          std::memcpy(wb4->data + pos4, (rec2.base + rec2.pstrand), rec2.lstrand);
          pos4 += rec2.lstrand;
          wb4->data[pos4++] = '\n';
          std::memcpy(wb4->data + pos4, (rec2.base + rec2.pqual), rec2.lqual);
          pos4 += rec2.lqual;
          wb4->data[pos4++] = '\n';
        }
        else
        {
          rstats.readsSurvivingForward++;
          std::memcpy(wb3->data + pos3, (rec1.base + rec1.pname), rec1.lname);
          pos3 += rec1.lname;
          wb3->data[pos3++] = '\n';
          std::memcpy(wb3->data + pos3, (rec1.base + rec1.pseq), rec1.lseq);
          pos3 += rec1.lseq;
          wb3->data[pos3++] = '\n';
          std::memcpy(wb3->data + pos3, (rec1.base + rec1.pstrand), rec1.lstrand);
          pos3 += rec1.lstrand;
          wb3->data[pos3++] = '\n';
          std::memcpy(wb3->data + pos3, (rec1.base + rec1.pqual), rec1.lqual);
          pos3 += rec1.lqual;
          wb3->data[pos3++] = '\n';
        }
      }
    }
    wb1->data[pos1] = '\0';
    wb2->data[pos2] = '\0';
    wb3->data[pos3] = '\0';
    wb4->data[pos4] = '\0';
    
    wb1->size = pos1;
    wb2->size = pos2;
    wb3->size = pos3;
    wb4->size = pos4;

    wbQueue1.Push(chunk_id, wb1);
    wbQueue2.Push(chunk_id, wb2);
    wbQueue3.Push(chunk_id, wb3);
    wbQueue4.Push(chunk_id, wb4);
    

    if(isLog)
    {
      wbDataPool -> Acquire(tb);
      uint32_t pos_ = 0;
      for(int i = 0; i < loaded; i++)
      {
        auto& rec1 = data[i];
        auto& rec2 = data[i + loaded];

        // rec1 
        std::memcpy(tb->data + pos_, (rec1.base + (rec1.pname + 1)), rec1.lname - 1);
        pos_ += rec1.lname - 1;
        tb->data[pos_++] = ' ';
        int length = 0;
        int startPos = 0;
        int endPos = 0;
        int trimTail = 0;
        if(rec1.lseq != 0)
        {
          length = rec1.lseq;
          startPos = rec1.pseq - rec1.porigin;
          endPos = startPos + length;
          trimTail = rec1.lorigin - endPos;
        }
        std::string tmp = std::to_string(length) + " " + std::to_string(startPos) + " " + std::to_string(endPos) + " " + std::to_string(trimTail) + "\n";
        std::memcpy(tb->data + pos_, tmp.c_str(), tmp.size());
        pos_ += tmp.size();

        // rec2
        std::memcpy(tb->data + pos_, (rec2.base + (rec2.pname + 1)), rec2.lname - 1);
        pos_ += rec2.lname - 1;
        tb->data[pos_++] = ' ';
        length = 0;
        startPos = 0;
        endPos = 0;
        trimTail = 0;
        if(rec2.lseq != 0)
        {
          length = rec2.lseq;
          startPos = rec2.pseq - rec2.porigin;
          endPos = startPos + length;
          trimTail = rec2.lorigin - endPos;
        }
        std::string tmp2 = std::to_string(length) + " " + std::to_string(startPos) + " " + std::to_string(endPos) + " " + std::to_string(trimTail) + "\n";
        std::memcpy(tb->data + pos_, tmp2.c_str(), tmp2.size());
        pos_ += tmp2.size();
      }

      tb->data[pos_] = '\0';
      tb->size = pos_;
      ASSERT(pos_ > MEM_PER_CHUNK);
      logQueue.Push(chunk_id, tb);
    }

    fastqPool->Release(fqPairChunk->chunk->left_part);
    fastqPool->Release(fqPairChunk->chunk->right_part);
  }
  wbQueue1.SetCompleted();
  wbQueue2.SetCompleted();
  wbQueue3.SetCompleted();
  wbQueue4.SetCompleted();
  logQueue.SetCompleted();
}

void rabbit::trim::writer_pe_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::string out_file, rabbit::Logger& logger){
  if(rabbit::trim::util::endsWith(out_file, ".gz"))
  {
#ifdef USE_IGZIP
      // use igzip to compress
#include "igzip_lib.h"
      int level_size_[10] = {
#ifdef ISAL_DEF_LVL0_DEFAULT
        ISAL_DEF_LVL0_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL1_DEFAULT
        ISAL_DEF_LVL1_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL2_DEFAULT
        ISAL_DEF_LVL2_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL3_DEFAULT
        ISAL_DEF_LVL3_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL4_DEFAULT
        ISAL_DEF_LVL4_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL5_DEFAULT
        ISAL_DEF_LVL5_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL6_DEFAULT
        ISAL_DEF_LVL6_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL7_DEFAULT
        ISAL_DEF_LVL7_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL8_DEFAULT
        ISAL_DEF_LVL8_DEFAULT,
#else
        0,
#endif
#ifdef ISAL_DEF_LVL9_DEFAULT
        ISAL_DEF_LVL9_DEFAULT,
#else
        0
#endif
      };
        FILE* out = fopen(out_file.c_str(), "wb");
        if(!out)
        {
          logger.errorln("Failed to open file ' " + out_file + " '");
          exit(1);
        }
        fflush(0);
        struct isal_zstream stream;
        struct isal_gzip_header gz_hdr;
        isal_gzip_header_init(&gz_hdr);
        isal_deflate_init(&stream);
        stream.flush = NO_FLUSH;
        const int LEVEL = rp.compressLevel > 3 ? 0 : rp.compressLevel;
        
        stream.level = LEVEL;
        stream.level_buf_size = level_size_[LEVEL];
        stream.level_buf = (uint8_t*)malloc(level_size_[LEVEL]);
        stream.gzip_flag = IGZIP_GZIP;
        // isal_write_gzip_header(&stream, &gz_hdr);
        stream.end_of_stream = 0;

        rabbit::int64 chunk_id;
        WriterBuffer* wb;
        WriterBuffer* wb2;
        uint8_t* igzp_output_buf = new uint8_t[MEM_PER_CHUNK];

        wbQueue.Pop(chunk_id, wb);
        while(true)
        {
          if(wbQueue.Pop(chunk_id, wb2))
          {
            stream.end_of_stream = 0;

            stream.avail_in = wb -> size;
            stream.next_in = (uint8_t*)(wb->data);
            stream.avail_out = MEM_PER_CHUNK;
            stream.next_out = igzp_output_buf;
            isal_deflate(&stream);
            fwrite(igzp_output_buf, 1, MEM_PER_CHUNK - stream.avail_out, out);
            wbDataPool -> Release(wb);
            wb = wb2;
          }
          else{
            stream.end_of_stream = 1;
            stream.flush = FULL_FLUSH;
            stream.avail_in = wb -> size;
            stream.next_in = (uint8_t*)(wb->data);
            stream.avail_out = MEM_PER_CHUNK;
            stream.next_out = igzp_output_buf;
            isal_deflate(&stream);
            // isal_deflate_end(&stream);
            fwrite(igzp_output_buf, 1, MEM_PER_CHUNK - stream.avail_out, out);
            wbDataPool -> Release(wb);
            break;

          }
        }
        fclose(out);
        delete [] igzp_output_buf;


#else
#include "zlib.h"
        // use zlib
        gzFile gz_out_stream = gzopen(out_file.c_str(), "w");
        if(gz_out_stream == NULL)
        {
          logger.errorln("Failed to open file ' " + out_file + " ' by using gzip");
          exit(1);
      }
      int ret = gzsetparams(gz_out_stream, rp.compressLevel, Z_DEFAULT_STRATEGY);
      
      if (ret != Z_OK) {
        logger.errorln("Failed to set gzip parameters");
        gzclose(gz_out_stream);
        exit(1);
      }
      rabbit::int64 chunk_id;
      WriterBuffer* wb;
      while(wbQueue.Pop(chunk_id, wb)){
				ret = gzwrite(gz_out_stream, wb->data, wb->size);
        if (wb -> size && ret <= 0) {
        	logger.errorln("Failed to write data");
          gzclose(gz_out_stream);
          exit(1);
        }
        wbDataPool -> Release(wb);
      }
      gzclose(gz_out_stream);
      
#endif

  }
  else
  {
    std::ofstream fout(out_file.c_str());
    if(fout.fail())
    {
      logger.errorln("Can not open output file : " + out_file);
      exit(1);
    }
    rabbit::int64 chunk_id;
    WriterBuffer* wb;
    while(wbQueue.Pop(chunk_id,wb))
    {
      fout << wb->data;
      wbDataPool -> Release(wb);
      
    }
    fout.close();
  }
}


// Ktrim consumer task
void rabbit::trim::consumer_pe_task2(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue1, WriterBufferDataQueue& wbQueue2, rabbit::log::TrimStat& rstats, std::vector<rabbit::trim::Trimmer*>& trimmers)
{
  rabbit::int64 chunk_id;
  rabbit::fq::FastqPairChunk* fqPairChunk = new rabbit::fq::FastqPairChunk;
  while(dq.Pop(chunk_id,fqPairChunk->chunk)){
    // std::vector<Reference> data;
    std::vector<neoReference> data;
    data.reserve(2e4);
    int loaded = rabbit::fq::chunkFormat(fqPairChunk->chunk->left_part , data, true);
    int loaded2  = rabbit::fq::chunkFormat(fqPairChunk->chunk->right_part, data, true);
    ASSERT(loaded == loaded2);
    rstats.readsInput += loaded;
    for(auto trimmer : trimmers){
      trimmer -> processRecords(data, true, false);
    }

    // copy data to WriterBuffer which is gotten form WriterBufferDataPool
    WriterBuffer* wb1 = NULL;
    WriterBuffer* wb2 = NULL;
    
    wbDataPool -> Acquire(wb1);
    wbDataPool -> Acquire(wb2);
    uint64_t pos1 = 0; // d1_p
    uint64_t pos2 = 0; // d2_p
    
    for(int i = 0; i < loaded; i++)
    {
      auto& rec1 = data[i];
      auto& rec2 = data[i + loaded];
      if(rec1.lorigin == 1) rstats.realHit++;
      if(rec1.lorigin == 2) rstats.tailHit++;
      if(rec1.lorigin == 3) 
      {
        rstats.realHit++;
        rstats.dimer++;
      }
      if(rec1.lseq == 0)
      {
        rstats.readsDropped++;
        continue;
      }
      rstats.readsSurvivingBoth++;
      std::memcpy(wb1->data + pos1, (rec1.base + rec1.pname), rec1.lname);
      pos1 += rec1.lname;
      wb1->data[pos1++] = '\n';
      std::memcpy(wb1->data + pos1, (rec1.base + rec1.pseq), rec1.lseq);
      pos1 += rec1.lseq;
      wb1->data[pos1++] = '\n';
      std::memcpy(wb1->data + pos1, (rec1.base + rec1.pstrand), rec1.lstrand);
      pos1 += rec1.lstrand;
      wb1->data[pos1++] = '\n';
      std::memcpy(wb1->data + pos1, (rec1.base + rec1.pqual), rec1.lqual);
      pos1 += rec1.lqual;
      wb1->data[pos1++] = '\n';
      // reverse
      std::memcpy(wb2->data + pos2, (rec2.base + rec2.pname), rec2.lname);
      pos2 += rec2.lname;
      wb2->data[pos2++] = '\n';
      std::memcpy(wb2->data + pos2, (rec2.base + rec2.pseq), rec2.lseq);
      pos2 += rec2.lseq;
      wb2->data[pos2++] = '\n';
      std::memcpy(wb2->data + pos2, (rec2.base + rec2.pstrand), rec2.lstrand);
      pos2 += rec2.lstrand;
      wb2->data[pos2++] = '\n';
      std::memcpy(wb2->data + pos2, (rec2.base + rec2.pqual), rec2.lqual);
      pos2 += rec2.lqual;
      wb2->data[pos2++] = '\n';
    }
    wb1->data[pos1] = '\0';
    wb2->data[pos2] = '\0';
    wb1->size = pos1;
    wb2->size = pos1;

    fastqPool->Release(fqPairChunk->chunk->left_part);
    fastqPool->Release(fqPairChunk->chunk->right_part);
    wbQueue1.Push(chunk_id, wb1);
    wbQueue2.Push(chunk_id, wb2);
  }
  wbQueue1.SetCompleted();
  wbQueue2.SetCompleted();
}


void rabbit::trim::pigzer_pe_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::pair<char*, int>& pigzLast, std::string out_file)
  // outfile is the output file name without ".gz"
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

    infos[8] = new char[out_file.length() + 1];
    std::memcpy(infos[8], out_file.c_str(), out_file.length());
    infos[8][out_file.length()] = '\0';
    main_pigz(cnt, infos, wbDataPool, wbQueue, pigzLast);
    delete [] pigzLast.first;
    pigzLast.first = NULL;

  }
