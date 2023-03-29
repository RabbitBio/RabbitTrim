#include "handler/HandlerSE.h"
#include <stdio.h>
#include <stdlib.h>
#include "pigz.h"

using namespace rabbit::trim;

int rabbit::trim::process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
  // the number of consumer 
  int consumer_num = rp.threads;

  rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
  rabbit::trim::FastqDataChunkQueue queue1(256,1);
  rabbit::trim::FastqDataChunkQueue phredQueue(256,1);
  rabbit::trim::WriterBufferDataQueue wbQueue (256, consumer_num);
  rabbit::trim::WriterBufferDataQueue logQueue(256, consumer_num);
  
  // output DataPool
  int wbDataPoolSize = 256;
  if(rp.trimLog.size()) wbDataPoolSize <<= 1;
  WriterBufferDataPool* wbDataPool = new WriterBufferDataPool(wbDataPoolSize, MEM_PER_CHUNK);

  // Trim Stats for each thread
  std::vector<rabbit::log::TrimStat> statsArr(consumer_num);

  
  std::atomic_bool isDeterminedPhred(true);
  if(rp.phred == 0) isDeterminedPhred = false;

  // producer
  std::thread producer(rabbit::trim::producer_se_task, std::ref(rp), std::ref(logger), fastqPool, std::ref(queue1), std::ref(phredQueue), std::ref(isDeterminedPhred));

  // determine phred
  if(rp.phred == 0)
  {
    int cur_pre_read_count = 0;
    int* qualHistogram = new int[256];
    
    rabbit::int64 tmp_chunk_id;
    rabbit::fq::FastqDataChunk* fqDataChunk;

    while(cur_pre_read_count < rabbit::trim::PREREAD_COUNT && phredQueue.Pop(tmp_chunk_id,fqDataChunk))
    {
      uint64_t pos = 0;
      const char* base = (char*)fqDataChunk -> data.Pointer();
      const uint64_t chunk_size = fqDataChunk -> size + 1;
      
      while(pos <= chunk_size)
      {
        int cnt = 0;
        while(pos <= chunk_size && cnt < 3)
        {
          if(base[pos] == '\n') cnt++;
          pos++;
        }

        while(pos <= chunk_size)
        {
          if(base[pos] == '\n')
          {
            cur_pre_read_count++;
            if(cur_pre_read_count == rabbit::trim::PREREAD_COUNT)
            {
              pos = chunk_size + 1;
              break;
            }

            pos++;
            break;
          }
          qualHistogram[base[pos] - 0]++;
          pos++;
        }
      }
      fastqPool -> Release(fqDataChunk);
    }
    isDeterminedPhred = true;
    //  release remain data chunk in phredQueue
    if(!phredQueue.IsEmpty())
    {
      phredQueue.SetCompleted();
      while(phredQueue.Pop(tmp_chunk_id,fqDataChunk))
      {
        fastqPool -> Release(fqDataChunk);
      }
    }

    int phred33Total = 0;
    int phred64Total = 0;
    
    for(int i = 33; i <= 58; i++)
    {
      phred33Total += qualHistogram[i];
    }
    for(int i = 80; i <= 104; i++)
    {
      phred64Total += qualHistogram[i];
    }

    if(phred33Total == 0 && phred64Total > 0)
    {
      rp.phred = 64;
      logger.infoln("Quality encoding detected as phred64" );
    }
    else
    {
      if(phred33Total > 0 && phred64Total == 0)
      {
        rp.phred = 33;
        logger.infoln("Quality encoding detected as phred33" );
      }
      else
      {
        logger.errorln("Unable to detect quality encoding");
        logger.errorln("You can try to use -p or --phred to specify phred");
        exit(1);
      }
    }
    
  }

  // trimmers
  TrimmerFactory *trimmerFactory = new TrimmerFactory(logger);
  std::vector<Trimmer*> trimmers;
  trimmerFactory -> makeTrimmers(rp, rp.steps, trimmers);
  
  

  // consumer
  std::atomic_ullong atomic_next_id;
  std::atomic_init(&atomic_next_id, 0ULL);
  std::thread **consumer_threads = new std::thread* [consumer_num]; 
  for(int tn = 0; tn < consumer_num; tn++){
    // asm ("":::"memory");
    consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_se_task, std::ref(rp), fastqPool, std::ref(queue1), wbDataPool, std::ref(wbQueue), std::ref(logQueue), std::ref(statsArr[tn]), std::ref(trimmers), tn, std::ref(atomic_next_id)));
  }

  
  // write thread
  std::thread* writer = NULL;
  // pigzer
  std::thread* pigzer = NULL;
  std::pair<char*, int> pigzLast;
  if(rabbit::trim::util::endsWith(rp.output, ".gz") && rp.usePigz){
    std::string out_file = rp.output.substr(0, rp.output.find(".gz"));
    std::ofstream fout(out_file.c_str());
    if(fout.fail()){
      logger.errorln("Can not open file " + out_file);
      exit(1);
    }
    fout.close();
    pigzLast.first = new char [ MEM_PER_CHUNK ];
    pigzLast.second = 0;
    pigzer = new std::thread(std::bind(rabbit::trim::pigzer_se_task, std::ref(rp),  wbDataPool, std::ref(wbQueue), std::ref(pigzLast), out_file));
  }
  else{
    // writer 
    writer = new std::thread(std::bind(rabbit::trim::writer_se_task, std::ref(rp), wbDataPool, std::ref(wbQueue), rp.output, std::ref(logger)));
  }

  // trim log
  std::thread* trimlog_thread = NULL;
  std::thread* trimlog_pigz = NULL;
  if(rp.trimLog.size()){
    if(rp.usePigz && rabbit::trim::util::endsWith(rp.trimLog, ".gz"))
    {
      std::string out_file = rp.trimLog.substr(0, rp.trimLog.size() - 3);
      std::ofstream fout(out_file);
      if(fout.fail())
      {
        logger.errorln("Failed to open file : '" + out_file + " '");
      }
      fout.close();
      std::pair<char*, int> pigzLast;
      pigzLast.first = new char [MEM_PER_CHUNK];
      pigzLast.second = 0;
      trimlog_pigz = new std::thread(std::bind(pigzer_se_task, std::ref(rp), wbDataPool, std::ref(logQueue), std::ref(pigzLast), out_file));
    }
    else
    {
      trimlog_thread = new std::thread(std::bind(writer_se_task, std::ref(rp), wbDataPool, std::ref(logQueue), rp.trimLog, std::ref(logger)));
    }
  }


  producer.join();
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn]->join();
  }
  if(writer) writer -> join();
  else pigzer -> join();

  if(trimlog_thread) trimlog_thread -> join();
  if(trimlog_pigz) trimlog_pigz -> join();

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
int rabbit::trim::producer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger& logger, rabbit::fq::FastqDataPool* fastqPool, FastqDataChunkQueue& dq, FastqDataChunkQueue& phredQueue, std::atomic_bool& isDeterminedPhred)
{
  int totalFileCount = rp.forwardFiles.size();
  std::vector<std::string> inputFiles = rp.forwardFiles;
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
      
      if(!isDeterminedPhred)
      {
        rabbit::fq::FastqDataChunk* tmp;
        fastqPool -> Acquire(tmp);
        tmp -> size = fqChunk -> chunk -> size;
        std::memcpy(tmp->data.Pointer(), fqChunk -> chunk -> data.Pointer(), tmp -> size + 1);
        rabbit::int64 tmp_chunk_id = n_chunks;
        phredQueue.Push(tmp_chunk_id, tmp);
      }

      dq.Push(n_chunks, fqChunk->chunk);
      n_chunks++;
    }
    delete fqFileReader;
  }
  dq.SetCompleted();
  phredQueue.SetCompleted();
  return 0;
}

  // comusmer task
void rabbit::trim::consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, WriterBufferDataQueue& logQueue, rabbit::log::TrimStat& rstats, const std::vector<rabbit::trim::Trimmer*>& trimmers, int threadId, std::atomic_ullong& atomic_next_id)
{
    bool isLog = rp.trimLog.size();
    rabbit::int64 chunk_id;
    rabbit::fq::FastqChunk* fqChunk = new rabbit::fq::FastqChunk;
    while(dq.Pop(chunk_id,fqChunk->chunk))
    {
      std::vector<neoReference> data;
      data.reserve(1e4);
      int loaded  = rabbit::fq::chunkFormat(fqChunk->chunk, data, true);
      rstats.readsInput += loaded;

      for(auto trimmer : trimmers){
        trimmer -> processRecords(data, threadId, false, false);
      }

      // copy data to WriterBuffer 
      WriterBuffer* wb = NULL;
      WriterBuffer* tb = NULL;
      wbDataPool -> Acquire(wb);
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
      if(isLog)
      {
        wbDataPool -> Acquire(tb);
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
        tb->data[pos_] = '\0';
        tb->size = pos_;
        while(atomic_next_id != chunk_id);
        int log_chunk_id = chunk_id;
        logQueue.Push(log_chunk_id, tb);
      }
      while(atomic_next_id != chunk_id);
      wbQueue.Push(chunk_id, wb);
      atomic_next_id++;
      fastqPool->Release(fqChunk->chunk);
    }
    wbQueue.SetCompleted();
    logQueue.SetCompleted();
 }

void rabbit::trim::writer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::string out_file, rabbit::Logger& logger)
{
    if(rabbit::trim::util::endsWith(out_file, ".gz")){
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
    else{
      std::ofstream fout(out_file.c_str());
      if(fout.fail()){
        logger.errorln("Failed to open file ' " + out_file + " '");
        exit(1);
      }
      rabbit::int64 chunk_id;
      WriterBuffer* wb;
      while(wbQueue.Pop(chunk_id, wb)){
        // write to file
        fout << wb->data;
        wbDataPool -> Release(wb);
      }
      fout.close();
    }
} 

void rabbit::trim::pigzer_se_task(rabbit::trim::RabbitTrimParam& rp, WriterBufferDataPool* wbDataPool, WriterBufferDataQueue& wbQueue, std::pair<char*, int>& pigzLast, std::string out_file)
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
