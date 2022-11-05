#include "handler/HandlerPE.h"

using namespace rabbit::trim;

int rabbit::trim::process_pe(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
  int consumer_num = rp.threads - 3;
  rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
  rabbit::trim::FastqDataPairChunkQueue queue1(256,1);
  rabbit::trim::PEWriterDataQueue queue2(256,consumer_num);
  std::vector<rabbit::trim::TrimStat> statsArr(consumer_num);

  // PairingValidator
  rabbit::PairingValidator* pairingValidator;
  if(rp.validatePairing)
    pairingValidator = new rabbit::PairingValidator(logger);

  // trimmers
  TrimmerFactory *trimmerFactory = new TrimmerFactory(logger);
  std::vector<Trimmer*> trimmers;
  trimmerFactory -> makeTrimmers(rp.steps, rp.phred, trimmers);

  // trim log
  std::ofstream foutTrimLog(rp.trimLog);
  if(rp.trimLog.size()){
    if(foutTrimLog.fail()){
      logger.errorln("Can not open file " + rp.trimLog);
      return 105;
    }
  }
  if(rp.trimLog.size()) foutTrimLog.close();

  // producer
  std::thread producer(rabbit::trim::producer_pe_task, std::ref(rp), std::ref(logger), fastqPool, std::ref(queue1));

  // consumer
  std::thread **consumer_threads = new std::thread* [consumer_num]; 
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_pe_task, std::ref(rp), fastqPool, std::ref(queue1), std::ref(queue2), std::ref(statsArr[tn]), std::ref(trimmers)));
  }

  // writer 
  std::thread writer(rabbit::trim::writer_pe_task, std::ref(rp), std::ref(queue2), std::ref(logger));

  producer.join();
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn]->join();
  }
  writer.join();

  // if(isFileExist(fname)){
  // 	// delete file
  // 	if(remove(fname) != 0){
  // 		fprintf( stderr, "\033[1;34mError: remove file %s error!\033[0m\n",fname);
  // 	}
  // }

  // write trim stats
  rabbit::trim::TrimStat * totalStat = new rabbit::trim::TrimStat(logger);
  totalStat -> merge(statsArr);
  totalStat -> printPE(rp.stats);
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
        fqFileReader = new rabbit::fq::FastqFileReader(rp.forwardFiles[fileCnt], *fastqPool, rp.reverseFiles[fileCnt], true);
      }else{
        fprintf(stderr,"\033[1;31mError: file %s and %s must have same file type!\033[0m\n",rp.forwardFiles[fileCnt].c_str(),rp.reverseFiles[fileCnt].c_str());
      }
    }else{
      if(!file_is_gz){
        fqFileReader = new rabbit::fq::FastqFileReader(rp.forwardFiles[fileCnt], *fastqPool, rp.reverseFiles[fileCnt], false);
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


// comusmer task
void rabbit::trim::consumer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, PEWriterDataQueue& dq2, rabbit::trim::TrimStat& rstats,
    std::vector<rabbit::trim::Trimmer*>& trimmers)
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
    // copy data to WriterBuffer
    PEWriterBuffer *wb = new PEWriterBuffer((unsigned int)MEM_PER_CHUNK);
    uint64_t pos1 = 0; // d1_p
    uint64_t pos2 = 0; // d2_p
    uint64_t pos3 = 0; // d1_u
    uint64_t pos4 = 0; // d2_u

    
    for(int i = 0; i < loaded; i++)
    {
      auto& rec1 = data[i];
      auto& rec2 = data[i + loaded];
      if(rec1.lseq == 0 && rec2.lseq == 0) continue;
      if(rec1.lseq && rec2.lseq)
      {
          rstats.readsSurvivingBoth++;
          std::memcpy(wb->d1_p + pos1, (rec1.base + rec1.pname), rec1.lname);
          pos1 += rec1.lname;
          wb->d1_p[pos1++] = '\n';
          std::memcpy(wb->d1_p + pos1, (rec1.base + rec1.pseq), rec1.lseq);
          pos1 += rec1.lseq;
          wb->d1_p[pos1++] = '\n';
          std::memcpy(wb->d1_p + pos1, (rec1.base + rec1.pstrand), rec1.lstrand);
          pos1 += rec1.lstrand;
          wb->d1_p[pos1++] = '\n';
          std::memcpy(wb->d1_p + pos1, (rec1.base + rec1.pqual), rec1.lqual);
          pos1 += rec1.lqual;
          wb->d1_p[pos1++] = '\n';
          // reverse
          std::memcpy(wb->d2_p + pos2, (rec2.base + rec2.pname), rec2.lname);
          pos2 += rec2.lname;
          wb->d2_p[pos2++] = '\n';
          std::memcpy(wb->d2_p + pos2, (rec2.base + rec2.pseq), rec2.lseq);
          pos2 += rec2.lseq;
          wb->d2_p[pos2++] = '\n';
          std::memcpy(wb->d2_p + pos2, (rec2.base + rec2.pstrand), rec2.lstrand);
          pos2 += rec2.lstrand;
          wb->d2_p[pos2++] = '\n';
          std::memcpy(wb->d2_p + pos2, (rec2.base + rec2.pqual), rec2.lqual);
          pos2 += rec2.lqual;
          wb->d2_p[pos2++] = '\n';
      }
      else
      {
        if(rec2.lseq)
        {
          rstats.readsSurvivingReverse++;
          std::memcpy(wb->d2_u + pos4, (rec2.base + rec2.pname), rec2.lname);
          pos4 += rec2.lname;
          wb->d2_u[pos4++] = '\n';
          std::memcpy(wb->d2_u + pos4, (rec2.base + rec2.pseq), rec2.lseq);
          pos4 += rec2.lseq;
          wb->d2_u[pos4++] = '\n';
          std::memcpy(wb->d2_u + pos4, (rec2.base + rec2.pstrand), rec2.lstrand);
          pos4 += rec2.lstrand;
          wb->d2_u[pos4++] = '\n';
          std::memcpy(wb->d2_u + pos4, (rec2.base + rec2.pqual), rec2.lqual);
          pos4 += rec2.lqual;
          wb->d2_u[pos4++] = '\n';
        }
        else
        {
          rstats.readsSurvivingForward++;
          std::memcpy(wb->d1_u + pos3, (rec1.base + rec1.pname), rec1.lname);
          pos3 += rec1.lname;
          wb->d1_u[pos3++] = '\n';
          std::memcpy(wb->d1_u + pos3, (rec1.base + rec1.pseq), rec1.lseq);
          pos3 += rec1.lseq;
          wb->d1_u[pos3++] = '\n';
          std::memcpy(wb->d1_u + pos3, (rec1.base + rec1.pstrand), rec1.lstrand);
          pos3 += rec1.lstrand;
          wb->d1_u[pos3++] = '\n';
          std::memcpy(wb->d1_u + pos3, (rec1.base + rec1.pqual), rec1.lqual);
          pos3 += rec1.lqual;
          wb->d1_u[pos3++] = '\n';
        }
      }
    }
    wb->d1_p[pos1++] = '\0';
    wb->d2_p[pos2++] = '\0';
    wb->d1_u[pos3++] = '\0';
    wb->d2_u[pos4++] = '\0';
    fastqPool->Release(fqPairChunk->chunk->left_part);
    fastqPool->Release(fqPairChunk->chunk->right_part);
    dq2.Push(chunk_id, wb);
  }
  dq2.SetCompleted();
}

void rabbit::trim::writer_pe_task(rabbit::trim::RabbitTrimParam& rp, PEWriterDataQueue& dq2, rabbit::Logger& logger){
  std::ofstream fout1((rp.output + ".read1_p.fq").c_str());
  std::ofstream fout2((rp.output + ".read2_p.fq").c_str());
  std::ofstream fout3((rp.output + ".read1_u.fq").c_str());
  std::ofstream fout4((rp.output + ".read2_u.fq").c_str());
  if(fout1.fail() || fout2.fail() || fout3.fail() || fout4.fail()){
    logger.errorln("Can not open file");
    exit(1);
  }
  rabbit::int64 chunk_id;
  PEWriterBuffer *wb = new PEWriterBuffer; 
  while(dq2.Pop(chunk_id,wb)){
    fout1 << wb->d1_p;
    fout2 << wb->d2_p;
    fout3 << wb->d1_u;
    fout4 << wb->d2_u;
  }
  fout1.close();
  fout2.close();
}
