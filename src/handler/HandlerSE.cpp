#include "handler/HandlerSE.h"

using namespace rabbit::trim;

int rabbit::trim::process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
  // the number of consumer 
  int consumer_num = rp.threads - 3;

  rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
  rabbit::trim::FastqDataChunkQueue queue1(256,1);
  rabbit::trim::WriterDataQueue queue2(256, consumer_num);

  // Trim Stats for each thread
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
  std::thread producer(rabbit::trim::producer_se_task, std::ref(rp), std::ref(logger), fastqPool, std::ref(queue1));


  // consumer
  std::thread **consumer_threads = new std::thread* [consumer_num]; 
  for(int tn = 0; tn < consumer_num; tn++){
    // asm ("":::"memory");
    consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_se_task, std::ref(rp), fastqPool, std::ref(queue1), std::ref(queue2), std::ref(statsArr[tn]), std::ref(trimmers)));
  }

  // writer 
  std::thread writer(rabbit::trim::writer_se_task, std::ref(rp), std::ref(queue2), std::ref(logger));

  producer.join();
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn]->join();
  }
  writer.join();

  // writer->join();

  // if(isFileExist(fname)){
  // 	// delete file
  // 	if(remove(fname) != 0){
  // 		fprintf( stderr, "\033[1;34mError: remove file %s error!\033[0m\n",fname);
  // 	}
  // }

  // write trim stats
  rabbit::trim::TrimStat * totalStat = new rabbit::trim::TrimStat(logger);
  totalStat -> merge(statsArr);
  totalStat -> printSE(rp.stats);
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
    unsigned int i = inputFiles[fileCnt].size() - 3;
    const char * p = inputFiles[fileCnt].c_str();
    if( p[i]=='.' && p[i+1]=='g' && p[i+2]=='z' ) {
      file_is_gz = true;
    }
    if(file_is_gz){
      fqFileReader = new rabbit::fq::FastqFileReader(inputFiles[fileCnt],*fastqPool,"",true);
    }else{
      fqFileReader = new rabbit::fq::FastqFileReader(inputFiles[fileCnt],*fastqPool,"",false);
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
void rabbit::trim::consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, WriterDataQueue& dq2, rabbit::trim::TrimStat& rstats,
    const std::vector<rabbit::trim::Trimmer*>& trimmers)
{
  bool isReverse = rp.reverseFiles.size();
  rabbit::int64 chunk_id;
  rabbit::fq::FastqChunk* fqChunk = new rabbit::fq::FastqChunk;
  while(dq.Pop(chunk_id,fqChunk->chunk)){
    // std::vector<Reference> data;
    std::vector<neoReference> data;
    int loaded  = rabbit::fq::chunkFormat(fqChunk->chunk, data, true);
    rstats.readsInput += loaded;
    for(auto trimmer : trimmers){
      trimmer -> processRecords(data, false, isReverse);
    }

    // copy data to WriterBuffer
    WriterBuffer* wb = new WriterBuffer((unsigned int)MEM_PER_CHUNK);
    uint64_t pos = 0;
    for(auto rec : data){
      if(rec.lseq == 0) continue;
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
    wb->data[pos++] = '\0';
    fastqPool->Release(fqChunk->chunk);
    dq2.Push(chunk_id, wb);
  }
  dq2.SetCompleted();
}

void rabbit::trim::writer_se_task(rabbit::trim::RabbitTrimParam& rp,  WriterDataQueue& dq2, rabbit::Logger& logger){
  std::ofstream fout((rp.output + ".fq").c_str());
  if(fout.fail()){
    logger.errorln("Can not open file " + rp.output);
    exit(1);
  }
  rabbit::int64 chunk_id;
  WriterBuffer* wb = new WriterBuffer;
  while(dq2.Pop(chunk_id, wb)){
    // write to file
    fout << wb->data;
  }
  fout.close();
} 
