#include "core/HandlerPE.h"

using namespace rabbit::trim;

int rabbit::trim::process_pe(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
  rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
  rabbit::trim::FastqDataPairChunkQueue queue1(256,1);

  int consumer_num = rp.threads - 3;

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
    consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_pe_task, std::ref(rp), fastqPool, std::ref(queue1), std::ref(statsArr[tn]), std::ref(trimmers)));
  }

  // writer 
  // std::thread* writer = new std::thread(std::bind(&writer_pe_task, rp, &consumer_task_finished));

  producer.join();
  for(int tn = 0; tn < consumer_num; tn++){
    consumer_threads[tn]->join();
  }
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
void rabbit::trim::consumer_pe_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataPairChunkQueue &dq, rabbit::trim::TrimStat& rstats,
    std::vector<rabbit::trim::Trimmer*>& trimmers)
{
  rabbit::int64 chunk_id;
  rabbit::fq::FastqPairChunk* fqPairChunk = new rabbit::fq::FastqPairChunk;
  while(dq.Pop(chunk_id,fqPairChunk->chunk)){
    std::vector<Reference> data;
    data.reserve(1e4);
    int loaded  = rabbit::fq::chunkFormat(fqPairChunk->chunk->right_part, data, true);
    int loaded2 = rabbit::fq::chunkFormat(fqPairChunk->chunk->left_part , data, true);
    fastqPool->Release(fqPairChunk->chunk->left_part);
    fastqPool->Release(fqPairChunk->chunk->right_part);
    ASSERT(loaded == loaded2);
    for(auto trimmer : trimmers){
      trimmer -> processRecords(data, true, false);
    }
    rstats.readsInput += loaded;
    data.clear();
  }
}
