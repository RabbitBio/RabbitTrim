#include "core/HandlerSE.h"

using namespace rabbit::trim;

int rabbit::trim::process_se(rabbit::trim::RabbitTrimParam& rp, rabbit::Logger &logger) {
	rabbit::fq::FastqDataPool * fastqPool = new rabbit::fq::FastqDataPool(256,MEM_PER_CHUNK);
	rabbit::trim::FastqDataChunkQueue queue1(256,1);
    // the number of consumer 
	int consumer_num = rp.threads;
    
    // Trim Stats for each thread
    std::vector<rabbit::trim::TrimStat> statsArr;

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
    statsArr.emplace_back(); 
		consumer_threads[tn] = new std::thread(std::bind(rabbit::trim::consumer_se_task, std::ref(rp), fastqPool, std::ref(queue1), std::ref(statsArr[tn]), trimmers));
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
void rabbit::trim::consumer_se_task(rabbit::trim::RabbitTrimParam& rp, rabbit::fq::FastqDataPool *fastqPool, FastqDataChunkQueue &dq, rabbit::trim::TrimStat& rstats,
                                    std::vector<rabbit::trim::Trimmer*>& trimmers)
{
	int consumer_num = rp.threads;
	rabbit::int64 chunk_id;
	rabbit::fq::FastqChunk* fqChunk = new rabbit::fq::FastqChunk;
    std::vector<Reference> data;
	while(dq.Pop(chunk_id,fqChunk->chunk)){
        int loaded  = rabbit::fq::chunkFormat(fqChunk->chunk, data, true);
		fastqPool->Release(fqChunk->chunk);
        for(auto trimmer : trimmers){
            trimmer -> processRecords(data, false, rp.reverseFiles.size());
        }
        rstats.readsInput += loaded;
	}
}
