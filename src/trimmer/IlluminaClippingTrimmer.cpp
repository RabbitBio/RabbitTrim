#include "trimmer/IlluminaClippingTrimmer.h"

using namespace rabbit::trim;


IlluminaClippingTrimmer::IlluminaClippingTrimmer(rabbit::Logger& logger_, int phred_, std::string fastaAdapterFile, int seedMaxMiss_, int minPalindromeLikelihood_, int minSequenceLikelihood_, int minPrefix_, bool  palindromeKeepBoth_, int consumerNum_) : logger(logger_){
    phred = phred_;
    // adapter file 
    std::ifstream fastaAdapter(fastaAdapterFile.c_str(), std::ifstream::in);
    if(!fastaAdapter.is_open()){
        logger.errorln("Can not open fastaAdapterFile : " + fastaAdapterFile);
        exit(0);
    }
    
    seedMaxMiss = seedMaxMiss_;
    minPalindromeLikelihood = minPalindromeLikelihood_;
    minSequenceLikelihood = minSequenceLikelihood_;
    minPrefix = minPrefix_;
    palindromeKeepBoth = palindromeKeepBoth_;
    consumerNum = consumerNum_;
    
    minSequenceOverlap = (int)(minSequenceLikelihood / LOG10_4);
    minSequenceOverlap = minSequenceOverlap > 15 ? 15 : minSequenceOverlap; // TODO 15 是依据什么呢？ 

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
        if(reversePrefix.count(*iter)) prefixSet.insert(*iter);
   }
   for(auto iter = prefixSet.begin(); iter != prefixSet.end(); iter++){
        std::string forwardName = *iter + SUFFIX_F;
        std::string reverseName = *iter + SUFFIX_R;

        std::string forwardRec = forwardSeqMap[forwardName];
        std::string reverseRec = reverseSeqMap[reverseName];
        // delete same elements in forwardSeqMap and reverseSeqMap
        forwardSeqMap.erase(forwardName);
        reverseSeqMap.erase(reverseName);
        IlluminaPrefixPair* onePrefixPair = new IlluminaPrefixPair(forwardRec, reverseRec, logger_, phred, minPrefix, seedMaxMiss,  minPalindromeLikelihood,  palindromeKeepBoth, consumerNum_);
        prefixPairs.emplace_back(onePrefixPair);
   }
   
   // mapClippingSet
   // 01 forwardSeqMap 
   for(auto iter = forwardSeqMap.begin(); iter != forwardSeqMap.end(); iter++){
        std::set<std::string> uniqueSeq;
        std::string sequence =  iter -> second;
        if(uniqueSeq.count(sequence))
            logger.warningln("Skipping duplicate Clipping Sequence: '" + sequence + "'");
        else{
            uniqueSeq.insert(sequence);
            if(sequence.size() < 16){
                IlluminaClippingSeq* clippingSeq = new IlluminaShortClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                forwardSeqs.emplace_back(clippingSeq);
            }else{
                if(sequence.size() < 24){
                    IlluminaClippingSeq* clippingSeq = new IlluminaMediumClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                    forwardSeqs.emplace_back(clippingSeq);
                }else{
                    IlluminaClippingSeq* clippingSeq = new IlluminaLongClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                    forwardSeqs.emplace_back(clippingSeq);
                }
            }
        }
         
   }
   // 02 reverseSeqMap 
   for(auto iter = reverseSeqMap.begin(); iter != reverseSeqMap.end(); iter++){
        std::set<std::string> uniqueSeq;
        std::string sequence =  iter -> second;
        if(uniqueSeq.count(sequence))
            logger.warningln("Skipping duplicate Clipping Sequence: '" + sequence + "'");
        else{
            uniqueSeq.insert(sequence);
            if(sequence.size() < 16){
                IlluminaClippingSeq* clippingSeq = new IlluminaShortClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                reverseSeqs.emplace_back(clippingSeq);
            }else{
                if(sequence.size() < 24){
                    IlluminaClippingSeq* clippingSeq = new IlluminaMediumClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                    reverseSeqs.emplace_back(clippingSeq);
                }else{
                    IlluminaClippingSeq* clippingSeq = new IlluminaLongClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                    reverseSeqs.emplace_back(clippingSeq);
                }
            }
        }
   }
   // 03 commonSeqMap 
   for(auto iter = commonSeqMap.begin(); iter != commonSeqMap.end(); iter++){
        std::set<std::string> uniqueSeq;
        std::string sequence =  iter -> second;
        if(uniqueSeq.count(sequence))
            logger.warningln("Skipping duplicate Clipping Sequence: '" + sequence + "'");
        else{
            uniqueSeq.insert(sequence);
            if(sequence.size() < 16){
                IlluminaClippingSeq* clippingSeq = new IlluminaShortClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                commonSeqs.emplace_back(clippingSeq);
            }else{
                if(sequence.size() < 24){
                    IlluminaClippingSeq* clippingSeq = new IlluminaMediumClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                    commonSeqs.emplace_back(clippingSeq);
                }else{
                    IlluminaClippingSeq* clippingSeq = new IlluminaLongClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap, consumerNum_);
                    commonSeqs.emplace_back(clippingSeq);
                }
            }
        }
   }

   logger.infoln("ILLUMINACLIP: Using " + std::to_string(prefixPairs.size()) + " prefix pairs, " + std::to_string(commonSeqs.size()) + " forward/reverse sequences, " + std::to_string(forwardSeqs.size()) +  " forward only sequences, " + std::to_string(reverseSeqs.size()) + " reverse only sequences");

}

IlluminaClippingTrimmer::~IlluminaClippingTrimmer(){
  std::cout << "=============" << std::endl;
  std::cout << "~IlluminaClippingTrimmer()" << std::endl;
  std::cout << "=============" << std::endl;

}
// single end
void IlluminaClippingTrimmer::processOneRecord(Reference& rec){}
void IlluminaClippingTrimmer::processOneRecord(neoReference& rec){}
void IlluminaClippingTrimmer::processSingleRecord(Reference& rec, bool isReverse){
    int toKeepLength = rec.length; 
    if(!isReverse){
        for(auto& iter : forwardSeqs){
            int toKeep = iter -> readsSeqCompare(rec);
            toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
        }
    }
    else{
        for(auto& iter : reverseSeqs){
            int toKeep = iter -> readsSeqCompare(rec);
            toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
        }
    }
    
    // common
    for(auto& iter : commonSeqs){
        int toKeep = iter -> readsSeqCompare(rec);
        toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
    }

    if(toKeepLength > 0){
        rec.length = toKeepLength;
    }else{
        rec.length = 0; // when offset less than zero
    }
}
void IlluminaClippingTrimmer::processSingleRecord(neoReference& rec, bool isReverse){
    int toKeepLength = rec.lseq; 
    if(!isReverse){
        for(auto& iter : forwardSeqs){
            int toKeep = iter -> readsSeqCompare(rec);
            toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
        }
    }
    else{
        for(auto& iter : reverseSeqs){
            int toKeep = iter -> readsSeqCompare(rec);
            toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
        }
    }
    
    // common
    for(auto& iter : commonSeqs){
        int toKeep = iter -> readsSeqCompare(rec);
        toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
    }

    if(toKeepLength > 0){
        rec.lseq = toKeepLength;
    }else{
        rec.lseq = 0; // when offset less than zero
    }
    rec.lqual = rec.lseq;
}

void IlluminaClippingTrimmer::processSingleRecord(neoReference& rec, int threadId, bool isReverse){
    int toKeepLength = rec.lseq; 
    for(auto& iter : forwardSeqs){
      int toKeep = iter -> readsSeqCompare(rec, threadId);
      toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
    }

    // common
    for(auto& iter : commonSeqs){
        int toKeep = iter -> readsSeqCompare(rec, threadId);
        toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
    }

    if(toKeepLength > 0){
        rec.lseq = toKeepLength;
    }else{
        rec.lseq = 0; // when offset less than zero
    }
    rec.lqual = rec.lseq;
}


// pair end
// 返回应该保留的序列的长度
// 没找到时返回INT_MAX 1 << 30
void IlluminaClippingTrimmer::processPairRecord(Reference& rec1, Reference& rec2){
    int toKeepForward = rec1.length;
    int toKeepReverse = rec2.length;
    int toKeep;
    for(auto& iter : prefixPairs){
        toKeep = iter -> palindromeReadsCompare(rec1, rec2);
        toKeepForward = (toKeep < toKeepForward) ? toKeep : toKeepForward;
        if(palindromeKeepBoth){
            toKeepReverse = (toKeep < toKeepReverse) ? toKeep : toKeepReverse;
        }else{
            toKeepReverse = 0;
        }
    }

    assert(toKeepForward >= 0);
    if(toKeepForward > 0){
        for(auto& iter : forwardSeqs){
            toKeep = iter -> readsSeqCompare(rec1);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
        for(auto& iter : commonSeqs){
            toKeep = iter -> readsSeqCompare(rec1);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
    }

    assert(toKeepReverse >= 0);
    if(toKeepReverse > 0){
        for(auto& iter : reverseSeqs){
            toKeep = iter -> readsSeqCompare(rec2);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
        for(auto& iter : commonSeqs){
            toKeep = iter -> readsSeqCompare(rec2);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
    }

    if(toKeepForward > 0){
        rec1.length = toKeepForward;
    }else{
        rec1.length = 0;
    }
    
    if(toKeepReverse > 0){
        rec2.length = toKeepReverse;
    }else{
        rec2.length = 0;
    }
    
}

void IlluminaClippingTrimmer::processPairRecord(neoReference& rec1, neoReference& rec2){
    int toKeepForward = rec1.lseq;
    int toKeepReverse = rec2.lseq;
    int toKeep;
    for(auto& iter : prefixPairs){
        toKeep = iter -> palindromeReadsCompare(rec1, rec2);
        if(toKeep != (1 << 30)){
          toKeepForward = (toKeep < toKeepForward) ? toKeep : toKeepForward;
          if(palindromeKeepBoth){
            toKeepReverse = (toKeep < toKeepReverse) ? toKeep : toKeepReverse;
          }else{
            toKeepReverse = 0;
          }
        }
    }

    assert(toKeepForward >= 0);
    if(toKeepForward > 0){
        for(auto& iter : forwardSeqs){
            toKeep = iter -> readsSeqCompare(rec1);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
        for(auto& iter : commonSeqs){
            toKeep = iter -> readsSeqCompare(rec1);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
    }

    assert(toKeepReverse >= 0);
    if(toKeepReverse > 0){
        for(auto& iter : reverseSeqs){
            toKeep = iter -> readsSeqCompare(rec2);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
        for(auto& iter : commonSeqs){
            toKeep = iter -> readsSeqCompare(rec2);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
    }

    if(toKeepForward > 0){
        rec1.lseq = toKeepForward;
    }else{
        rec1.lseq = 0;
    }
    
    if(toKeepReverse > 0){
        rec2.lseq = toKeepReverse;
    }else{
        rec2.lseq = 0;
    }
    rec1.lqual = rec1.lseq;
    rec2.lqual = rec2.lseq;
    
}

void IlluminaClippingTrimmer::processPairRecord(neoReference& rec1, neoReference& rec2, int threadId){
    int toKeepForward = rec1.lseq;
    int toKeepReverse = rec2.lseq;
    int toKeep;
    for(auto& iter : prefixPairs){
        toKeep = iter -> palindromeReadsCompare(rec1, rec2, threadId);
        if(toKeep != (1 << 30)){
          toKeepForward = (toKeep < toKeepForward) ? toKeep : toKeepForward;
          if(palindromeKeepBoth){
            toKeepReverse = (toKeep < toKeepReverse) ? toKeep : toKeepReverse;
          }else{
            toKeepReverse = 0;
          }
        }
    }

    assert(toKeepForward >= 0);
    if(toKeepForward > 0){
        for(auto& iter : forwardSeqs){
            toKeep = iter -> readsSeqCompare(rec1, threadId);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
        for(auto& iter : commonSeqs){
            toKeep = iter -> readsSeqCompare(rec1, threadId);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
    }

    assert(toKeepReverse >= 0);
    if(toKeepReverse > 0){
        for(auto& iter : reverseSeqs){
            toKeep = iter -> readsSeqCompare(rec2, threadId);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
        for(auto& iter : commonSeqs){
            toKeep = iter -> readsSeqCompare(rec2, threadId);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
    }

    if(toKeepForward > 0){
        rec1.lseq = toKeepForward;
    }else{
        rec1.lseq = 0;
    }
    
    if(toKeepReverse > 0){
        rec2.lseq = toKeepReverse;
    }else{
        rec2.lseq = 0;
    }
    rec1.lqual = rec1.lseq;
    rec2.lqual = rec2.lseq;
    
}

void IlluminaClippingTrimmer::processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse){
    if(isPair){
        int n = recs.size() / 2;
        ASSERT(n * 2 == recs.size());
        for(int i = 0; i < n; i++){
            Reference& rec1 = recs[i];
            Reference& rec2 = recs[i + n];
            processPairRecord(rec1, rec2);
        }
    }
    else{
        for(Reference& rec : recs){
            processSingleRecord(rec, isReverse);
        }
        
    }
}

void IlluminaClippingTrimmer::processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse){
    if(isPair){
        int n = recs.size() / 2;
        ASSERT(n * 2 == recs.size());
        for(int i = 0; i < n; i++){
            neoReference& rec1 = recs[i];
            neoReference& rec2 = recs[i + n];
            processPairRecord(rec1, rec2);
        }
    }
    else{
      for(neoReference& rec : recs){
        processSingleRecord(rec, isReverse);
      }

      // remove virtual function
      // rabbit::trim::IlluminaShortClippingSeq* iter = dynamic_cast<rabbit::trim::IlluminaShortClippingSeq*>(forwardSeqs[0]);
      // for(neoReference& rec : recs){
      //   int toKeepLength = rec.lseq; 
      //   int toKeep = iter -> readsSeqCompare(rec);
      //   toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
      //   rec.lseq = std::max(0, toKeepLength);
      // }


      // if(!isReverse){
      //   for(auto& iter : forwardSeqs){
      //     for(neoReference& rec : recs){
      //       int toKeepLength = rec.lseq; 
      //       int toKeep = iter -> readsSeqCompare(rec);
      //       toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
      //       rec.lseq = std::max(0, toKeepLength);
      //       rec.lqual = rec.lseq;
      //     }

      //   }
      // }
      // else{
      //   for(auto& iter : reverseSeqs){
      //     for(neoReference& rec : recs){
      //       int toKeepLength = rec.lseq; 
      //       int toKeep = iter -> readsSeqCompare(rec);
      //       toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
      //       rec.lseq = std::max(0, toKeepLength);
      //       rec.lqual = rec.lseq;
      //     }
      //   }
      // }

      // // common
      // for(auto& iter : commonSeqs){
      //   for(neoReference& rec : recs){
      //     int toKeepLength = rec.lseq; 
      //     int toKeep = iter -> readsSeqCompare(rec);
      //     toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
      //     rec.lseq = std::max(0, toKeepLength);
      //     rec.lqual = rec.lseq;
      //   }
      // }

    }
}

void IlluminaClippingTrimmer::processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse){
    if(isPair){
        int n = recs.size() / 2;
        ASSERT(n * 2 == recs.size());
        for(int i = 0; i < n; i++){
            neoReference& rec1 = recs[i];
            neoReference& rec2 = recs[i + n];
            // processPairRecord(rec1, rec2);
            processPairRecord(rec1, rec2, threadId);
        }
    }
    else{
      for(neoReference& rec : recs){
        // processSingleRecord(rec, isReverse);
        processSingleRecord(rec, threadId, isReverse);
      }
    }
}

void IlluminaClippingTrimmer::printCnt(){
  logger.infoln("repeat offset nums : " + std::to_string(dynamic_cast<rabbit::trim::IlluminaShortClippingSeq*>(forwardSeqs[0]) ->cnt));
  logger.infoln("total offset nums : " + std::to_string(dynamic_cast<rabbit::trim::IlluminaShortClippingSeq*>(forwardSeqs[0]) ->total));
}
