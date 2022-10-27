#include "IlluminaClippingTrimmer.h"

using namespace rabbit;


IlluminaClippingTrimmer::IlluminaClippingTrimmer(rabbit::Logger& logger_, int phred_, std::string fastaAdapterFile, int seedMaxMiss_, int minPalindromeLikelihood_, int minSequenceLikelihood_, int minPrefix_ = 1, int palindromeKeepBoth_ = false) : logger(logger_){
    phred = phred_;
    // adapter file 
    std::ifstream fastaAdapter(fastaAdapterFile.to_str(), std::ifstream::in);
    if(!fastaAdapter.is_open()){
        logger.errorln("Can not open fastaAdapterFile : " + fastaAdapterFile);
        exit(0);
    }
    
    seedMaxMiss = seedMaxMiss_;
    minPalindromeLikelihood = minPalindromeLikelihood_;
    minSequenceLikelihood = minSequenceLikelihood_;
    minPrefix = minPrefix_;
    palindromeKeepBoth = palindromeKeepBoth_;
    
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
        ClearHeadTailSpace(curLine);
        while(curLine.size() && curLine.at(0) != '>'){
            if(!std::getline(fastaAdapter, curLine)) break;
            ClearHeadTailSpace(curLine);
        }
        if(curLine.size() && curLine.at(0) == '>'){
            std::string fullName = curLine.substr(1);
            std:vector<std::string> tokens = split(fullName, "[\\| ]");
            name = tokens[0];
            
            sequence = "";
            while(std::getline(fastaAdapter, curLine)){
                ClearHeadTailSpace(curLine);
                if(curLine.size() == 0) continue;
                if(curLine.at(0) == '>') break;
                if(curLine.at(0) == ';') continue;
                sequence += curLine;
            }
        }
        if(sequence.size() == 0) continue;

        if(endsWith(name, SUFFIX_F)){
            forwardSeqMap.insert(std::make_pair(name, sequence));
            if(startsWith(name, PREFIX)){
                forwardPrefix.insert(name.substr(0, name.size() - SUFFIX_F.size()));
            }
        }
        else{
            if(endsWith(name, SUFFIX_R)){
                reverseSeqMap.insert(std::make_pair(name, sequence));
                if(startsWith(name, PREFIX)){
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
        IlluminaPrefixPair* onePrefixPair = new IlluminaPrefixPair(forwardRec, reverseRec, logger_, phred, minPrefix, seedMaxMiss,  minPalindromeLikelihood,  palindromeKeepBoth);
        prefixPairs.emplace_back(*onePrefixPair);
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
                IlluminaClippingSeq* clippingSeq = new IlluminaShortClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                forwardSeqs.emplace_back(*clippingSeq);
            }else{
                if(sequence.size() < 24){
                    IlluminaClippingSeq* clippingSeq = new IlluminaMediumClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                    forwardSeqs.emplace_back(*clippingSeq);
                }else{
                    IlluminaClippingSeq* clippingSeq = new IlluminaLongClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                    forwardSeqs.emplace_back(*clippingSeq);
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
                IlluminaClippingSeq* clippingSeq = new IlluminaShortClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                reverseSeqs.emplace_back(*clippingSeq);
            }else{
                if(sequence.size() < 24){
                    IlluminaClippingSeq* clippingSeq = new IlluminaMediumClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                    reverseSeqs.emplace_back(*clippingSeq);
                }else{
                    IlluminaClippingSeq* clippingSeq = new IlluminaLongClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                    reverseSeqs.emplace_back(*clippingSeq);
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
                IlluminaClippingSeq* clippingSeq = new IlluminaShortClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                commonSeqs.emplace_back(*clippingSeq);
            }else{
                if(sequence.size() < 24){
                    IlluminaClippingSeq* clippingSeq = new IlluminaMediumClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                    commonSeqs.emplace_back(*clippingSeq);
                }else{
                    IlluminaClippingSeq* clippingSeq = new IlluminaLongClippingSeq(logger, phred, sequence, seedMaxMiss,  minSequenceLikelihood,  minSequenceOverlap);
                    commonSeqs.emplace_back(*clippingSeq);
                }
            }
        }
   }

   logger.infoln("ILLUMINACLIP: Using " + prefixPairs.size() + " prefix pairs, " + commonSeq.size() + " forward/reverse sequences, " + forwardSeqs.size() +  " forward only sequences, " + reverseSeqs.size() + " reverse only sequences");

}
// single end
void IlluminaClippingTrimmer::processOneRecord(Reference& rec){}
void IlluminaClippingTrimmer::processSingleRecord(Reference& rec, bool isReverse = false){
    int toKeepLength = rec.length; 
    if(!isReverse){
        for(auto iter = forwardSeqs.begin(); iter != forwardSeqs.end(); iter++){
            int toKeep = iter -> readsSeqCompare(rec);
            toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
        }
    }
    else{
        for(auto iter = reverseSeqs.begin(); iter != reverseSeqs.end(); iter++){
            int toKeep = iter -> readsSeqCompare(rec);
            toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
        }
    }
    
    // common
    for(auto iter = commonSeqs.begin(); iter != commonSeqs.end(); iter++){
        int toKeep = iter -> readsSeqCompare(rec);
        toKeepLength = toKeep < toKeepLength ? toKeep : toKeepLength;
    }

    if(toKeepLength > 0){
        rec.length = toKeepLength;
    }else{
        rec.length = 0; // when offset less than zero
    }
}


// pair end
// 返回应该保留的序列的长度
// 没找到时返回INT_MAX
void IlluminaClippingTrimmer::processPairRecord(Reference& rec1, Reference& rec2){
    int toKeepForward = rec1.length;
    int toKeepReverse = rec2.length;
    for(auto iter = prefixPairs.begin(); iter != prefixPairs.end(); iter++){
        int toKeep = iter -> palindromeReadsCompare(rec1, rec2);
        toKeepForward = (toKeep < toKeepForward) ? toKeep : toKeepForward;
        if(palindromeKeepBoth){
            toKeepReverse = (toKeep < toKeepReverse) ? toKeep : toKeepReverse;
        }else{
            toKeepReverse = 0;
        }
    }

    assert(toKeepForward >= 0);
    if(toKeepForward > 0){
        for(auto iter = forwardSeqs.begin(); iter != forwardSeqs.end(); iter++){
            toKeep = iter -> readsSeqCompare(rec);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
        for(auto iter = commonSeqs.begin(); iter != commonSeqs.end(); iter++){
            toKeep = iter -> readsSeqCompare(rec);
            toKeepForward = toKeep < toKeepForward ? toKeep : toKeepForward;
        }
    }

    assert(toKeepReverse >= 0);
    if(toKeepReverse > 0){
        for(auto iter = reverseSeqs.begin(); iter != reverseSeqs.end(); iter++){
            toKeep = iter -> readsSeqCompare(rec);
            toKeepReverse = toKeep < toKeepReverse ? toKeep : toKeepReverse;
        }
        for(auto iter = commonSeqs.begin(); iter != commonSeqs.end(); iter++){
            toKeep = iter -> readsSeqCompare(rec);
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


void IlluminaClippingTrimmer::processRecords(std::vector<Reference&> recs, bool isPair, bool isReverse){
    if(isPair){
        int n = recs.size();
        ASSERT(n % 2 == 0);
        for(int i = 0; i < n; i += 2){
            Reference& rec1 = recs[i];
            Reference& rec2 = recs[i + 1];
            processPairRecord(rec1, rec2);
        }
    }
    else{
        for(Reference& rec : recs){
            processSingleRecord(rec, isReverse);
        }
        
    }
}