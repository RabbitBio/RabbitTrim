#ifndef SEED_CLIPPING_TRIMMER_H
#define SEED_CLIPPING_TRIMMER_H

#include "trimmer/Trimmer.h"
#include "Reference.h"
#include <vector>
#include <string>
#include <set>
#include <cstring>
#include <math.h>

namespace rabbit
{
    namespace trim
    {
        class SeedClippingTrimmer : public Trimmer
        {
        private:
            const int DIMER_INSERT = 1;
            int phred;
            double mismatch;
            int minQual;
            int window;
            int minLen;
            std::string seqA;
            std::string seqB;

        public:
            SeedClippingTrimmer(double mismatch_, std::string seqA_, std::string seqB_, int minLen_, int minQual_, int window_, int phred_);
            ~SeedClippingTrimmer();
            void processOneRecord(Reference& rec);
            void processRecords(std::vector<Reference>& recs, bool isPair, bool isReverse);
            void processSingleRecord(neoReference& rec);
            void processPairRecord(neoReference& rec1, neoReference& rec2);
            void processOneRecord(neoReference& rec);
            void processRecords(std::vector<neoReference>& recs, bool isPair, bool isReverse);
            void processRecords(std::vector<neoReference>& recs, int threadId, bool isPair, bool isReverse);
            inline bool is_revcomp(char a, char b);
        };
        
        inline bool SeedClippingTrimmer::is_revcomp(char a, char b)
        {
            switch (a)
            {
            case 'A': return b == 'T';
            case 'T': return b == 'A';
            case 'C': return b == 'G';
            case 'G': return b == 'C';
            default : return false;
            }
            return false;

        }
        
    } // namespace trim
    
} // namespace rabbit

#endif
