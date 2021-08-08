#include "FileManipulation.h"
#include "SequenceRange.h"
#include "MatchLocations.h"
#include "SequenceMatching.h"
#include "SuffixArray.hpp"
#include "LCPArrayKasai.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>

int main(int argc, const char *argv[]) {
    //CHECK ARGUMENTS
    if (argc < 4){
        std::cout << "Incorrect Arguments \n" <<std::endl;
        std::cout << "./Sequence_Matching <seq1.txt> <seq2.txt> ... <seqN.txt> <min match length>" << std::endl;
        return EXIT_FAILURE;
    }
    size_t numSequences = argc-2;
    if (numSequences > 32){
        std::cout << "Too Many Sequences \n" <<std::endl;
        std::cout << "Number of Sequences <= 32 \n" <<std::endl;
        return EXIT_FAILURE;
    }

    //READ SEQUENCES
    size_t index;
    size_t minimumMatchSize = (std::stoi(argv[argc-1]));
    size_t previousSeqEnd=0;
    std::shared_ptr<std::string> seqStringArray[numSequences];
    std::shared_ptr<SequenceRange> seqRangeArray[numSequences];
    std::cout << "Reading Files..." << std::endl;
    for (index = 0; index < numSequences; ++index) {
        seqStringArray[index] = Load_Sequence((char *) argv[index + 1]);
        if (seqStringArray[index] == nullptr){
            return EXIT_FAILURE;
        }
        SequenceRange item(previousSeqEnd,previousSeqEnd + seqStringArray[index]->length());
        seqRangeArray[index] = std::make_shared<SequenceRange>(item);
        previousSeqEnd += seqStringArray[index]->length() + 1; //Add +1 to include sentinel after each sequence

    }

    //COMBINE SEQUENCES AND ADD SENTINELS
    std::shared_ptr<std::string> seqStringCombined = Combine_Sequences(seqStringArray, numSequences);

    //CREATE SUFFIX ARRAY AND LCP ARRAYS
    std::cout << "Determining Suffix/LCP arrays..." << std::endl;
    SuffixArray<int> seqSuffixArray(*seqStringCombined);
    std::vector<int> SAvector = seqSuffixArray.GetSuffixArray();
    LCPArrayKasai<int> seqLCPArray(SAvector, *seqStringCombined);
    std::vector<int>LCPVector = seqLCPArray.GetLCPArray();

    //DETERMINE MATCHES
    std::cout << "Determining Matches..." << std::endl;
    std::shared_ptr<std::unordered_map<std::string, MatchLocations>> matchesMap = Determine_Matches(LCPVector, SAvector, seqRangeArray, seqStringArray, minimumMatchSize, numSequences);

    //DETERMINE SUBMATCHING
    std::cout << "Determining Sub-Matches..." << std::endl;
    Determine_Submatching(*matchesMap, minimumMatchSize, numSequences);


    //DETERMINE SIMILARITY METRICS
    std::cout << "Determining Similarity Metrics..." << std::endl;
    std::shared_ptr<std::vector<double>> similarityMetricVector = Determine_SimilarityMetrics(*matchesMap, *seqStringCombined, seqStringArray, numSequences);


    //Write to outfile
    std::cout << "Writing Results..." << std::endl; //Todo Consider using
    if(!Write_Matches(*matchesMap,*similarityMetricVector,numSequences, "Results_SUB.txt")){
        return EXIT_FAILURE;
    }

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}

