#include "FileManipulation.h"
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
    size_t lengthSum = 0;
    std::vector<std::shared_ptr<std::string>> seqStringVector;
    std::vector<size_t> seqRangeVector;
    std::cout << "Reading Files..." << std::endl;
    for (index = 0; index < numSequences; ++index) {
        seqStringVector.push_back(Load_Sequence((char *) argv[index + 1]));
        if (seqStringVector[index] == nullptr){
            return EXIT_FAILURE;
        }
        seqRangeVector.push_back(lengthSum + seqStringVector[index]->length());
        lengthSum = seqStringVector[index]->length() + 1; //Add for sentinel
    }

    //COMBINE SEQUENCES AND ADD SENTINELS
    std::shared_ptr<std::string> seqStringCombined = Combine_Sequences(seqStringVector, numSequences);

    //CREATE SUFFIX ARRAY AND LCP ARRAYS
    std::cout << "Determining Suffix/LCP arrays..." << std::endl;
    SuffixArray<int> seqSuffixArray(*seqStringCombined); //Todo modify Sais.c to use size_t
    std::vector<int> SAVector = seqSuffixArray.GetSuffixArray();
    LCPArrayKasai<int> seqLCPArray(SAVector, *seqStringCombined);
    std::vector<int> LCPVector = seqLCPArray.GetLCPArray();

    //Todo consider subtracting lengths of previous sequence and sentinel so that suffix array indecies are well positioned.

    //CREATE SUFFIX ARRAY INDEX TO SEQUENCE MAPPING
    std::shared_ptr<std::vector<size_t>> indexVector = Determine_Index_Mapping(SAVector, seqRangeVector);

    //DETERMINE MATCHES
    std::cout << "Determining Matches..." << std::endl;
    std::shared_ptr<std::unordered_map<std::string, MatchLocations>> matchesMap = Determine_Matches(LCPVector, SAVector, *indexVector, seqStringVector, minimumMatchSize, numSequences);


    //DETERMINE SIMILARITY METRICS
    std::cout << "Determining Similarity Metrics..." << std::endl;
    std::shared_ptr<std::vector<double>> similarityMetricVector = Determine_SimilarityMetrics(*matchesMap, *seqStringCombined, seqStringVector, numSequences);


    //Write to outfile
    std::cout << "Writing Results..." << std::endl;
    if(!Write_Matches(*matchesMap,*similarityMetricVector,numSequences, "Results_LCP.txt")){
        return EXIT_FAILURE;
    }

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}

