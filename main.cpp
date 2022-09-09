#include "FileManipulation.h"
#include "MatchLocations.h"
#include "SequenceMatching.h"
#include "SuffixArray.hpp"
#include "LCPArrayKasai.hpp"
#include <iostream>
#include <vector>
#include <boost/unordered_map.hpp>

int main(int argc, const char *argv[]) {
    //CHECK ARGUMENTS
    if (argc < 5){
        std::cout << "Incorrect Arguments \n" <<std::endl;
        std::cout << "./Sequence_Matching <seq1.txt> <seq2.txt> ... <seqN.txt> <min match length> <max match length>" << std::endl;
        return EXIT_FAILURE;
    }
    int numSequences = argc-3;
    if (numSequences > 32){ //Todo Change Sentinel Characters to allow for more sequencecs to be inputted.
        std::cout << "Too Many Sequences \n" <<std::endl;
        std::cout << "Number of Sequences <= 32 \n" <<std::endl;
        return EXIT_FAILURE;
    }

    //READ SEQUENCES
    int index;
    int minimumMatchSize = (std::stoi(argv[argc-2]));
    int maximumMatchSize;
    if(strcmp(argv[argc-1], "INF") != 0){
        maximumMatchSize = (std::stoi(argv[argc-1]));
    }
    else{
        maximumMatchSize = INT_MAX;
    }

    if(minimumMatchSize > maximumMatchSize){
        std::cout << "Minimum Match Size <= Maximum Match Size \n" << std::endl;
        return EXIT_FAILURE;
    }

    int lengthSum = 0;
    std::vector<std::shared_ptr<std::string>> seqStringVector;
    std::vector<int> seqRangeVector;
    std::cout << "Reading Files..." << std::endl;
    for (index = 0; index < numSequences; ++index) {
        seqStringVector.push_back(Load_Sequence((char *) argv[index + 1]));
        if (seqStringVector[index] == nullptr){
            return EXIT_FAILURE;
        }
        seqRangeVector.push_back(lengthSum + static_cast<int>(seqStringVector[index]->length()));
        lengthSum += static_cast<int>(seqStringVector[index]->length()) + 1; //Add for sentinel
    }

    //COMBINE SEQUENCES AND ADD SENTINELS
    std::shared_ptr<std::string> seqStringCombined = Combine_Sequences(seqStringVector, numSequences);

    //CREATE SUFFIX ARRAY AND LCP ARRAYS
    std::cout << "Determining Suffix/LCP arrays..." << std::endl;
    SuffixArray<int> seqSuffixArray(*seqStringCombined); //Todo modify Sais.c to use int
    std::vector<int> SAVector = seqSuffixArray.GetSuffixArray();
    LCPArrayKasai<int> seqLCPArray(SAVector, *seqStringCombined);
    std::vector<int> LCPVector = seqLCPArray.GetLCPArray();

    //CREATE SUFFIX ARRAY INDEX TO SEQUENCE MAPPING
    std::shared_ptr<std::vector<int>> indexVector = Determine_Index_Mapping(SAVector, seqRangeVector);


    //DETERMINE MATCHES
    std::cout << "Determining Matches..." << std::endl;

    std::shared_ptr<boost::unordered_map<std::string, MatchLocations>> matchesMap = Determine_Possible_Matches_Parent(
            LCPVector,
            SAVector,
            *indexVector,
            minimumMatchSize,
            maximumMatchSize,
            numSequences,
            seqStringVector,
            seqRangeVector);


    //DETERMINE SIMILARITY METRICS
    std::cout << "Determining Similarity Metrics..." << std::endl;
    std::shared_ptr<std::vector<double>> similarityMetricVector = Determine_SimilarityMetrics(*matchesMap,
                                                                                              *seqStringCombined,
                                                                                              seqStringVector,
                                                                                              seqRangeVector,
                                                                                              numSequences);
    //Write to outfile
    std::cout << "Writing Results..." << std::endl;
    if(!Write_Matches(*matchesMap,*similarityMetricVector,numSequences, "Results.txt")){
        return EXIT_FAILURE;
    }

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}

