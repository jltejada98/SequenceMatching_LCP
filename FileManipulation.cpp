//
// Created by Jose Luis Tejada on 16/7/21.
//

#include "FileManipulation.h"

std::shared_ptr<std::string> Load_Sequence(const char *filename){
    std::ifstream File;
    std::string sequenceString;
    File.exceptions(::std::ios_base::failbit | ::std::ios_base::badbit);
    try{
        File.open(filename);
        File.seekg(0, std::ios::end);
        sequenceString.reserve(File.tellg());
        File.seekg(0, std::ios::beg);
        sequenceString.assign((std::istreambuf_iterator<char>(File)), std::istreambuf_iterator<char>());
        File.close();
    }
    catch (const std::ifstream::failure &e){
        std::cerr << "Exception" << e.what() << e.code() <<"opening/reading/closing file" <<std::endl;
        return nullptr;
    }
    catch (const std::exception &e) {
        std::cerr << "Exception" << e.what() <<std::endl;
        return nullptr;
    }

    return std::make_shared<std::string>(sequenceString);
}


std::shared_ptr<std::string> Combine_Sequences(const std::shared_ptr<std::string> *seqStringArray, size_t numSequences) {
    char sentinel = '!'; //Sentinel not included in language dictionary.
    size_t index;
    std::string combinedString;

    for(index = 0; index < numSequences; ++index){
        combinedString.append(*seqStringArray[index]);
        combinedString.push_back(sentinel);
        ++sentinel;
    }

    return std::make_shared<std::string>(combinedString);
}

bool Write_Matches(std::unordered_map<std::string, MatchLocations> &matchesMap,std::vector<double> &similarityMetrics,
                   size_t &numSequences, const std::string &outFilename){
    std::ofstream File;
    File.exceptions(~::std::ios_base::goodbit);
    std::vector<std::shared_ptr<std::unordered_set<size_t>>> matchIndeciesVector;

    try{
        File.open(outFilename); //Open for writing

        //Write Similarity Metrics
        size_t i;
        for (i = 0; i < numSequences; ++i) {
            File << similarityMetrics[i] << std::endl;
        }
        File << similarityMetrics[numSequences] << std::endl; //Combined metric

        //Write Matches
        for (auto &key : matchesMap){
            File << key.first << std::endl;
            matchIndeciesVector = *key.second.getMatchVector();
            for (auto &seqSet: matchIndeciesVector){
                for (auto &seqIndex : *seqSet){
                    File << seqIndex << ",";
                }
                File << std::endl;
            }
        }

        File.close();
    }
    catch (const std::ifstream::failure &e){
        std::cerr << "Exception" << e.what() << e.code() <<" opening/reading/closing file" <<std::endl;
        return false;
    }
    catch (const std::exception &e) {
        std::cerr << "Exception" << e.what() <<std::endl;
        return false;
    }

    return true;
}






