//
// Created by Jose Luis Tejada on 16/7/21.
//

#ifndef SEQUENCEMATCHING_LCP_FILEMANIPULATION_H
#define SEQUENCEMATCHING_LCP_FILEMANIPULATION_H

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "MatchLocations.h"

std::shared_ptr<std::string> Load_Sequence(const char *filename);
std::shared_ptr<std::string> Combine_Sequences(const std::shared_ptr<std::string> *seqStringArray, size_t numSequences);
bool Write_Matches(std::unordered_map<std::string, MatchLocations> &matchesMap,std::vector<double> &similarityMetrics,
                   size_t &numSequences, const std::string& outFilename);



#endif //SEQUENCEMATCHING_LCP_FILEMANIPULATION_H
