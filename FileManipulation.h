//
// Created by Jose Luis Tejada on 16/7/21.
//

#ifndef SEQUENCEMATCHING_LCP_FILEMANIPULATION_H
#define SEQUENCEMATCHING_LCP_FILEMANIPULATION_H

#include <fstream>
#include <iostream>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include "MatchLocations.h"

std::shared_ptr<std::string> Load_Sequence(const char *filename);
std::shared_ptr<std::string> Combine_Sequences(const std::vector<std::shared_ptr<std::string>> &seqStringArray, int numSequences);
bool Write_Matches(boost::unordered_map<std::string, MatchLocations> &matchesMap,std::vector<double> &similarityMetrics,
                   int &numSequences, const std::string& outFilename);



#endif //SEQUENCEMATCHING_LCP_FILEMANIPULATION_H
