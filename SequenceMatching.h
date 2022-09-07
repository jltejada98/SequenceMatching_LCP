//
// Created by Jose Luis Tejada on 5/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
#define SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H

#include <vector>
#include <thread>
#include <unordered_map>
#include <map>
#include <iostream>
#include <regex>

#include "MatchValidity.h"
#include "MatchLocations.h"

std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray);


std::shared_ptr<std::unordered_map<std::string, MatchLocations>> Determine_Possible_Matches_Parent(
        std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
        int maximumMatchSize, int numSequences, std::vector<std::shared_ptr<std::string>> &seqStringVector,
        std::vector<int> &seqRangeVector);

void
Determine_Valid_Matches_Child(std::unordered_map<std::string, MatchValidity> &matchesMap, std::vector<int> &LCPVector,
                              std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
                              int maximumMatchSize, std::vector<std::shared_ptr<std::string>> &seqStringVector,
                              std::vector<int> &seqRangeVector, size_t startIndex, size_t endIndex);


std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_StringPartitions(const std::string &key, const int &keyLen, const int &minLength, const int &maxLength);


void Determine_Match_Locations_Child(std::unordered_map<std::string, MatchLocations> &matchesMap,
                                     std::vector<std::string> &validMatches, std::vector<int> &indexVector,
                                     std::vector<std::shared_ptr<std::string>> &seqStringVector, int numSequences,
                                     size_t startIndex, size_t endIndex);

std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap, std::string &seqStringCombined,
                            std::vector<std::shared_ptr<std::string>> &seqStringVector,
                            std::vector<int> &seqRangeVector, int &numSequences);



#endif //SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
