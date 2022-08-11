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

#include "MatchLocations.h"

std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray);

std::shared_ptr<std::unordered_map<std::string, MatchLocations>> Determine_Matches_Parent(
        std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
        int maximumMatchSize, int numSequences, std::vector<std::shared_ptr<std::string>> &seqStringVector,
        std::vector<int> &seqRangeVector);

void
Determine_Matches_Child(std::unordered_map<std::string, MatchLocations> &matchesMap, std::vector<int> &LCPVector,
                        std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
                        int maximumMatchSize, int numSequences,
                        std::vector<std::shared_ptr<std::string>> &seqStringVector,
                        std::vector<int> &seqRangeVector, size_t startIndex, size_t endIndex);

std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_Partitions(const std::string &key, const int &keyLen, const int &minLength, const int &maxLength,
                     std::shared_ptr<std::vector<int>> &partitionsShiftList);

std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap, std::string &seqStringCombined,
                            std::vector<std::shared_ptr<std::string>> &seqStringVector,
                            std::vector<int> &seqRangeVector, int &numSequences);


//Added functions for testing
std::shared_ptr<std::unordered_map<std::string, MatchLocations>> Determine_Possible_Matches_Parent(
        std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
        int maximumMatchSize, int numSequences, std::vector<std::shared_ptr<std::string>> &seqStringVector,
        std::vector<int> &seqRangeVector);

void Determine_Possible_Matches_Child(std::unordered_map<std::string, MatchLocations> &matchesMap, std::vector<int> &LCPVector,
                                      std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
                                      int maximumMatchSize, int numSequences, std::vector<std::shared_ptr<std::string>> &seqStringVector,
                                      std::vector<int> &seqRangeVector, size_t startIndex, size_t endIndex);


std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_StringPartitions(const std::string &key, const int &keyLen, const int &minLength, const int &maxLength);


#endif //SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
