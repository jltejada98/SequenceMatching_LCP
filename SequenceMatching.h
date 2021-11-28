//
// Created by Jose Luis Tejada on 5/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
#define SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H

#include <vector>
#include <unordered_map>
#include <map>
#include <iostream>

#include "MatchLocations.h"
#include "PossibleMatches.h"

std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray);

std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches(std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector,
                  std::vector<std::shared_ptr<std::string>> &seqStringVector, int &minimumMatchSize,
                  int &numSequences, std::string &seqStringCombined);

std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_Partitions(const std::string &key, const int &keyLen, const int &minLength,
                     std::shared_ptr<std::vector<int>> &partitionsShiftList);

std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap, std::string &seqStringCombined,
                            std::vector<std::shared_ptr<std::string>> &seqStringArray,
                            int &numSequences);


#endif //SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
