//
// Created by Jose Luis Tejada on 5/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
#define SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <regex>

#include "MatchLocations.h"
#include "SequenceRange.h"

std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches(std::vector<int> &LCPVector,std::vector<int> &SAvector, std::shared_ptr<SequenceRange> *seqRangeArray,
                  std::shared_ptr<std::string> *seqStringArray,size_t &minimumMatchSize, size_t &numSequences);

std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap,std::string &seqStringCombined,
                            std::shared_ptr<std::string> *seqStringArray,
                            size_t &numSequences);

void Determine_Submatching(std::unordered_map<std::string, MatchLocations> &matchesMap, size_t &minimumMatchSize,
                           size_t &numSequences);


#endif //SEQUENCEMATCHING_LCP_SEQUENCEMATCHING_H
