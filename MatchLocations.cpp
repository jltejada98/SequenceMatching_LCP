//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

MatchLocations::MatchLocations(int & numSequences) {
    int index;
    std::shared_ptr<std::unordered_set<size_t>> unorderedSet;
    matchVector.reserve(numSequences);
    for (index = 0; index < numSequences; ++index) {
        unorderedSet = std::make_shared<std::unordered_set<size_t>>();
        matchVector.push_back(*unorderedSet);
    }
}

void MatchLocations::InsertMatch(size_t SAIndex, size_t &seqIndex) {
    matchVector[seqIndex].insert(SAIndex);
}

std::shared_ptr<std::vector<std::unordered_set<size_t>>> MatchLocations::getMatchVector() {
    return std::make_shared<std::vector<std::unordered_set<size_t>>>(matchVector);
}

void MatchLocations::mergeMatches(const std::shared_ptr<std::vector<std::unordered_set<size_t>>> &matchVectorToMerge) {
    matchVector = *matchVectorToMerge;

//    for (int i = 0; i < matchVectorToMerge->size(); ++i) {
//        matchVector[i].insert(matchVectorToMerge->at(i).begin(), matchVectorToMerge->at(i).end());
//    }
//    matchVectorToMerge->clear(); //After merging clear contents of merged set.
}
