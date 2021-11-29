//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

MatchLocations::MatchLocations(int & numSequences) {
    int index;
    std::shared_ptr<std::unordered_set<int>> unorderedSet;
    matchVector.reserve(numSequences);
    for (index = 0; index < numSequences; ++index) {
        unorderedSet = std::make_shared<std::unordered_set<int>>();
        matchVector.push_back(*unorderedSet);
    }
    uniqueSeqIndices.clear();
}

void MatchLocations::InsertMatch(int SAIndex, int &seqIndex) {
    matchVector.at(seqIndex).insert(SAIndex);
    uniqueSeqIndices.insert(seqIndex); //Only keeps unique seqIndices.
}

int MatchLocations::numSeqIncluded() {
    return uniqueSeqIndices.size();
}

std::shared_ptr<std::vector<std::unordered_set<int>>> MatchLocations::getMatchVector() {
    return std::make_shared<std::vector<std::unordered_set<int>>>(matchVector);
}
