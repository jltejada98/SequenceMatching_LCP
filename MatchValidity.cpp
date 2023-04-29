//
// Created by Jose Luis Tejada on 8/11/22.
//

#include "MatchValidity.h"

MatchValidity::MatchValidity() {
    uniqueSeqIndices.clear();
}


size_t MatchValidity::numSeqIncluded() {
    return uniqueSeqIndices.size();
}

void MatchValidity::InsertSeqIndex(int &seqIndex) {
    uniqueSeqIndices.insert(seqIndex);
}

void MatchValidity::mergeSeqIndex(const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies) {
    for (auto &mergeSeqIndex: *mergeSeqIndecies) { //Add the sequence indecies
        uniqueSeqIndices.insert(mergeSeqIndex);
    }
}

std::shared_ptr<std::unordered_set<int>> MatchValidity::getUniqueSeqIndices() {
    return std::make_shared<std::unordered_set<int>>(uniqueSeqIndices);
}

