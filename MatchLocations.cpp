//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

MatchLocations::MatchLocations() {
    uniqueSeqIndices.clear();
    LCPIndecies.clear();
}

MatchLocations::MatchLocations(int numSequences) {
    int index;
    std::shared_ptr<std::unordered_set<int>> unorderedSet;
    matchVector.reserve(numSequences);
    for (index = 0; index < numSequences; ++index) {
        unorderedSet = std::make_shared<std::unordered_set<int>>();
        matchVector.push_back(*unorderedSet);
    }
    uniqueSeqIndices.clear();
    LCPIndecies.clear();
}

void MatchLocations::InsertUniqueSeq(int seqIndex) {
    uniqueSeqIndices.insert(seqIndex);
}

void MatchLocations::Insert_LCPIndex(size_t lcpIndex) {
    LCPIndecies.insert(lcpIndex);
}

int MatchLocations::NumSeqIncluded() {
    return uniqueSeqIndices.size();
}

std::shared_ptr<std::unordered_set<int>> MatchLocations::GetUniqueSeqIndices() {
    return std::make_shared<std::unordered_set<int>>(uniqueSeqIndices);
}

std::shared_ptr<std::unordered_set<size_t>> MatchLocations::GetLCPIndecies() {
    return std::make_shared<std::unordered_set<size_t>>(LCPIndecies);
}

void MatchLocations::MergeUniqueSeq_LCPIndex(const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies,
                                             const std::shared_ptr<std::unordered_set<size_t>> &mergeLCPIndecies) {
    for (auto &mergeSeqIndex: *mergeSeqIndecies) { //Add the sequence indecies
        uniqueSeqIndices.insert(mergeSeqIndex);
    }
    for(auto &mergeLCPIndex: *mergeLCPIndecies){ //Add LCP Indecies.
        LCPIndecies.insert(mergeLCPIndex);
    }
}

void MatchLocations::InsertLocation(int locationIndex, int seqIndex) {
    matchVector[seqIndex].insert(locationIndex);
}

void MatchLocations::MergeMatches(const std::shared_ptr<std::vector<std::unordered_set<int>>> &matchVectorToMerge) {
    for (int i = 0; i < NumSeqIncluded(); ++i) {
        matchVector[i].insert(matchVectorToMerge->at(i).begin(), matchVectorToMerge->at(i).end());
    }
    matchVectorToMerge->clear();
}

std::shared_ptr<std::vector<std::unordered_set<int>>> MatchLocations::GetMatchVector() {
    return std::make_shared<std::vector<std::unordered_set<int>>>(matchVector);
}





