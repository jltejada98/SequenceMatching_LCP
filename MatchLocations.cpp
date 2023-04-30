//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

MatchLocations::MatchLocations() {
    uniqueSeqIndices.clear();
    LCPIndecies.clear();
}

MatchLocations::MatchLocations(int & numSequences) {
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

void MatchLocations::InsertUniqueSeq_LCPIndex(int &seqIndex, int &lcpIndex) {
    uniqueSeqIndices.insert(seqIndex);
    LCPIndecies.insert(lcpIndex);
}

int MatchLocations::numSeqIncluded() {
    return uniqueSeqIndices.size();
}

std::shared_ptr<std::unordered_set<int>> MatchLocations::getUniqueSeqIndices() {
    return std::make_shared<std::unordered_set<int>>(uniqueSeqIndices);
}

std::shared_ptr<std::unordered_set<int>> MatchLocations::getLCPIndecies() {
    return std::make_shared<std::unordered_set<int>>(LCPIndecies);
}

void MatchLocations::mergeUniqueSeq_LCPIndex(const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies,
                                             const std::shared_ptr<std::unordered_set<int>> &mergeLCPIndecies) {
    for (auto &mergeSeqIndex: *mergeSeqIndecies) { //Add the sequence indecies
        uniqueSeqIndices.insert(mergeSeqIndex);
    }
    for(auto &mergeLCPIndex: *mergeLCPIndecies){ //Add LCP Indecies.
        LCPIndecies.insert(mergeLCPIndex);
    }
}






void MatchLocations::InsertMatch(int &SAIndex, int &seqIndex) {
    matchVector.at(seqIndex).insert(SAIndex);
    uniqueSeqIndices.insert(seqIndex); //Only keeps unique seqIndices.
}


std::shared_ptr<std::vector<std::unordered_set<int>>> MatchLocations::getMatchVector() {
    return std::make_shared<std::vector<std::unordered_set<int>>>(matchVector);
}



void MatchLocations::mergeMatches(const std::shared_ptr<std::vector<std::unordered_set<int>>> &matchVectorToMerge,
                                  const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies) {
    for (auto &mergeSeqIndex: *mergeSeqIndecies) { //Add the sequence indecies
        uniqueSeqIndices.insert(mergeSeqIndex);
    }

    for (int i = 0; i < numSeqIncluded(); ++i) {
        matchVector[i].insert(matchVectorToMerge->at(i).begin(), matchVectorToMerge->at(i).end());
    }
    matchVectorToMerge->clear(); //After merging clear contents of merged set.
}





