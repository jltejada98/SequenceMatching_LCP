//
// Created by Jose Luis Tejada on 9/8/21.
//

#include "PossibleMatches.h"

PossibleMatches::PossibleMatches() {

}

void PossibleMatches::InsertMatch(int SAIndex, int &seqIndex) {
    SAIndexMap.insert(std::make_pair(SAIndex, seqIndex));
    uniqueSeqIndices.insert(seqIndex); //Only keeps unique seqIndices.
}

int PossibleMatches::UniqueSetIndices() {
    return uniqueSeqIndices.size();
}

std::shared_ptr<std::vector<std::unordered_set<int>>> PossibleMatches::TransferMatches() {
    int index;
    std::vector<std::unordered_set<int>> matchIndexVector;
    std::shared_ptr<std::unordered_set<int>> unorderedSet;
    matchIndexVector.reserve(uniqueSeqIndices.size()); //Assume that transfer matches only called when indices included in all sequences
    for (index = 0; index < uniqueSeqIndices.size(); ++index) {
        unorderedSet = std::make_shared<std::unordered_set<int>>();
        matchIndexVector.push_back(*unorderedSet);
    }
    for (auto & match: SAIndexMap){
        matchIndexVector.at(match.second).insert(match.first);
    }

    return std::make_shared<std::vector<std::unordered_set<int>>>(matchIndexVector);
}

void PossibleMatches::AddSubMatch(std::unordered_map<int, int> & map) {
    for(auto &match: map){
        InsertMatch(const_cast<int &>(match.first), match.second); //Inserts only unique elements.
    }
}

std::unordered_map<int, int> & PossibleMatches::getSAIndexMap() {
    return SAIndexMap;
}



