//
// Created by Jose Luis Tejada on 9/8/21.
//

#include "PossibleMatches.h"

PossibleMatches::PossibleMatches() {

}

void PossibleMatches::InsertMatch(int &SAIndex, size_t &seqIndex) {
    SAIndexMap.insert(std::make_pair(SAIndex, seqIndex));
    uniqueSeqIndices.insert(seqIndex); //Only keeps unique seqIndices.
}

size_t PossibleMatches::UniqueSetIndices() {
    return uniqueSeqIndices.size();
}

std::shared_ptr<std::vector<std::unordered_set<size_t>>> PossibleMatches::TransferMatches() {
    size_t index;
    std::vector<std::unordered_set<size_t>> matchIndexVector;
    std::shared_ptr<std::unordered_set<size_t>> unorderedSet;
    matchIndexVector.reserve(uniqueSeqIndices.size()); //Assume that transfer matches only called when indices included in all sequences
    for (index = 0; index < uniqueSeqIndices.size(); ++index) {
        unorderedSet = std::make_shared<std::unordered_set<size_t>>();
        matchIndexVector.push_back(*unorderedSet);
    }
    for (auto & match: SAIndexMap){
        matchIndexVector.at(match.second).insert(match.first);
    }

    return std::make_shared<std::vector<std::unordered_set<size_t>>>(matchIndexVector);
}

void PossibleMatches::AddSubMatch(std::unordered_map<int, size_t> & map) {
    for(auto &match: map){
        InsertMatch(const_cast<int &>(match.first), match.second); //Inserts only unique elements.
    }
}

std::unordered_map<int, size_t> & PossibleMatches::getSAIndexMap() {
    return SAIndexMap;
}



