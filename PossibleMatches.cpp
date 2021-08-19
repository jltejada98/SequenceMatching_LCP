//
// Created by Jose Luis Tejada on 9/8/21.
//

#include "PossibleMatches.h"

PossibleMatches::PossibleMatches() {

}

void PossibleMatches::InsertMatch(int &SAIndex, size_t &seqIndex) {
    //Todo consider using a mapping from SAIndex to SeqIndex (Consider slow iteration when transfering sets)
    SAIndices.push_back(SAIndex);
    seqIndices.push_back(seqIndex);
    uniqueSAIndices.insert(SAIndex); //Used for fast lookup when submatching
    uniqueSeqIndices.insert(seqIndex); //Only keeps unique seqIndices
}

size_t PossibleMatches::UniqueSetIndices() {
    return uniqueSeqIndices.size();
}

std::shared_ptr<std::vector<std::unordered_set<size_t>>> PossibleMatches::MatchesForMap() {
    size_t index;
    size_t SASize = SAIndices.size();
    std::vector<std::unordered_set<size_t>> matchIndexVector;
    std::shared_ptr<std::unordered_set<size_t>> unorderedSet;

    matchIndexVector.reserve(uniqueSeqIndices.size()); //Assume that transfer matches only called when indices included in all sequences
    for (index = 0; index < uniqueSeqIndices.size(); ++index) {
        unorderedSet = std::make_shared<std::unordered_set<size_t>>();
        matchIndexVector.push_back(*unorderedSet);
    }
    for(index = 0; index < SASize; ++index){
        matchIndexVector.at(seqIndices[index]).insert(SAIndices[index]);
    }

    return std::make_shared<std::vector<std::unordered_set<size_t>>>(matchIndexVector);
}

void PossibleMatches::AddSubMatch(std::vector<int> &SA, std::vector<size_t> &seq) {
    size_t SAsize = SA.size();
    for (size_t i = 0; i < SAsize; ++i) {
        if(!uniqueSAIndices.count(SA.at(i))){ //Check if Suffix Array Already exists to avoid duplication
            InsertMatch(SA.at(i), seq.at(i));
        }
    }
}

std::vector<int> &PossibleMatches::getSA() {
    return SAIndices;
}

std::vector<size_t> &PossibleMatches::getSeq() {
    return seqIndices;
}




