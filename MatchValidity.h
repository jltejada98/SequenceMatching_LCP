//
// Created by Jose Luis Tejada on 8/11/22.
//

#ifndef SEQUENCEMATCHING_CPP_MATCHVALIDITY_H
#define SEQUENCEMATCHING_CPP_MATCHVALIDITY_H

#include <unordered_set>

//Minimal Data Structure used to determine validity of matches
//i.e: If a given match is contained in all sequences

class MatchValidity{
private:
    std::unordered_set<int> uniqueSeqIndices;
public:
    MatchValidity();
    size_t numSeqIncluded();
    void InsertSeqIndex(int &seqIndex);
    void mergeSeqIndex(const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies);
    std::shared_ptr<std::unordered_set<int>> getUniqueSeqIndices();
};








#endif //SEQUENCEMATCHING_CPP_MATCHVALIDITY_H
