//
// Created by Jose Luis Tejada on 9/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H
#define SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H


#include <unordered_map>
#include <vector>
#include <unordered_set>

class PossibleMatches {
private:
    std::unordered_map<int, size_t> SAIndexMap; //Key : SAIndex, Val: seqIndex
    std::unordered_set<size_t> uniqueSeqIndices;
public:
    PossibleMatches();
    void InsertMatch(int SAIndex, size_t &seqIndex);
    void AddSubMatch(std::unordered_map<int, size_t> & map);
    std::unordered_map<int, size_t> & getSAIndexMap();
    size_t UniqueSetIndices();
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> TransferMatches();
};


#endif //SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H
