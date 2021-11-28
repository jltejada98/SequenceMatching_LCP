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
    std::unordered_map<int, int> SAIndexMap; //Key : SAIndex, Val: seqIndex
    std::unordered_set<int> uniqueSeqIndices;
public:
    PossibleMatches();
    void InsertMatch(int SAIndex, int &seqIndex);
    void AddSubMatch(std::unordered_map<int, int> & map);
    std::unordered_map<int, int> & getSAIndexMap();
    int UniqueSetIndices();
    std::shared_ptr<std::vector<std::unordered_set<int>>> TransferMatches();
};


#endif //SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H
