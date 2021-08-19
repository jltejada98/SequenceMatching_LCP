//
// Created by Jose Luis Tejada on 9/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H
#define SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H

#include <vector>
#include <unordered_set>

class PossibleMatches {
private:
    std::vector<int> SAIndices;
    std::vector<size_t> seqIndices;
    std::unordered_set<int> uniqueSAIndices;
    std::unordered_set<size_t> uniqueSeqIndices;
public:
    PossibleMatches();
    void InsertMatch(int &SAIndex, size_t &seqIndex);
    void AddSubMatch(std::vector<int> & SA, std::vector<size_t> & seq);
    std::vector<int> & getSA();
    std::vector<size_t> & getSeq();
    size_t UniqueSetIndices();
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> MatchesForMap();
};


#endif //SEQUENCEMATCHING_LCP_POSSIBLEMATCHES_H
