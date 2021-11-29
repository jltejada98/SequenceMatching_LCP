//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
#define SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H

#include <unordered_set>
#include <vector>

class MatchLocations {
private:
    std::vector<std::unordered_set<int>> matchVector; //Vec Pos: SeqIndex, Unordered set: SAIndex
    std::unordered_set<int> uniqueSeqIndices;
public:
    explicit MatchLocations(int &numSequences);
    void InsertMatch(int SAIndex, int &seqIndex);
    int numSeqIncluded();
    std::shared_ptr<std::vector<std::unordered_set<int>>> getMatchVector();
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
