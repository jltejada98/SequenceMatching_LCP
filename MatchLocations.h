//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
#define SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H

#include <unordered_set>
#include <vector>

class MatchLocations {
private:
    std::vector<std::unordered_set<size_t>> matchVector; //Vec Pos: SeqIndex, Unordered set: SAIndex
public:
    explicit MatchLocations(int &numSequences);
    void InsertMatch(size_t SAIndex, size_t &seqIndex);
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> getMatchVector();
    void mergeMatches(const std::shared_ptr<std::vector<std::unordered_set<size_t>>> &matchVectorToMerge);
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
