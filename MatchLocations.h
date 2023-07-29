//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
#define SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H

#include <unordered_set>
#include <vector>

class MatchLocations {
private:
    std::unordered_set<int> uniqueSeqIndices;
    std::unordered_set<size_t> LCPIndecies; //Consider using size_t
    std::vector<std::unordered_set<int>> matchVector; //Vec Pos: SeqIndex, Unordered set: SAIndex
public:
    explicit MatchLocations(int numSequences);
    explicit MatchLocations();
    void InsertUniqueSeq_LCPIndex(int seqIndex, size_t lcpIndex);
    void MergeUniqueSeq_LCPIndex(const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies,
                                 const std::shared_ptr<std::unordered_set<size_t>> &mergeLCPIndecies);
    int NumSeqIncluded();
    std::shared_ptr<std::unordered_set<int>> GetUniqueSeqIndices();
    std::shared_ptr<std::unordered_set<size_t>> GetLCPIndecies();
    void InsertLocation(int locationIndex, int seqIndex);
    void MergeMatches(const std::shared_ptr<std::vector<std::unordered_set<int>>> &matchVectorToMerge);
    std::shared_ptr<std::vector<std::unordered_set<int>>> GetMatchVector();
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
