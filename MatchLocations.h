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
    std::unordered_set<int> LCPIndecies;
    std::vector<std::unordered_set<int>> matchVector; //Vec Pos: SeqIndex, Unordered set: SAIndex
public:
    explicit MatchLocations(int &numSequences);
    explicit MatchLocations();
    void InsertUniqueSeq_LCPIndex(int &seqIndex, int &lcpIndex);
    void mergeUniqueSeq_LCPIndex(const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies,
                                 const std::shared_ptr<std::unordered_set<int>> &mergeLCPIndecies);
    int numSeqIncluded();
    std::shared_ptr<std::unordered_set<int>> getUniqueSeqIndices();
    std::shared_ptr<std::unordered_set<int>> getLCPIndecies();



    //To modify/Delete
    void InsertMatch(int &SAIndex, int &seqIndex);
    std::shared_ptr<std::vector<std::unordered_set<int>>> getMatchVector();
    void mergeMatches(const std::shared_ptr<std::vector<std::unordered_set<int>>> &matchVectorToMerge,
                      const std::shared_ptr<std::unordered_set<int>> &mergeSeqIndecies);
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
