//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_SEQUENCERANGE_H
#define SEQUENCEMATCHING_LCP_SEQUENCERANGE_H

#include <unordered_set>

class SequenceRange {
private:
    size_t seqStart;
    size_t seqEnd;
    bool seqChecked = false;
    std::shared_ptr<std::unordered_set<size_t>> seqIndecies = nullptr;
public:
    explicit SequenceRange(size_t end, size_t start);
    bool getCheck() const;
    bool checkInsertIndex(size_t index);
    bool checkSetAllocation();
    void insertSet(std::shared_ptr<std::unordered_set<size_t>> seqInd);
    std::shared_ptr<std::unordered_set<size_t>> getSet();
    std::shared_ptr<std::unordered_set<size_t>> transferSet();
    void reset();
};


#endif //SEQUENCEMATCHING_LCP_SEQUENCERANGE_H
