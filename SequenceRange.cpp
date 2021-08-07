//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "SequenceRange.h"
#include <utility>

SequenceRange::SequenceRange(size_t start,size_t end) {
    seqStart = start;
    seqEnd = end;
}

void SequenceRange::insertSet(std::shared_ptr<std::unordered_set<size_t>> seqInd) {
    seqIndecies = std::move(seqInd);
}

bool SequenceRange::checkInsertIndex(size_t index) {
    if (index < seqEnd){
        seqChecked = true;
        seqIndecies->insert(index-seqStart); //Shift to original positions
        return true;
    }
    return false;
}

bool SequenceRange::getCheck() const {
    return seqChecked;
}

void SequenceRange::reset() {
    seqChecked = false;
    if(seqIndecies != nullptr){
        seqIndecies->clear(); //Clear values, use allocated sets for next match checking.
    }
}

std::shared_ptr<std::unordered_set<size_t>> SequenceRange::transferSet() {
    std::shared_ptr<std::unordered_set<size_t>> passSet = seqIndecies;
    seqIndecies = nullptr;
    seqChecked = false; //Reset checks for next set allocation.
    return passSet;
}

bool SequenceRange::checkSetAllocation() {
    if(seqIndecies == nullptr){
        return false;
    }
    return true;
}

std::shared_ptr<std::unordered_set<size_t>> SequenceRange::getSet() {
    return seqIndecies;
}
