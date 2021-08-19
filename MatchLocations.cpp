//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

MatchLocations::MatchLocations() {
    matchVector = nullptr;
}

void MatchLocations::insertMatchVector(std::shared_ptr<std::vector<std::unordered_set<size_t>>> vector) {
    matchVector = vector;
}

void MatchLocations::addMatches(std::shared_ptr<std::vector<std::unordered_set<size_t>>> vector) {
    size_t vectorSize = vector->size();
    for (int seqInd = 0; seqInd < vectorSize; ++seqInd) { //Insert all matches from set into corresponding vector.
        matchVector->at(seqInd).insert(vector->at(seqInd).begin(), vector->at(seqInd).end());
    }

}

std::shared_ptr<std::vector<std::unordered_set<size_t>>> MatchLocations::getMatchVector() {
    return matchVector;
}

