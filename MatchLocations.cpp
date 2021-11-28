//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

#include <utility>

MatchLocations::MatchLocations() {
    matchVector = nullptr;
}

void MatchLocations::insertMatchVector(std::shared_ptr<std::vector<std::unordered_set<int>>> vector) {
    matchVector = std::move(vector);
}

void MatchLocations::addMatches(const std::shared_ptr<std::vector<std::unordered_set<int>>>& vector) {
    int vectorSize = vector->size();
    for (int seqInd = 0; seqInd < vectorSize; ++seqInd) { //Insert all matches from set into corresponding vector.
        matchVector->at(seqInd).insert(vector->at(seqInd).begin(), vector->at(seqInd).end()); //Insertion only keeps unique indecies
    }

}

std::shared_ptr<std::vector<std::unordered_set<int>>> MatchLocations::getMatchVector() {
    return matchVector;
}

