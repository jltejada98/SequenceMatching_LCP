//
// Created by Jose Luis Tejada on 1/8/21.
//

#include "MatchLocations.h"

MatchLocations::MatchLocations() {

}

//Assumes insertion in sequence order
void MatchLocations::insertSet(std::shared_ptr<std::unordered_set<size_t>> indexSet) {
    matchVector.push_back(indexSet);
}

std::shared_ptr<std::vector<std::shared_ptr<std::unordered_set<size_t>>>> MatchLocations::getMatchVector() {
    return std::make_shared<std::vector<std::shared_ptr<std::unordered_set<size_t>>>>(matchVector);
}

void MatchLocations::insertIndex(size_t setIndex, size_t index) {
    matchVector.at(setIndex)->insert(index);
}

