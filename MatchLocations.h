//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
#define SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H

#include <string>
#include <vector>
#include <unordered_set>

class MatchLocations {
private:
    std::vector<std::shared_ptr<std::unordered_set<size_t>>> matchVector;
public:
    MatchLocations();
    void insertSet(std::shared_ptr<std::unordered_set<size_t>> indexSet);
    void insertIndex(size_t setIndex, size_t index);
    std::shared_ptr<std::vector<std::shared_ptr<std::unordered_set<size_t>>>> getMatchVector();
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
