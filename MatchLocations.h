//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
#define SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H

#include <unordered_set>
#include <vector>

class MatchLocations {
private:
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> matchVector;
public:
    MatchLocations();
    void insertMatchVector(std::shared_ptr<std::vector<std::unordered_set<size_t>>> vector);
    void addMatches(std::shared_ptr<std::vector<std::unordered_set<size_t>>> vector);
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> getMatchVector();
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
