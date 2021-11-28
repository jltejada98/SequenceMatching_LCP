//
// Created by Jose Luis Tejada on 1/8/21.
//

#ifndef SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
#define SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H

#include <unordered_set>
#include <vector>

class MatchLocations {
private:
    std::shared_ptr<std::vector<std::unordered_set<int>>> matchVector;
public:
    MatchLocations();
    void insertMatchVector(std::shared_ptr<std::vector<std::unordered_set<int>>> vector);
    void addMatches(const std::shared_ptr<std::vector<std::unordered_set<int>>>& vector);
    std::shared_ptr<std::vector<std::unordered_set<int>>> getMatchVector();
};


#endif //SEQUENCEMATCHING_LCP_MATCHLOCATIONS_H
