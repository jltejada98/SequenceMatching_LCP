//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"

std::shared_ptr<std::vector<size_t>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<size_t> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    size_t SAIndex;
//    size_t SAShift;
    size_t rangeIndex;
    std::vector<size_t> indexVector;

    for(SAIndex = 0; SAIndex < SAvector.size(); ++SAIndex) {
        rangeIndex = 0;
        //Dont need to check boundaries since index guaranteed to be less than last item.
        while(SAvector[SAIndex] > seqRangeArray[rangeIndex]){  //Equality added for sentinels
            ++rangeIndex;
        }
        indexVector.push_back(rangeIndex);
        //Shift suffix array back to original sequence positions
//        SAShift = rangeIndex > 0 ? seqRangeArray[rangeIndex] : 0;
//        SAvector[SAIndex] -= SAShift;
    }

    return std::make_shared<std::vector<size_t>>(indexVector);
}


std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches(std::vector<int> &LCPVector, std::vector<int> &SAVector,
                  std::vector<size_t> &indexVector, std::vector<std::shared_ptr<std::string>> &seqStringVector, size_t &minimumMatchSize,
                  size_t &numSequences){
    size_t LCPVectorSize = LCPVector.size();
    size_t index = numSequences; //Guaranteed for LCP[0..numSequences] == 0 due to sentinel higher lexographical order.
    size_t matchLength;
    size_t matchIndexForString;
    std::unordered_map<std::string, MatchLocations> matchesMap;
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> matchVector;
    std::shared_ptr<std::string> newMatchString;
    std::shared_ptr<MatchLocations> newMatchLocation;
    std::pair<std::string,MatchLocations> newMatchPair;


    //Todo Consider using ordered map for adding substring indices across all runMaps in vector.
    std::vector<std::shared_ptr<std::map<size_t, PossibleMatches>>> runMapVector;
    std::shared_ptr<std::map<size_t, PossibleMatches>> runMap;
    std::pair<size_t, PossibleMatches> runPair;
    std::shared_ptr<PossibleMatches> runLocation;
    std::shared_ptr<int> runLength;

    //New Algorithhm
    // Determine Run (Qualified by all LCP values >= minLength)
        // Within Run, determine strictly increasing run
        // Add matches to runMap and check for existance based on LCP value (key in map)
    //After determining all strictly increasing runs, add submatching by adding indecies to all matches with ength >= current length
    //Continue adding matches based on runMap.



    while(index < LCPVectorSize){
        if (LCPVector.at(index) >= minimumMatchSize){
            //Todo Consider using threads to find next run?
            //Determine runs and possible matches
            while((index  < LCPVectorSize) && (LCPVector.at(index) >= minimumMatchSize)){ //Determine end of run
                //Create new runMap
                runMap = std::make_shared<std::map<size_t, PossibleMatches>>();
                runLocation = std::make_shared<PossibleMatches>();
                //SA comparisons start at n-1. (LCP[index] -> comparing SA[index-1] and SA[index])
                runLocation->InsertMatch(SAVector.at(index - 1), indexVector.at(index - 1));
                runLocation->InsertMatch(SAVector.at(index), indexVector.at(index));
                runPair = std::make_pair(LCPVector.at(index), *runLocation);
                runMap->insert(runPair);
                //Due to SA[index-1] being included in previous strictly increasing run, we must include it in the previous runMap to perform submatching
                if(!runMapVector.empty()){
                    runMapVector.back()->insert(runPair);
                }
                ++index;

                while((index < LCPVectorSize) && (LCPVector.at(index) >= LCPVector.at(index - 1))){ //Determine end of strictly increasing run
                    if(runMap->count(LCPVector.at(index))){ //If matchLength exists
                        runMap->at(LCPVector.at(index)).InsertMatch(SAVector.at(index), indexVector.at(index));
                    }
                    else{ //Create new possible match
                        runLocation = std::make_shared<PossibleMatches>();
                        //LCP[index] -> comparing SA[index-1] and SA[index]
                        runLocation->InsertMatch(SAVector.at(index - 1), indexVector.at(index - 1));
                        runLocation->InsertMatch(SAVector.at(index), indexVector.at(index));
                        runPair = std::make_pair(LCPVector.at(index), *runLocation);
                        runMap->insert(runPair);
                    }
                    ++index;
                }
                runMapVector.push_back(runMap);
            }

            // Todo Consider using threads to determine submatching
            //Determine submatching amongst strictly increasing runs.
            for(auto &map: runMapVector){ //In each strictly increasing map, we can perform submatching for all of the strings
                auto currMapItem = map->begin();
                auto currMapEnd = map->end();
                while(currMapItem != currMapEnd){
                      auto nextMapItem = std::next(currMapItem);
                      while(nextMapItem != currMapEnd){
                          //Add Matches from the current map Item to the next one:
                          //Map is sorted and it is guaranteed that in a strictly increasing run, it includes the match.
                          currMapItem->second.AddSubMatch(nextMapItem->second.getSA(),nextMapItem->second.getSeq());
//                          nextMapItem->second.AddSubMatch(currMapItem->second.getSA(), currMapItem->second.getSeq());
                          ++nextMapItem;
                      }
                      ++currMapItem;
                }
            }

            std::cout << "Complete" << std::endl;

            //Determine if each of the possible match lengths has indices from all sequences
//            for (auto &run: runMap) {
//                if(run.second.UniqueSetIndices() >= numSequences){ //Indices from all sequences
//                    //Insert match onto matchesMap
//                    matchLength = run.first;
//                    matchVector = run.second.TransferMatches();
//                    matchIndexForString = *(matchVector->at(0).begin());
//                    newMatchString = std::make_shared<std::string>(seqStringVector.at(0)->substr(matchIndexForString, matchLength));
//                    //No need to check if match exists since lexicographical ordering of LCP array guarantees this
//                    newMatchLocation = std::make_shared<MatchLocations>();
//                    newMatchLocation->insertMatchVector(matchVector);
//                    newMatchPair = std::make_pair(*newMatchString, *newMatchLocation);
//                    matchesMap.insert(newMatchPair);
//                }
//            }
//            //Perform setup for next run
//            runMap.clear();
        }
        else{
            ++index;
        }
    }


    return std::make_shared<std::unordered_map<std::string, MatchLocations>>(matchesMap);
}


std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap, std::string &seqStringCombined,
                            std::vector<std::shared_ptr<std::string>> &seqStringArray,
                            size_t &numSequences){
    size_t index;
    size_t seqIndex;
    size_t matchLength;
    size_t overallCoverageSize = 0;
    std::vector<double> similarityMetricVector;
    std::vector<std::unordered_set<size_t>> seqVectorSet;
    std::vector<std::unordered_set<size_t>> matchIndicesVector;
    for (seqIndex = 0; seqIndex < numSequences; ++seqIndex) {
        seqVectorSet.emplace_back(std::unordered_set<size_t>());
    }
    for (auto &match: matchesMap){
        matchIndicesVector = *match.second.getMatchVector();
        matchLength = match.first.length();
        for (seqIndex = 0; seqIndex < numSequences ; ++seqIndex) {
            //For each index in each sequence set.
            for (auto &seqSetIndex : matchIndicesVector[seqIndex]){
                for (index = seqSetIndex; index < (seqSetIndex+matchLength); ++index){
                    seqVectorSet[seqIndex].insert(index);
                }
            }
        }
    }
    for(seqIndex = 0; seqIndex < numSequences; ++seqIndex){
        overallCoverageSize += seqVectorSet[seqIndex].size();
        similarityMetricVector.emplace_back((double)seqVectorSet.size()/(double)seqStringArray[seqIndex]->length());
    }
    similarityMetricVector.emplace_back((double) overallCoverageSize /(double)(seqStringCombined.length() - numSequences));

    return std::make_shared<std::vector<double>>(similarityMetricVector);
}
