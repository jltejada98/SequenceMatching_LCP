//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"

std::shared_ptr<std::vector<size_t>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<size_t> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    size_t SAIndex;
    int SAShift;
    size_t rangeIndex;
    std::vector<size_t> indexVector;

    for(SAIndex = 0; SAIndex < SAvector.size(); ++SAIndex) {
        rangeIndex = 0;
        //Dont need to check boundaries since index guaranteed to be less than last item.
        while(SAvector[SAIndex] > seqRangeArray[rangeIndex]){  //Equality added for sentinels
            ++rangeIndex;
        }
        indexVector.push_back(rangeIndex);
        //Shift suffix array back to original sequence positions (Taking into account sentinels before it)
        SAShift = rangeIndex > 0 ? seqRangeArray[rangeIndex-1] + rangeIndex : 0;
        SAvector.at(SAIndex) -= SAShift;
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

            std::cout << "Completed Run" << std::endl;
            // Todo Consider using threads to determine submatching

            //SUBMATCHING ALGORITHM
            //Perform submatching between SIRs by comparing the currentSmallest LCP value with the smallest LCP value of each of the SIR.
            std::map<size_t, PossibleMatches>::iterator currMapItem;
            std::map<size_t, PossibleMatches>::iterator currMapEnd;
            size_t runMapVectorSize = runMapVector.size();
            size_t checkMap;
            std::map<size_t, PossibleMatches>::iterator checkMapItem;
            std::map<size_t, PossibleMatches>::iterator checkMapEnd;
            for (size_t currentMap = 0; currentMap < runMapVectorSize; ++currentMap) {
                //Perform submatching within a strictly increasing run.
                currMapItem = runMapVector.at(currentMap)->begin();
                currMapEnd = runMapVector.at(currentMap)->end();
                while(currMapItem != currMapEnd){ //Todo Check if iteration correct
                    auto nextMapItem = std::next(currMapItem);
                    while(nextMapItem != currMapEnd){
                        //Add Matches from the current map Item to the next one:
                        //Map is sorted and it is guaranteed that in a strictly increasing run, it includes the match.
                        currMapItem->second.AddSubMatch(nextMapItem->second.getSAIndexMap());
                        ++nextMapItem;
                    }
                    ++currMapItem;
                }
                //Perform submatching on other strictly increasing runs by comparing:
                //minLCP(CurrentRun) < minLCP(i)  -> Then add match indices from minLCP[CurrentRun] to run I
                checkMap = 0;
                currMapItem = runMapVector.at(currentMap)->begin();
                currMapEnd = runMapVector.at(currentMap)->end();
                while((checkMap < runMapVectorSize) && (currMapItem->first < runMapVector.at(checkMap)->begin()->first)){ //Check min(LCP[currMap]) <= min(LCP[i])
                    if(checkMap != currentMap){
                        //Add submatches for currMapItem, against  all of checkMap
                        // Todo Determine if need to increment  currMapItem
                        // Todo Determine if only need to compare against
                        checkMapItem = runMapVector.at(checkMap)->begin();
                        checkMapEnd = runMapVector.at(checkMap)->end();
                        while(checkMapItem != checkMapEnd){
                            currMapItem->second.AddSubMatch(checkMapItem->second.getSAIndexMap());
                            ++checkMapItem;
                        }
                    }
                    ++checkMap;
                }
            }

            std::cout << "Completed Submatching" << std::endl;

            //Determine if each of the possible match lengths has indices from all sequences
            for(auto &map : runMapVector){
                for(auto &match: *map){
                    if(match.second.UniqueSetIndices() >= numSequences){ //Includes Indices from all sequences
                        matchLength = match.first;
                        matchVector = match.second.TransferMatches();
                        matchIndexForString = *(matchVector->at(0).begin());
                        newMatchString = std::make_shared<std::string>(seqStringVector.at(0)->substr(matchIndexForString, matchLength));
                        if(matchesMap.count(*newMatchString)){ //MatchString already exists, merge matches
                            matchesMap.at(*newMatchString).addMatches(matchVector);
                        }
                        else{
                            newMatchLocation = std::make_shared<MatchLocations>();
                            newMatchLocation->insertMatchVector(matchVector);
                            newMatchPair = std::make_pair(*newMatchString, *newMatchLocation);
                            matchesMap.insert(newMatchPair);
                        }
                    }
                }
                map->clear();
            }

            //Clear runMapVector Todo determine if need to clear individual
            runMapVector.clear();
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
