//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"

std::shared_ptr<std::vector<size_t>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<size_t> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    size_t SAIndex;
    size_t rangeIndex;
    std::vector<size_t> indexVector;

    for(SAIndex = 0; SAIndex < SAvector.size(); ++SAIndex) {
        rangeIndex = 0;
        //Dont need to check boundaries since index guaranteed to be less than last item.
        while(SAvector[SAIndex] > seqRangeArray[rangeIndex]){
            ++rangeIndex;
        }
        indexVector.push_back(rangeIndex);
    }

    return std::make_shared<std::vector<size_t>>(indexVector);
}


std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches(std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<size_t> &indexVector,
                  std::vector<std::shared_ptr<std::string>> &seqStringVector, size_t &minimumMatchSize,
                  size_t &numSequences, std::string &seqStringCombined) {
    size_t LCPVectorSize = LCPVector.size();
    size_t index = numSequences; //Guaranteed for LCP[0..numSequences] == 0 due to sentinel higher lexographical order.
    std::unordered_map<std::string, MatchLocations> matchesMap;
    std::shared_ptr<std::vector<std::unordered_set<size_t>>> matchVector;
    std::shared_ptr<std::string> newMatchString;
    std::shared_ptr<MatchLocations> newMatchLocation;
    std::pair<std::string,MatchLocations> newMatchPair;
    //Possible matches map
    std::shared_ptr<std::map<std::string, PossibleMatches>> runMap;
    std::pair<std::string, PossibleMatches> runPair;
    std::shared_ptr<PossibleMatches> runLocation;
    std::vector<int> partitionShiftList;
    std::shared_ptr<std::vector<std::shared_ptr<std::string>>> partitions;
    std::shared_ptr<std::vector<int>> partitionsShiftList = std::make_shared<std::vector<int>>(partitionShiftList);


    while(index < LCPVectorSize){
        if (LCPVector.at(index) >= minimumMatchSize){
            //Todo Consider using threads to find next run?
            //Determine runs and possible matches
            runMap = std::make_shared<std::map<std::string, PossibleMatches>>();
            while((index  < LCPVectorSize) && (LCPVector.at(index) >= minimumMatchSize)){ //Determine end of run
                //For each match, partition by taking prefix of suffix
                partitions = Determine_Partitions(seqStringCombined.substr(SAVector.at(index), LCPVector.at(index)), LCPVector.at(index), minimumMatchSize, partitionsShiftList);
                //Check all partitions against map
                for (size_t i = 0; i < partitions->size(); ++i) { //Iterate all partitions.
                    if(runMap->count(*partitions->at(i))){ //If partition exists, add match location
                        runMap->at(*partitions->at(i)).InsertMatch(SAVector.at(index) + partitionsShiftList->at(i), indexVector.at(index));
                    }
                    else{ //Create new possible match
                        runLocation = std::make_shared<PossibleMatches>();
                        runLocation->InsertMatch(SAVector.at(index - 1) + partitionsShiftList->at(i), indexVector.at(index - 1));
                        runLocation->InsertMatch(SAVector.at(index) + partitionsShiftList->at(i), indexVector.at(index));
                        runPair = std::make_pair(*partitions->at(i), *runLocation);
                        runMap->insert(runPair);
                    }
                }
                partitions->clear();
                partitionsShiftList->clear();
                ++index;
            }

            std::cout << "RunMap Complete" << std::endl;

            //Determine if each of the possible matches has indices from all sequences
            for(auto &match: *runMap){
                if(match.second.UniqueSetIndices() >= numSequences){
                    matchVector = match.second.TransferMatches();
                    if(matchesMap.count(match.first)){ //MatchString already exists, merge matches
                        matchesMap.at(match.first).addMatches(matchVector);
                    }
                    else{
                        newMatchString = std::make_shared<std::string>(match.first);
                        newMatchLocation = std::make_shared<MatchLocations>();
                        newMatchLocation->insertMatchVector(matchVector);
                        newMatchPair = std::make_pair(*newMatchString, *newMatchLocation);
                        matchesMap.insert(newMatchPair);
                    }
                }
            }
            runMap->clear();
        }
        else{
            ++index;
        }
    }


    return std::make_shared<std::unordered_map<std::string, MatchLocations>>(matchesMap);
}

std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_Partitions(const std::string &key, const int &keyLen, const size_t &minLength,
                     std::shared_ptr<std::vector<int>> &partitionsShiftList){
    std::vector<std::shared_ptr<std::string>> partitionsStringList;

    for(int p = keyLen; p >= minLength; --p){
        for (int i = 0; (i+p) <= keyLen ; ++i){
            partitionsShiftList->push_back(i);
            partitionsStringList.push_back(std::make_shared<std::string>(std::string(key.substr(i, p))));
        }
    }

    return std::make_shared<std::vector<std::shared_ptr<std::string>>>(partitionsStringList);
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
