//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"

std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    int SAIndex;
    int rangeIndex;
    std::vector<int> indexVector;

    for(SAIndex = 0; SAIndex < SAvector.size(); ++SAIndex) {
        rangeIndex = 0;
        //Dont need to check boundaries since index guaranteed to be less than last item.
        while(SAvector[SAIndex] > seqRangeArray[rangeIndex]){
            ++rangeIndex;
        }
        indexVector.push_back(rangeIndex);
    }

    return std::make_shared<std::vector<int>>(indexVector);
}


std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches(std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector,
                  std::vector<std::shared_ptr<std::string>> &seqStringVector, int &minimumMatchSize,
                  int &numSequences, std::string &seqStringCombined) {
    int LCPVectorSize = static_cast<int>(LCPVector.size());
    int index = numSequences; //Guaranteed for LCP[0..numSequences] == 0 due to sentinel higher lexographical order.
    std::unordered_map<std::string, MatchLocations> matchesMap;
    std::shared_ptr<std::vector<std::unordered_set<int>>> matchVector;
    std::shared_ptr<std::string> newMatchString;
    std::shared_ptr<MatchLocations> newMatchLocation;
    std::pair<std::string,MatchLocations> newMatchPair;
    //Todo Merge possible matches and matches map classes


    //Possible matches map
    std::shared_ptr<std::map<std::string, PossibleMatches>> possibleMatchesMap;
    std::pair<std::string, PossibleMatches> possibleMatchPair;
    std::shared_ptr<PossibleMatches> possibleMatchLocation;
    std::vector<int> partitionShiftList;
    std::shared_ptr<std::vector<std::shared_ptr<std::string>>> partitions;
    std::shared_ptr<std::vector<int>> partitionsShiftList = std::make_shared<std::vector<int>>(partitionShiftList);
    possibleMatchesMap = std::make_shared<std::map<std::string, PossibleMatches>>(); //Determine runs and possible matches


    while(index < LCPVectorSize){
        if (LCPVector.at(index) >= minimumMatchSize){
            //Todo Consider using threads to find next run?
            while((index  < LCPVectorSize) && (LCPVector.at(index) >= minimumMatchSize)){ //Determine end of run
                //For each match, partition by taking prefix of suffix
                partitions = Determine_Partitions(seqStringCombined.substr(SAVector.at(index), LCPVector.at(index)), LCPVector.at(index), minimumMatchSize, partitionsShiftList);
                //Check all partitions against map
                for (int i = 0; i < partitions->size(); ++i) { //Iterate all partitions.
                    if(possibleMatchesMap->count(*partitions->at(i))){ //If partition exists, add match location
                        possibleMatchesMap->at(*partitions->at(i)).InsertMatch(SAVector.at(index - 1) + partitionsShiftList->at(i), indexVector.at(index - 1));
                        possibleMatchesMap->at(*partitions->at(i)).InsertMatch(SAVector.at(index) + partitionsShiftList->at(i), indexVector.at(index));
                    }
                    else{ //Create new possible match
                        possibleMatchLocation = std::make_shared<PossibleMatches>();
                        possibleMatchLocation->InsertMatch(SAVector.at(index - 1) + partitionsShiftList->at(i), indexVector.at(index - 1));
                        possibleMatchLocation->InsertMatch(SAVector.at(index) + partitionsShiftList->at(i), indexVector.at(index));
                        possibleMatchPair = std::make_pair(*partitions->at(i), *possibleMatchLocation);
                        possibleMatchesMap->insert(possibleMatchPair);
                    }
                }
                partitions->clear();
                partitionsShiftList->clear();
                ++index;
            }
        }
        ++index;
    }

    //Determine if each of the possible matches has indices from all sequences
    for(auto &match: *possibleMatchesMap){
        if(match.second.UniqueSetIndices() >= numSequences){ //No need to check for existance, since it is guaranteed to exist.
            matchVector = match.second.TransferMatches();
            newMatchString = std::make_shared<std::string>(match.first);
            newMatchLocation = std::make_shared<MatchLocations>();
            newMatchLocation->insertMatchVector(matchVector);
            newMatchPair = std::make_pair(*newMatchString, *newMatchLocation);
            matchesMap.insert(newMatchPair);
        }
    }

    possibleMatchesMap->clear();


    return std::make_shared<std::unordered_map<std::string, MatchLocations>>(matchesMap);
}

std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_Partitions(const std::string &key, const int &keyLen, const int &minLength,
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
                            int &numSequences){
    int index;
    int seqIndex;
    int matchLength;
    int overallCoverageSize = 0;
    std::vector<double> similarityMetricVector;
    std::vector<std::unordered_set<int>> seqVectorSet;
    std::vector<std::unordered_set<int>> matchIndicesVector;
    for (seqIndex = 0; seqIndex < numSequences; ++seqIndex) {
        seqVectorSet.emplace_back(std::unordered_set<int>());
    }
    for (auto &match: matchesMap){
        matchIndicesVector = *match.second.getMatchVector();
        matchLength = static_cast<int>(match.first.length());
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
        overallCoverageSize += static_cast<int>(seqVectorSet[seqIndex].size());
        similarityMetricVector.emplace_back((double)seqVectorSet.size()/(double)seqStringArray[seqIndex]->length());
    }
    similarityMetricVector.emplace_back((double) overallCoverageSize /(double)(seqStringCombined.length() - numSequences));

    return std::make_shared<std::vector<double>>(similarityMetricVector);
}
