//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"

std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    int shift;
    int SAIndex;
    int rangeIndex;
    std::vector<int> indexVector;

    for(SAIndex = 0; SAIndex < SAvector.size(); ++SAIndex) {
        rangeIndex = 0;
        //Dont need to check boundaries since index guaranteed to be less than last item.
        while(SAvector[SAIndex] > seqRangeArray[rangeIndex]){
            ++rangeIndex;
        }
        shift = rangeIndex > 0 ? seqRangeArray[rangeIndex-1]+1 : 0;
        SAvector[SAIndex] -= shift;
        indexVector.push_back(rangeIndex);
    }

    return std::make_shared<std::vector<int>>(indexVector);
}

std::shared_ptr<std::unordered_map<std::string, MatchLocations>> Determine_Matches_Parent(
        std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
        int maximumMatchSize, int numSequences, std::vector<std::shared_ptr<std::string>> &seqStringVector,
        std::vector<int> &seqRangeVector){
    int i;
    int threadMapIndex;
    size_t numSplits = std::thread::hardware_concurrency();
    std::thread threadArray[numSplits];
    size_t splitSize = SAVector.size() / numSplits;
    size_t splitRemainder = SAVector.size() % numSplits; //Determine remainder, and add to first task.
    size_t startIndex = 0; //Guaranteed for LCP[0..numSequences] == 0 due to sentinel higher lexographical order (Only for 1st thread)
    size_t endIndex = splitSize + splitRemainder;
    std::vector<std::shared_ptr<std::unordered_map<std::string, MatchLocations>>> threadMapVector;

    //Initialize Vector
    for (i = 0; i < numSplits; ++i) { //Check if works
        std::shared_ptr<std::unordered_map<std::string, MatchLocations>> matchesMap = std::make_shared<std::unordered_map<std::string, MatchLocations>>();
        threadMapVector.push_back(matchesMap);
    }

    //Create 1st Thread
    threadArray[0] = std::thread(Determine_Matches_Child, std::ref(*threadMapVector[0]), std::ref(LCPVector),
                                 std::ref(SAVector), std::ref(indexVector), minimumMatchSize, maximumMatchSize,
                                 numSequences, std::ref(seqStringVector), std::ref(seqRangeVector), startIndex, endIndex);
//    Determine_Matches_Child(*threadMapVector[0], LCPVector, SAVector, indexVector, minimumMatchSize, maximumMatchSize, numSequences, seqStringVector, seqRangeVector, startIndex, endIndex);  //Single Threaded for Testing


    //Create subsequent threads
    for (i = 1; i < numSplits; ++i) {
        startIndex = endIndex;
        endIndex = startIndex+splitSize;
        threadArray[i] = std::thread(Determine_Matches_Child, std::ref(*threadMapVector[i]), std::ref(LCPVector),
                                     std::ref(SAVector), std::ref(indexVector), minimumMatchSize, maximumMatchSize,
                                     numSequences, std::ref(seqStringVector), std::ref(seqRangeVector), startIndex, endIndex);
//        Determine_Matches_Child(*threadMapVector[i], LCPVector, SAVector, indexVector, minimumMatchSize, maximumMatchSize, numSequences, seqStringVector, seqRangeVector, startIndex, endIndex); //Single Threaded for Testing
    }

    //Wait for all threads to complete
    //Todo Consider adding std::this_thread::yield(); Either in for loop or outside
    for (i=0; i<numSplits; i++){
        threadArray[i].join();
    }

    std::cout << "Merging Match Threads" << std::endl;

    //Merge matchesMap
    std::unordered_map<std::string, MatchLocations> matchesMap;
    std::vector<std::string> matchesMapKeys;

    for (threadMapIndex = 0; threadMapIndex < numSplits; ++threadMapIndex) {
        std::cout << "Merging Thread: " << threadMapIndex << std::endl;
        for(auto &threadMapMatch: *threadMapVector.at(threadMapIndex)) { //For each match in threadMap
            if(matchesMap.count(threadMapMatch.first)){ //If match already exists, add sequence indecies
                matchesMap.at(threadMapMatch.first).mergeMatches(threadMapMatch.second.getMatchVector(), threadMapMatch.second.getUniqueSeqIndices());
            }
            else{ //Create new pair and push to map
                matchesMapKeys.push_back(threadMapMatch.first);
                std::pair<std::string, MatchLocations> newMatchPair = std::make_pair(threadMapMatch.first, threadMapMatch.second);
                matchesMap.insert(newMatchPair);
            }
        }
        threadMapVector.at(threadMapIndex)->clear();
    }

    std::cout <<  "Removing non-qualifying matches" << std::endl;
    //Check if has items from all sequences.
    for(auto & key: matchesMapKeys){
        if(matchesMap.at(key).numSeqIncluded() < numSequences){
            matchesMap.erase(key);
        }
    }

    return std::make_shared<std::unordered_map<std::string, MatchLocations>>(matchesMap);
}




void Determine_Matches_Child(std::unordered_map<std::string, MatchLocations> &matchesMap, std::vector<int> &LCPVector,
                             std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
                             int maximumMatchSize, int numSequences,
                             std::vector<std::shared_ptr<std::string>> &seqStringVector,
                             std::vector<int> &seqRangeVector, size_t startIndex, size_t endIndex) {
    size_t index = startIndex;
    size_t hardEndIndex = SAVector.size();
    std::shared_ptr<std::vector<std::unordered_set<int>>> matchVector;
    std::shared_ptr<std::string> newMatchString;
    std::vector<std::string> matchesMapKeys;

    std::shared_ptr<MatchLocations> newMatchLocation;
    std::vector<int> partitionShiftList;
    std::shared_ptr<std::vector<std::shared_ptr<std::string>>> partitions;
    std::shared_ptr<std::vector<int>> partitionsShiftList = std::make_shared<std::vector<int>>(partitionShiftList);

    while(index < endIndex){
        if (LCPVector.at(index) >= minimumMatchSize){
            while((index  < hardEndIndex) && (LCPVector.at(index) >= minimumMatchSize)){ //Determine end of run (Go Past end if match overlaps)
                //For each match, partition by taking prefix of suffix
                partitions = Determine_Partitions(seqStringVector[indexVector[index]]->substr(SAVector.at(index), LCPVector.at(index)),
                                                  LCPVector.at(index), minimumMatchSize, maximumMatchSize,
                                                  partitionsShiftList);
                //Check all partitions against map
                for (int i = 0; i < partitions->size(); ++i) { //Iterate all partitions.
                    if(matchesMap.count(*partitions->at(i))){ //If partition exists, add match location
                        matchesMap.at(*partitions->at(i)).InsertMatch(SAVector.at(index - 1) + partitionsShiftList->at(i), indexVector.at(index - 1));
                        matchesMap.at(*partitions->at(i)).InsertMatch(SAVector.at(index) + partitionsShiftList->at(i), indexVector.at(index));
                    }
                    else{ //Create new possible match
                        newMatchLocation = std::make_shared<MatchLocations>(numSequences);
                        newMatchLocation->InsertMatch(SAVector.at(index - 1) + partitionsShiftList->at(i), indexVector.at(index - 1));
                        newMatchLocation->InsertMatch(SAVector.at(index) + partitionsShiftList->at(i), indexVector.at(index));
                        std::pair<std::string, MatchLocations> newMatchPair = std::make_pair(*partitions->at(i), *newMatchLocation);
                        matchesMap.insert(newMatchPair);
                        matchesMapKeys.push_back(*partitions->at(i));
                    }
                }
                partitions->clear();
                partitionsShiftList->clear();
                ++index;
            }
        }
        ++index;
    }

}

std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_Partitions(const std::string &key, const int &keyLen, const int &minLength, const int &maxLength,
                     std::shared_ptr<std::vector<int>> &partitionsShiftList) {
    std::vector<std::shared_ptr<std::string>> partitionsStringList;
    int partitionLength = keyLen < maxLength ? keyLen : maxLength;
    for(; partitionLength >= minLength; --partitionLength){
        for (int partitionIndex = 0; (partitionIndex + partitionLength) <= keyLen ; ++partitionIndex){
            partitionsShiftList->push_back(partitionIndex);
            partitionsStringList.push_back(std::make_shared<std::string>(std::string(key.substr(partitionIndex, partitionLength))));
        }
    }

    return std::make_shared<std::vector<std::shared_ptr<std::string>>>(partitionsStringList);
}


std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap, std::string &seqStringCombined,
                            std::vector<std::shared_ptr<std::string>> &seqStringVector,
                            std::vector<int> &seqRangeVector, int &numSequences) {
    int seqIndex;
    int matchSet;
    int matchIndices;
    int overallCoverageSize = 0;
    std::vector<double> similarityMetricVector;
    std::vector<std::unordered_set<int>> matchVector;
    std::vector<std::unordered_set<int>> matchCoverageVector;
    for (seqIndex = 0; seqIndex < numSequences; ++seqIndex) {
        matchCoverageVector.emplace_back(std::unordered_set<int>());
    }

    for(auto &match: matchesMap){
        matchVector = *(match.second.getMatchVector());
        for(matchSet = 0; matchSet < numSequences; ++matchSet){
            for(auto &matchIndex: matchVector.at(matchSet)){
                for(matchIndices = matchIndex; matchIndices < (matchIndex + match.first.length()); ++matchIndices){
                    matchCoverageVector[matchSet].insert(matchIndices);
                }
            }
        }
    }

    for(seqIndex = 0; seqIndex < numSequences; ++seqIndex){
        overallCoverageSize += static_cast<int>(matchCoverageVector[seqIndex].size());
        similarityMetricVector.push_back((double)matchCoverageVector[seqIndex].size()/(double)(seqStringVector[seqIndex]->length())); //Subtracting number of sentinels before current sequence
    }

    similarityMetricVector.push_back((double) overallCoverageSize / (double) (seqStringCombined.length()-numSequences)); //Subtracting total number of sentinels


    return std::make_shared<std::vector<double>>(similarityMetricVector);
}
