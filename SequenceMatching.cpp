//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"

std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    int shift;
    int rangeIndex;
    size_t SAIndex;
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

//Testing new method of determining matches: Doesn't use memory for non-valid matches, removes them before determining locations
//1) For all new child threads
//   -Determine if a given possible match satisfies min,max
//   -Determine all partitions, and store their sequence index (i.e: Which sequence they belong to)
//   -Store the corresponding LCP vector.
//   -For all partitions, determine if already exists, if so, add their sequence index.
//2) Merge all matches found by threads, specifically their sequence index.
//3) Determine if a given match includes sequences from all indecies.
//    - If not remove the match
//4) Create new child threads to iterate through valid remaining matches using the LCP indecies previously stored to expand and
//   generate all matches.


std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches_Parent(std::vector<int> &LCPVector, std::vector<int> &SAVector, std::vector<int> &indexVector,
                         int minimumMatchSize, int maximumMatchSize, int numSequences,
                         std::vector<std::shared_ptr<std::string>> &seqStringVector) {
    size_t i;
    size_t threadMapIndex;
    size_t numSplits = std::thread::hardware_concurrency();
    std::thread threadArray[numSplits];
    size_t splitSize = SAVector.size() / numSplits;
    size_t splitRemainder = SAVector.size() % numSplits; //Determine remainder, and add to first task.
    size_t startIndex = 0; //Guaranteed for LCP[0..numSequences] == 0 due to sentinel higher lexographical order (Only for 1st thread)
    size_t endIndex = splitSize + splitRemainder;
    std::vector<std::shared_ptr<std::unordered_map<std::string, MatchLocations>>> threadMapVector;


    //Todo change print messages
    //////////////////////////////////////////////////////////////////////////
    //1. Remove invalid matches, by determining number of sequences present.//
    //////////////////////////////////////////////////////////////////////////

    //Initialize Vector
    for (i = 0; i < numSplits; ++i) { //Check if works
        std::shared_ptr<std::unordered_map<std::string, MatchLocations>> matchesMap = std::make_shared<std::unordered_map<std::string, MatchLocations>>();
        threadMapVector.push_back(matchesMap);
    }

    //Create 1st Thread
//    threadArray[0] = std::thread(Determine_Valid_Matches_Child, std::ref(*threadMapVector[0]), std::ref(LCPVector),
//                                 std::ref(SAVector), std::ref(indexVector), minimumMatchSize, maximumMatchSize,
//                                 std::ref(seqStringVector), std::ref(seqRangeVector), startIndex, endIndex);
    Determine_Valid_Matches_Child(*threadMapVector[0], LCPVector, SAVector, indexVector, minimumMatchSize,
                                  maximumMatchSize, seqStringVector, startIndex,
                                  endIndex);  //Single Threaded for Testing

    //Create subsequent threads
    for (i = 1; i < numSplits; ++i) {
        startIndex = endIndex;
        endIndex = startIndex+splitSize;
//        threadArray[i] = std::thread(Determine_Valid_Matches_Child, std::ref(*threadMapVector[i]), std::ref(LCPVector),
//                                     std::ref(SAVector), std::ref(indexVector), minimumMatchSize, maximumMatchSize,
//                                     startIndex,endIndex);
                Determine_Valid_Matches_Child(*threadMapVector[i], LCPVector, SAVector, indexVector, minimumMatchSize,
                                              maximumMatchSize, seqStringVector, startIndex,
                                              endIndex); //Single Threaded for Testing
    }

    //Wait for all threads to complete
    //Todo Consider adding std::this_thread::yield(); Either in for loop or outside
//    for (i=0; i<numSplits; i++){
//        threadArray[i].join();
//    }

    std::cout << "Finished Determining Valid Matches." << std::endl;

    //Merge all thread validity maps into parentMatchValidity
    std::vector<std::string> keysToCheckValidity;

    //Using threadMapVector[0] as the parent.
    //Add keys from first map to vector
    for(auto const &item: *threadMapVector[0]){
        keysToCheckValidity.push_back(item.first);
    }

    for (threadMapIndex = 1; threadMapIndex < numSplits; ++threadMapIndex) {
        std::cout << "Merging Validity Thread: " << threadMapIndex << std::endl;
        for(auto &threadMapMatch: *threadMapVector[threadMapIndex]) { //For each match in threadMap
            if(threadMapVector[0]->count(threadMapMatch.first)){ //If match already exists, add sequence indecies
                threadMapVector[0]->at(threadMapMatch.first).MergeUniqueSeq_LCPIndex(
                        threadMapMatch.second.GetUniqueSeqIndices(), threadMapMatch.second.GetLCPIndecies());
            }
            else{ //Create new pair and push to map, remove from source map.
                keysToCheckValidity.push_back(threadMapMatch.first); //Creating match
                std::pair<std::string, MatchLocations> newMatchPair = std::make_pair(threadMapMatch.first, threadMapMatch.second);
                threadMapVector[0]->insert(newMatchPair);
            }
        }
        threadMapVector[threadMapIndex]->clear();
    }

    std::cout <<  "Removing invalid matches..." << std::endl;
    //Check if key has items from all sequences.
    std::unordered_set<std::string> validMatches; //Contains valid matches to be used in splitting matches
    std::unordered_set<int> validLCPIndeciesSet; //Contains indecies to be used in determining locations
    for(auto & key: keysToCheckValidity){
        if(threadMapVector[0]->at(key).NumSeqIncluded() == numSequences){
            validMatches.insert(key);
            validLCPIndeciesSet.insert(threadMapVector[0]->at(key).GetLCPIndecies()->begin(),
                                    threadMapVector[0]->at(key).GetLCPIndecies()->end());
        }
    }
    //Cleanup tasks
    keysToCheckValidity.clear();
    threadMapVector[0]->clear(); //ThreadMapVector Completely Cleared, ready for determining valid locations.

    //Add indecies from set to vector.
    std::vector<int> validLCPIndecies(validLCPIndeciesSet.begin(),validLCPIndeciesSet.end());
    validLCPIndeciesSet.clear();

    //////////////////////////////////////////////////////////
    //2. Use threads to determine locations of valid matches//
    //////////////////////////////////////////////////////////
    //Split length of validLCPIndecies amongst each thread. Each thread then proceeds to
    //Check all of the matches in its share and add the match
    splitSize = validLCPIndecies.size() / numSplits;
    splitRemainder = validLCPIndecies.size() % numSplits; //Determine remainder, and add to first task.
    startIndex = 0;
    endIndex = splitSize + splitRemainder;
    //Todo verify that all LCP indecies are being visited

    //Create 1st Thread
//    threadArray[0] = std::thread(Determine_Match_Locations_Child, std::ref(*threadMapVector[0]), std::ref(LCPVector),
//                                 std::ref(SAVector), std::ref(indexVector), std::ref(validMatches), std::ref(validLCPIndecies),
//                                 minimumMatchSize, maximumMatchSize, std::ref(seqStringVector), std::ref(seqRangeVector), startIndex, endIndex);

    //Single Thread Version
//    Determine_Match_Locations_Child(*threadMapVector[0], LCPVector, SAVector, indexVector, validMatches, validLCPIndecies, minimumMatchSize, maximumMatchSize, numSequences, seqStringVector, startIndex, endIndex);

    //All In One
    Determine_Match_Locations_Child(*threadMapVector[0], LCPVector, SAVector, indexVector, validMatches, validLCPIndecies, minimumMatchSize, maximumMatchSize, numSequences, seqStringVector, startIndex, validLCPIndecies.size());


    //Create subsequent threads
//    for (i = 1; i < numSplits; ++i) {
//        startIndex = endIndex;
//        endIndex = startIndex+splitSize;
//        //Single Thread Version
//        Determine_Match_Locations_Child(*threadMapVector[i], LCPVector, SAVector, indexVector, validMatches,
//                                        validLCPIndecies, minimumMatchSize, maximumMatchSize, numSequences,
//                                        seqStringVector, startIndex, endIndex);//Single Threaded for Testing
//    }


    //Wait for all threads to complete
    //Todo Consider adding std::this_thread::yield(); Either in for loop or outside
//    for (i=0; i<numSplits; i++){
//        threadArray[i].join();
//    }

    //Cleanup Taks
    validMatches.clear();
    validLCPIndecies.clear();


    std::cout << "Finished Determining Match Locations." << std::endl;

    //Merge all thread validity maps into parentMatchValidity
    //Using threadMapVector[0] as the parent.
    for (threadMapIndex = 1; threadMapIndex < numSplits; ++threadMapIndex) {
        std::cout << "Merging Location Thread: " << threadMapIndex << std::endl;
        for(auto &threadMapMatch: *threadMapVector[threadMapIndex]) { //For each match in threadMap
            if(threadMapVector[0]->count(threadMapMatch.first)){ //If match already exists, add sequence indecies
                threadMapVector[0]->at(threadMapMatch.first).MergeMatches(threadMapMatch.second.GetMatchVector());
            }
            else{ //Create new pair and push to map, remove from source map.
                std::pair<std::string, MatchLocations> newMatchPair = std::make_pair(threadMapMatch.first, threadMapMatch.second);
                threadMapVector[0]->insert(newMatchPair);
            }
        }
        threadMapVector[threadMapIndex]->clear();
    }


    return std::make_shared<std::unordered_map<std::string, MatchLocations>>(*threadMapVector[0]);
}

void Determine_Match_Locations_Child(std::unordered_map<std::string, MatchLocations> &matchesMap,
                                     std::vector<int> &LCPVector, std::vector<int> &SAVector,
                                     std::vector<int> &indexVector, std::unordered_set<std::string> &validMatches,
                                     std::vector<int> &validLCPIndecies, int minimumMatchSize, int maximumMatchSize,
                                     int numSequences, std::vector<std::shared_ptr<std::string>> &seqStringVector,
                                     size_t startIndex, size_t endIndex) {

    std::shared_ptr<MatchLocations> newMatchLocation;
    std::pair<std::string, MatchLocations> newMatchPair;
    std::shared_ptr<std::vector<std::shared_ptr<std::string>>> partitions;
    std::shared_ptr<std::vector<int>> partitionsShiftList = std::make_shared<std::vector<int>>();
    size_t actualLCPIndex;

    //Iterate across LCP Indecies
    for(size_t validIndexToCheck = startIndex; validIndexToCheck < endIndex; ++validIndexToCheck) {
        //LCP Index that is being queried
        actualLCPIndex = validLCPIndecies[validIndexToCheck];

        //When creating partitions, check if the partition created needs to be created.
        partitions = Determine_LocationPartitions(seqStringVector[indexVector[actualLCPIndex]]->substr(SAVector[actualLCPIndex], LCPVector[actualLCPIndex]),
                                                               LCPVector[actualLCPIndex], minimumMatchSize, maximumMatchSize, partitionsShiftList);

        for(size_t partitionIndex = 0; partitionIndex < partitions->size(); ++partitionIndex){
            if(validMatches.count(*partitions->at(partitionIndex))){ //Determine if any of the partitions are valid
                if(matchesMap.count(*partitions->at(partitionIndex))){
                    matchesMap[*partitions->at(partitionIndex)].InsertLocation(
                            SAVector[actualLCPIndex - 1] + partitionsShiftList->at(partitionIndex),
                            indexVector[actualLCPIndex - 1]);
                    matchesMap[*partitions->at(partitionIndex)].InsertLocation(
                            SAVector[actualLCPIndex] + partitionsShiftList->at(partitionIndex),
                            indexVector[actualLCPIndex]);
                }
                else{
                    newMatchLocation = std::make_shared<MatchLocations>(numSequences);
                    newMatchLocation->InsertLocation(
                            SAVector[actualLCPIndex - 1] + partitionsShiftList->at(partitionIndex),
                            indexVector[actualLCPIndex - 1]);
                    newMatchLocation->InsertLocation(SAVector[actualLCPIndex] + partitionsShiftList->at(partitionIndex),
                                                     indexVector[actualLCPIndex]);
                    std::pair<std::string, MatchLocations> newMatchPair = std::make_pair(*partitions->at(partitionIndex),*newMatchLocation);
                    matchesMap.insert(newMatchPair);
                }
            }
        }
        partitions->clear();
        partitionsShiftList->clear();
    }
}

void
Determine_Valid_Matches_Child(std::unordered_map<std::string, MatchLocations> &matchesMap, std::vector<int> &LCPVector,
                              std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
                              int maximumMatchSize, std::vector<std::shared_ptr<std::string>> &seqStringVector,
                              size_t startIndex, size_t endIndex) {
    size_t index = startIndex;
    size_t hardEndIndex = SAVector.size();

    std::shared_ptr<MatchLocations> newMatchValidity;
    std::pair<std::string, MatchLocations> newMatchPair;
    std::shared_ptr<std::vector<std::shared_ptr<std::string>>> partitions;

    while(index < endIndex){
        if (LCPVector[index] >= minimumMatchSize){
            while((index  < hardEndIndex) && (LCPVector[index] >= minimumMatchSize)){ //Determine end of run (Go Past end if match overlaps)
                //For each substring, partition string by taking prefix of suffix
                //A given substring will always belong to one string as idenfied by LCP value, since LCP will not match beyond a string due to sentinels.
                //Therefore, all of a substrings partitions will also belong to its parent string.
                partitions = Determine_StringPartitions(seqStringVector[indexVector[index]]->substr(SAVector[index], LCPVector[index]),
                                                        LCPVector[index], minimumMatchSize, maximumMatchSize);
                //Check all partitions against map
                for (auto & partition : *partitions) { //Iterate all partitions.
                    if(matchesMap.count(*partition)){ //If partition exists, add match location
                        //Add both sequence instances (i.e: Where each sequence belongs)
                        matchesMap[*partition].InsertUniqueSeq(indexVector[index - 1]);
                        matchesMap[*partition].InsertUniqueSeq(indexVector[index]);
                        //Add only singular LCP index as all necessary information is given from single index
                        //I.E: We can extract both locations of a valid match from SA[index] and SA[index-1]
                        matchesMap[*partition].Insert_LCPIndex(index);
                    }
                    else{ //Create new possible match
                        newMatchValidity = std::make_shared<MatchLocations>();
                        //Add both sequence instances (i.e: Where each sequence belongs)
                        newMatchValidity->InsertUniqueSeq(indexVector[index - 1]);
                        newMatchValidity->InsertUniqueSeq(indexVector[index]);
                        //Add only singular LCP index as all necessary information is given from single index
                        //I.E: We can extract both locations of a valid match from SA[index] and SA[index-1]
                        newMatchValidity->Insert_LCPIndex(index);
                        newMatchPair = std::make_pair(*partition, *newMatchValidity);
                        matchesMap.insert(newMatchPair);
                    }
                }
                partitions->clear();
                ++index;
            }
        }
        ++index;
    }

}






std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_StringPartitions(const std::string &key, const int &keyLen, const int &minLength, const int &maxLength) {
    //Determines all possible substrings of a given key, that are length = [min,max]
    std::vector<std::shared_ptr<std::string>> partitionsStringList;
    int partitionLength = keyLen < maxLength ? keyLen : maxLength;
    for(; partitionLength >= minLength; --partitionLength){
        for (int partitionIndex = 0; (partitionIndex + partitionLength) <= keyLen ; ++partitionIndex){
            partitionsStringList.push_back(std::make_shared<std::string>(std::string(key.substr(partitionIndex, partitionLength))));
        }
    }

    return std::make_shared<std::vector<std::shared_ptr<std::string>>>(partitionsStringList);
}



std::shared_ptr<std::vector<std::shared_ptr<std::string>>>
Determine_LocationPartitions(const std::string &key, const int &keyLen, const int &minLength, const int &maxLength,
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
        matchVector = *(match.second.GetMatchVector());
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
