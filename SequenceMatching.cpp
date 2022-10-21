//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"
#include "boost/xpressive/xpressive.hpp"


std::shared_ptr<std::vector<int>> Determine_Index_Mapping(std::vector<int> &SAvector, std::vector<int> &seqRangeArray){
    //Determines mapping between SA index and sequence it belongs to.
    int shift;
    int SAIndex;
    int rangeIndex;
    std::vector<int> indexVector;

    for(SAIndex = 0; SAIndex < SAvector.size(); ++SAIndex) {
        rangeIndex = 0;
        //Don't need to check boundaries since index guaranteed to be less than last item.
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
//   -For all partitions, determine if already exists, if so, add their sequence index.
//2) Merge all matches found by threads, specifically their sequence index.
//3) Determine if a given match includes sequences from all indecies.
//    - If not remove the match
//4) Create new child threads to iterate through valid remaining matches using regex to find all instances of a given
//   match and store their indecies.


std::shared_ptr<std::unordered_map<std::string, MatchLocations>> Determine_Possible_Matches_Parent(
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
    std::vector<std::shared_ptr<std::unordered_map<std::string, MatchValidity>>> validityThreadVector;

    //Todo change print messages
    //////////////////////////////////////////////////////////////////////////
    //1. Remove invalid matches, by determining number of sequences present.//
    //////////////////////////////////////////////////////////////////////////
    //Initialize Vector

    std::cout << "Determining Valid Matches..." << std::endl;

    for (i = 0; i < numSplits; ++i) {
        std::shared_ptr<std::unordered_map<std::string, MatchValidity>> matchValidityMap = std::make_shared<std::unordered_map<std::string, MatchValidity>>();
        validityThreadVector.push_back(matchValidityMap);
    }

    //Create 1st Thread
    std::cout << "Creating Validity Thread: 0" << std::endl;
    threadArray[0] = std::thread(Determine_Valid_Matches_Child, std::ref(*validityThreadVector[0]), std::ref(LCPVector),
                                 std::ref(SAVector), std::ref(indexVector), minimumMatchSize, maximumMatchSize,
                                 std::ref(seqStringVector), std::ref(seqRangeVector), startIndex, endIndex);

//    Determine_Valid_Matches_Child(*validityThreadVector[0], LCPVector, SAVector, indexVector, minimumMatchSize,
//                                  maximumMatchSize, seqStringVector, seqRangeVector, startIndex,
//                                  endIndex);  //Single Threaded for Testing

    //Create subsequent threads
    for (i = 1; i < numSplits; ++i) {
        startIndex = endIndex;
        endIndex = startIndex+splitSize;
        std::cout << "Creating Validity Thread: "<< i << std::endl;
        threadArray[i] = std::thread(Determine_Valid_Matches_Child, std::ref(*validityThreadVector[i]), std::ref(LCPVector),
                                     std::ref(SAVector), std::ref(indexVector), minimumMatchSize, maximumMatchSize,
                                     std::ref(seqStringVector), std::ref(seqRangeVector), startIndex, endIndex);
//        Determine_Valid_Matches_Child(*validityThreadVector[i], LCPVector, SAVector, indexVector, minimumMatchSize,
//                                      maximumMatchSize, seqStringVector, seqRangeVector, startIndex,
//                                      endIndex); //Single Threaded for Testing
    }


    //Wait for all threads to complete
    for (i=0; i<numSplits; i++){
        threadArray[i].join();
    }

    std::cout << "Finished Determining Valid Matches." << std::endl;

    //Merge all thread validity maps into parentMatchValidity
    std::unordered_map<std::string, MatchValidity> parentMatchValidity;
    std::vector<std::string> keysToCheckValidity;

    for (threadMapIndex = 0; threadMapIndex < numSplits; ++threadMapIndex) {
        std::cout << "Merging Validity Thread: " << threadMapIndex << std::endl;
        for(auto &threadMapMatch: *validityThreadVector[threadMapIndex]) { //For each match in threadMap
            if(parentMatchValidity.count(threadMapMatch.first)){ //If match already exists, add sequence indecies
                parentMatchValidity[threadMapMatch.first].mergeSeqIndex(threadMapMatch.second.getUniqueSeqIndices());
            }
            else{ //Create new pair and push to map
                keysToCheckValidity.push_back(threadMapMatch.first); //Creating ma
                std::pair<std::string, MatchValidity> newMatchPair = std::make_pair(threadMapMatch.first, threadMapMatch.second);
                parentMatchValidity.insert(newMatchPair);
            }
        }
        validityThreadVector[threadMapIndex]->clear();
    }
    validityThreadVector.clear();

    std::cout <<  "Removing invalid matches..." << std::endl;
    //Check if key has items from all sequences.
    std::vector<std::string> validMatches; //Contains valid matches to be used in splitting matches
    for(auto & key: keysToCheckValidity){
        if(parentMatchValidity[key].numSeqIncluded() == numSequences){
            validMatches.push_back(key);
        }
    }
    //No Longer needed, only need match strings to search locations.
    parentMatchValidity.clear();
    keysToCheckValidity.clear();

    //////////////////////////////////////////////////////////
    //2. Use threads to determine locations of valid matches//
    //////////////////////////////////////////////////////////
    //Split length of validMatches amongst each thread. Each thread then proceeds to
    //Check all of the matches in its share and add the
    splitSize = validMatches.size() / numSplits;
    splitRemainder = validMatches.size() % numSplits; //Determine remainder, and add to first task.
    startIndex = 0;
    endIndex = splitSize + splitRemainder;

    std::vector<std::shared_ptr<std::unordered_map<std::string, MatchLocations>>> threadMapVector;

    std::cout << "Determining Match Locations" << std::endl;

    //Initialize Vector
    for (i = 0; i < numSplits; ++i) {
        std::shared_ptr<std::unordered_map<std::string, MatchLocations>> matchLocationMap = std::make_shared<std::unordered_map<std::string, MatchLocations>>(numSequences);
        threadMapVector.push_back(matchLocationMap);
    }

    //Create 1st Thread
    std::cout << "Creating Location Thread: 0" << std::endl;
    threadArray[0] = std::thread(Determine_Match_Locations_Child, std::ref(*threadMapVector[0]),
                                 std::ref(validMatches), std::ref(indexVector), std::ref(seqStringVector),
                                 numSequences,startIndex, endIndex);

//    Determine_Match_Locations_Child(*threadMapVector[0],validMatches,indexVector,seqStringVector,numSequences,startIndex, endIndex); //Single Threaded for Testing

    //Create subsequent threads
    for (i = 1; i < numSplits; ++i) {
        startIndex = endIndex;
        endIndex = startIndex + splitSize;
        std::cout << "Creating Location Thread: " << i << std::endl;
        threadArray[i] = std::thread(Determine_Match_Locations_Child, std::ref(*threadMapVector[i]),
                                     std::ref(validMatches), std::ref(indexVector), std::ref(seqStringVector),
                                     numSequences,startIndex, endIndex);

//        Determine_Match_Locations_Child(*threadMapVector[i],validMatches,indexVector,seqStringVector,numSequences,startIndex, endIndex); //Single Threaded for Testing
    }

    //Wait for all threads to complete
    for (i=0; i<numSplits; i++){
        threadArray[i].join();
    }

    std::cout << "Finished Determining Match Locations." << std::endl;

    //Merge all thread match locations into parent match location map (In this case we are using the map stored in
    //threadMapVector[0] as the parent to avoid re-allocation of memory
    for (i = 1; i < numSplits; ++i) {
        std::cout << "Merging Match Locations: " << i-1 << std::endl;
        threadMapVector[0]->insert(threadMapVector[i]->begin(), threadMapVector[i]->end());
        threadMapVector[i]->clear();
    }

    return threadMapVector[0];
}


void
Determine_Valid_Matches_Child(std::unordered_map<std::string, MatchValidity> &matchesMap, std::vector<int> &LCPVector,
                              std::vector<int> &SAVector, std::vector<int> &indexVector, int minimumMatchSize,
                              int maximumMatchSize, std::vector<std::shared_ptr<std::string>> &seqStringVector,
                              std::vector<int> &seqRangeVector, size_t startIndex, size_t endIndex) {
    size_t index = startIndex;
    size_t hardEndIndex = SAVector.size();
    std::shared_ptr<std::vector<std::unordered_set<int>>> matchVector;
    std::shared_ptr<std::string> newMatchString;
    std::vector<std::string> matchesMapKeys;

    std::shared_ptr<MatchValidity> newMatchValidity;
    std::pair<std::string, MatchValidity> newMatchPair;
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
                        matchesMap[*partition].InsertSeqIndex(indexVector[index - 1]);
                        matchesMap[*partition].InsertSeqIndex(indexVector[index]);
                    }
                    else{ //Create new possible match
                        newMatchValidity = std::make_shared<MatchValidity>();
                        newMatchValidity->InsertSeqIndex(indexVector[index - 1]);
                        newMatchValidity->InsertSeqIndex(indexVector[index]);
                        newMatchPair = std::make_pair(*partition, *newMatchValidity);
                        matchesMap.insert(newMatchPair);
                        matchesMapKeys.push_back(*partition);
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


void Determine_Match_Locations_Child(std::unordered_map<std::string, MatchLocations> &matchesMap,
                                     std::vector<std::string> &validMatches, std::vector<int> &indexVector,
                                     std::vector<std::shared_ptr<std::string>> &seqStringVector, int numSequences,
                                     size_t startIndex, size_t endIndex) {
    std::string matchString;
    size_t foundAtIndex;
    size_t matchesMapIndex = startIndex;
    boost::xpressive::sregex regexKey;
    boost::xpressive::sregex_iterator foundIterator();




    std::shared_ptr<MatchLocations> newMatchLocations;

    while(matchesMapIndex < endIndex){ //Iterate through all matches in partition
        matchString = validMatches[matchesMapIndex];
        regexKey = boost::xpressive::sregex::compile("(?=(" + matchString + ").)"); //Added positive lookahead to find overlapping matchStrings.
        newMatchLocations = std::make_shared<MatchLocations>(numSequences);
        auto newMatchPair = std::make_pair(matchString, *newMatchLocations);
        matchesMap.insert(newMatchPair);



        for (size_t seqIndex = 0; seqIndex < numSequences; ++seqIndex) { //For each sequence in query, determine instances of match
            foundIterator = std::sregex_iterator(seqStringVector[seqIndex]->begin(), seqStringVector[seqIndex]->end(), regexKey);
            for(std::sregex_iterator i = foundIterator; i != regexEnd; i++){
                foundAtIndex = i->position();
                matchesMap.at(matchString).InsertMatch(foundAtIndex,seqIndex);
            }
        }
        ++matchesMapIndex;
    }

}





std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap, std::string &seqStringCombined,
                            std::vector<std::shared_ptr<std::string>> &seqStringVector,
                            std::vector<int> &seqRangeVector, int &numSequences) {
    int seqIndex;
    int matchSet;
    size_t matchIndices;
    int overallCoverageSize = 0;
    std::vector<double> similarityMetricVector;
    std::vector<std::unordered_set<size_t>> matchVector;
    std::vector<std::unordered_set<size_t>> matchCoverageVector;
    for (seqIndex = 0; seqIndex < numSequences; ++seqIndex) {
        matchCoverageVector.emplace_back(std::unordered_set<size_t>());
    }

    for(auto &match: matchesMap){
        matchVector = *(match.second.getMatchVector());
        for(matchSet = 0; matchSet < numSequences; ++matchSet){
            for(auto &matchIndex: matchVector[matchSet]){
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
