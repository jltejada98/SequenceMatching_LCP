//
// Created by Jose Luis Tejada on 5/8/21.
//

#include "SequenceMatching.h"


std::shared_ptr<std::unordered_map<std::string, MatchLocations>>
Determine_Matches(std::vector<int> &LCPVector,std::vector<int> &SAvector, std::shared_ptr<SequenceRange> *seqRangeArray,
                  std::shared_ptr<std::string> *seqStringArray,size_t &minimumMatchSize, size_t &numSequences){
    size_t LCPVectorSize = LCPVector.size();
    size_t windowIndex;
    size_t index = 1; //LCP[0] == 0 (Always)
    size_t seqRangeIndex;
    size_t matchLength;
    size_t matchIndexForString;
    std::unordered_map<std::string, MatchLocations> matchesMap;
    std::shared_ptr<std::unordered_set<size_t>> matchSet;
    std::shared_ptr<std::string> newMatchString;
    std::shared_ptr<MatchLocations> newMatchLocation;
    std::pair<std::string,MatchLocations> newMatchPair;
    std::shared_ptr<std::vector<std::shared_ptr<std::unordered_set<size_t>>>> existingMatchVector;
    std::shared_ptr<std::unordered_set<size_t>> transferredSet;

    //Todo (Consider doing with threads) Split LCP index into n parts and process in parallel (See Previous Implementation)
    // (Would need locks to ensure only 1 thread write to map at a time.)
    while(index < LCPVectorSize){
        if (LCPVector[index] >= minimumMatchSize){
            windowIndex = index + 1;
            //Determine number of possible matches
            while(LCPVector[windowIndex] >= LCPVector[index]){ //LCP values must be equal for match to be the same (Adjacent)
                ++windowIndex;
            }

            //Determine if enough possible matches were found.
            if ((windowIndex - index) >= numSequences){
                //Create Index Sets
                if(!seqRangeArray[0]->checkSetAllocation()){ //Check if Index Sets are allocated
                    for (seqRangeIndex = 0; seqRangeIndex < numSequences ; ++seqRangeIndex) {
                        matchSet = std::make_shared<std::unordered_set<size_t>>();
                        seqRangeArray[seqRangeIndex]->insertSet(matchSet);
                    }
                }

                //Determine if the read indices include matches from all strings.
                for (size_t saIndex = (index - 1); saIndex < windowIndex; ++saIndex) { //Start at index-1 Since LCP Compares Adjacent Values
                    seqRangeIndex = 0;
                    while(!seqRangeArray[seqRangeIndex]->checkInsertIndex(SAvector[saIndex]) ){ //Not necesary to check index length (Guaranteed to exist)
                        ++seqRangeIndex;
                    }
                }

                //Check if all seqRangeArrays are true
                seqRangeIndex = 0;
                while((seqRangeArray[seqRangeIndex]->getCheck()) && (seqRangeIndex < numSequences)){
                    ++seqRangeIndex;
                }

                if (seqRangeIndex < numSequences){ //Match is not included in all sequences.
                    //Reset range checks and clear sets
                    for (seqRangeIndex = 0;  seqRangeIndex < numSequences ; ++seqRangeIndex) {
                        seqRangeArray[seqRangeIndex]->reset();
                    }
                }
                else{//Add Matches and Locations to Matches Map.
                    //Obtain element index from first sequence index.
                    matchLength = LCPVector[index];
                    matchIndexForString = *seqRangeArray[0]->getSet()->begin();
                    newMatchString = std::make_shared<std::string>(std::string(seqStringArray[0]->substr(matchIndexForString, matchLength)));

                    if(matchesMap.count(*newMatchString)){ //Check if match already exists
                        existingMatchVector = matchesMap.at(*newMatchString).getMatchVector();
                        for (seqRangeIndex = 0; seqRangeIndex < numSequences; ++seqRangeIndex) {
                            transferredSet = seqRangeArray[seqRangeIndex]->transferSet();
                            existingMatchVector->at(seqRangeIndex)->insert(transferredSet->begin(), transferredSet->end());
                        }
                    }
                    else{ //Insert new match
                        newMatchLocation = std::make_shared<MatchLocations>();
                        //Transfer set pointers
                        for (seqRangeIndex = 0;  seqRangeIndex < numSequences ; ++seqRangeIndex) {
                            newMatchLocation->insertSet(seqRangeArray[seqRangeIndex]->transferSet()); //Transferring sets resets range checks.
                        }
                        newMatchPair = std::make_pair(*newMatchString, *newMatchLocation);
                        matchesMap.insert(newMatchPair);
                    }
                }
            }
            else{ //Reset range checks and clear sets
                for (seqRangeIndex = 0;  seqRangeIndex < numSequences ; ++seqRangeIndex) {
                    seqRangeArray[seqRangeIndex]->reset();
                }
            }
        }

        ++index;
    }

    return std::make_shared<std::unordered_map<std::string, MatchLocations>>(matchesMap);
}


std::shared_ptr<std::vector<double>>
Determine_SimilarityMetrics(std::unordered_map<std::string, MatchLocations> &matchesMap,std::string &seqStringCombined,
                            std::shared_ptr<std::string> *seqStringArray,
                            size_t &numSequences){
    size_t index;
    size_t seqIndex;
    size_t matchLength;
    size_t overallCoverageSize = 0;
    std::vector<double> similarityMetricVector;
    std::vector<std::unordered_set<size_t>> seqVectorSet;
    std::vector<std::shared_ptr<std::unordered_set<size_t>>> matchIndeciesVector;
    for (seqIndex = 0; seqIndex < numSequences; ++seqIndex) {
        seqVectorSet.emplace_back(std::unordered_set<size_t>());
    }
    for (auto &match: matchesMap){
        matchIndeciesVector = *match.second.getMatchVector();
        matchLength = match.first.length();
        for (seqIndex = 0; seqIndex < numSequences ; ++seqIndex) {
            //For each index in each sequence set.
            for (auto &seqSetIndex : *matchIndeciesVector[seqIndex]){
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


void Determine_Submatching(std::unordered_map<std::string, MatchLocations> &matchesMap, size_t &minimumMatchSize,
                           size_t &numSequences){
    //Todo (Consider doing with threads) Split MatchMap into n parts and process in parallel
    std::string matchString;
    MatchLocations matchLocations;
    size_t matchLength;
    size_t numPartitions;

    for(auto &matchToCheck: matchesMap){

        if(matchToCheck.first.length() > minimumMatchSize){ //Submatching condition
            matchString = matchToCheck.first;
            matchLocations = matchToCheck.second;
            matchLength = matchString.length();
            auto matchVector = *matchLocations.getMatchVector();
            numPartitions = ((matchLength - minimumMatchSize) * (matchLength - minimumMatchSize + 3)) / 2; //Closed form for number of partitions


            if (numPartitions <= matchesMap.size()){ //Partion match into all possible substrings and check if exists in map.
                //Determine partitions
                size_t partitionShift;
                std::vector<std::shared_ptr<std::string>> partitionsStringList;
                std::vector<std::shared_ptr<size_t>> partitionsShiftList;
                for(size_t p = matchLength - 1; p >= minimumMatchSize; --p){
                    for (size_t i = 0; (i+p) <= matchLength ; ++i){
                        partitionsShiftList.push_back(std::make_shared<size_t>(i));
                        partitionsStringList.push_back(std::make_shared<std::string>(std::string(matchString.substr(i, p))));
                    }
                }
                //Check if partions exist in map, if so add indecies of matchLocations into submatch
                for (size_t partitionIndex = 0; partitionIndex < partitionsShiftList.size(); ++partitionIndex) {
                    if(matchesMap.count(*partitionsStringList.at(partitionIndex)) >= 1){
                        partitionShift = *partitionsShiftList.at(partitionIndex);
                        for (size_t setIndex = 0; setIndex < numSequences; ++setIndex) {
                            for(auto &index: *matchVector.at(setIndex)){
                                matchesMap.at(*partitionsStringList.at(partitionIndex)).insertIndex(setIndex, index+partitionShift);
                            }
                        }
                    }
                }

            }
            else{ //Check all elements in hashtable against the keyToCheck.
                size_t foundAtIndex;
                std::regex regexKey;
                std::sregex_iterator foundIterator;
                auto end = std::sregex_iterator();

                for(auto &key: matchesMap){ //Iterate all keys
                    if ((key.first != matchString) && (matchString.length() > key.first.length()) && (matchString.find(key.first) != std::string::npos)){ //Determines if key is contained in matchString
                        regexKey = std::regex(std::string(key.first));
                        foundIterator = std::sregex_iterator(matchString.begin(), matchString.end(), regexKey);
                        for(std::sregex_iterator i = foundIterator; i != end; i++){ //For each instance of substring found.
                            std::smatch subtringMatch = *i;
                            foundAtIndex = i->position();
                            for (size_t setIndex = 0; setIndex < numSequences; ++setIndex) {
                                for(auto &index: *matchVector.at(setIndex)){
                                    matchesMap.at(matchString).insertIndex(foundAtIndex, index+foundAtIndex);
                                }
                            }
                        }
                    }
                }

            }

        }
    }



}
