cmake_minimum_required(VERSION 3.17)
project(SequenceMatching_LCP)

set(CMAKE_CXX_STANDARD 11)

add_executable(SequenceMatching_LCP main.cpp
        Sais.c
        FileManipulation.cpp
        FileManipulation.h
        MatchLocations.cpp
        MatchLocations.h SequenceMatching.cpp SequenceMatching.h PossibleMatches.cpp PossibleMatches.h)