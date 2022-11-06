cmake_minimum_required(VERSION 3.17)
project(SequenceMatching_LCP)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

#######################
######RE2 Library######
#######################

# ABI version http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html
set(SONAME 10)

set(RE2_SOURCES
        re2/bitstate.cc
        re2/compile.cc
        re2/dfa.cc
        re2/filtered_re2.cc
        re2/mimics_pcre.cc
        re2/nfa.cc
        re2/onepass.cc
        re2/parse.cc
        re2/perl_groups.cc
        re2/prefilter.cc
        re2/prefilter_tree.cc
        re2/prog.cc
        re2/re2.cc
        re2/regexp.cc
        re2/set.cc
        re2/simplify.cc
        re2/stringpiece.cc
        re2/tostring.cc
        re2/unicode_casefold.cc
        re2/unicode_groups.cc
        util/rune.cc
        util/strutil.cc
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/re2)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/util)

add_library(re2 ${RE2_SOURCES})
target_compile_features(re2 PUBLIC cxx_std_11)
target_include_directories(re2 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
set_target_properties(re2 PROPERTIES SOVERSION ${SONAME} VERSION ${SONAME}.0.0)
add_library(re2::re2 ALIAS re2)


add_executable(SequenceMatching_LCP main.cpp
        Sais.c
        FileManipulation.cpp
        FileManipulation.h
        MatchLocations.cpp
        MatchLocations.h
        SequenceMatching.cpp
        SequenceMatching.h
        MatchValidity.cpp
        MatchValidity.h)

target_link_libraries(SequenceMatching_LCP re2)