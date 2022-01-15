CXX = g++
CC = gcc
CFLAGS = -c -std=c11 -Wall -O3
CXXFLAGS = -c -std=c++11 -Wall -O3
CSOURCES = Sais.c
CPPSOURCES = main.cpp FileManipulation.cpp MatchLocations.cpp SequenceMatching.cpp
COBJECTS = $(CSOURCES:.c=.o)                              
CPPOBJECTS = $(CPPSOURCES:.cpp=.o) 
EXECUTABLE = SequenceMatching_LCP

.PHONY : all
all: $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS) $(CPPOBJECTS)
	$(CXX) $(COBJECTS) $(CPPOBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY : clean
clean:
		rm -rf $(COBJECTS)
		rm -rf $(CPPOBJECTS)
		rm -rf $(EXECUTABLE)