# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SequenceMatching_LCP.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/SequenceMatching_LCP.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/SequenceMatching_LCP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SequenceMatching_LCP.dir/flags.make

CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/flags.make
CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o: /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/main.cpp
CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o -MF CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o.d -o CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o -c /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/main.cpp

CMakeFiles/SequenceMatching_LCP.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SequenceMatching_LCP.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/main.cpp > CMakeFiles/SequenceMatching_LCP.dir/main.cpp.i

CMakeFiles/SequenceMatching_LCP.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SequenceMatching_LCP.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/main.cpp -o CMakeFiles/SequenceMatching_LCP.dir/main.cpp.s

CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o: CMakeFiles/SequenceMatching_LCP.dir/flags.make
CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o: /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/Sais.c
CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o: CMakeFiles/SequenceMatching_LCP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o"
	/Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o -MF CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o.d -o CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o -c /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/Sais.c

CMakeFiles/SequenceMatching_LCP.dir/Sais.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SequenceMatching_LCP.dir/Sais.c.i"
	/Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/Sais.c > CMakeFiles/SequenceMatching_LCP.dir/Sais.c.i

CMakeFiles/SequenceMatching_LCP.dir/Sais.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SequenceMatching_LCP.dir/Sais.c.s"
	/Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/Sais.c -o CMakeFiles/SequenceMatching_LCP.dir/Sais.c.s

CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/flags.make
CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o: /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/FileManipulation.cpp
CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o -MF CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o.d -o CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o -c /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/FileManipulation.cpp

CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/FileManipulation.cpp > CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.i

CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/FileManipulation.cpp -o CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.s

CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/flags.make
CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o: /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/MatchLocations.cpp
CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o -MF CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o.d -o CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o -c /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/MatchLocations.cpp

CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/MatchLocations.cpp > CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.i

CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/MatchLocations.cpp -o CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.s

CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/flags.make
CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o: /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/SequenceMatching.cpp
CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o: CMakeFiles/SequenceMatching_LCP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o -MF CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o.d -o CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o -c /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/SequenceMatching.cpp

CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/SequenceMatching.cpp > CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.i

CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/SequenceMatching.cpp -o CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.s

# Object files for target SequenceMatching_LCP
SequenceMatching_LCP_OBJECTS = \
"CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o" \
"CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o" \
"CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o" \
"CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o" \
"CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o"

# External object files for target SequenceMatching_LCP
SequenceMatching_LCP_EXTERNAL_OBJECTS =

SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/main.cpp.o
SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/Sais.c.o
SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/FileManipulation.cpp.o
SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/MatchLocations.cpp.o
SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/SequenceMatching.cpp.o
SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/build.make
SequenceMatching_LCP: CMakeFiles/SequenceMatching_LCP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable SequenceMatching_LCP"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SequenceMatching_LCP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SequenceMatching_LCP.dir/build: SequenceMatching_LCP
.PHONY : CMakeFiles/SequenceMatching_LCP.dir/build

CMakeFiles/SequenceMatching_LCP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SequenceMatching_LCP.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SequenceMatching_LCP.dir/clean

CMakeFiles/SequenceMatching_LCP.dir/depend:
	cd /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug /Users/joseluistejada/Documents/Other/Development/C++/SequenceMatching_LCP/cmake-build-debug/CMakeFiles/SequenceMatching_LCP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SequenceMatching_LCP.dir/depend

