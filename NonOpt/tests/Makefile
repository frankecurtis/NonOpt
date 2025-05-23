# Copyright (C) 2025 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Frank E. Curtis

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Set sources, etc.
headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)

# Executable(s)
EXES = $(sources:.cpp=)

# Libraries
NonOptLIB = "$(NONOPTDIR)"/NonOpt/src/libNonOpt.a
NonOptProblemsLIB = "$(NONOPTDIR)"/NonOpt/problems/libNonOptProblems.a

# Includes
INCLUDES = -I "$(NONOPTDIR)"/NonOpt/src -I "$(NONOPTDIR)"/NonOpt/problems

# Rule for all
all: $(EXES)

# Rule for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(NonOptLIB) $(NonOptProblemsLIB) -ldl -lblas -llapack

# Dependencies for executable
$(EXES): $(NonOptLIB) $(NonOptProblemsLIB)

# Rules for libraries
$(NonOptLIB):
	$(MAKE) --directory="$(NONOPTDIR)"/NonOpt/src
$(NonOptProblemsLIB):
	$(MAKE) --directory="$(NONOPTDIR)"/NonOpt/problems

# Dependencies for libraries
$(NonOptLIB): $(wildcard "$(NONOPTDIR)"/NonOpt/src/*.hpp) $(wildcard "$(NONOPTDIR)"/NonOpt/src/*.cpp)
$(NonOptProblemsLIB): $(wildcard "$(NONOPTDIR)"/NonOpt/problems/*.hpp) $(wildcard "$(NONOPTDIR)"/NonOpt/problems/*.cpp)

# Rule for objects
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Objects depend on headers and sources
$(objects): $(headers)
$(objects): $(sources)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(EXES)
