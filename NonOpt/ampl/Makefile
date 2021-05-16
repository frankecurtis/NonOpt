# Copyright (C) 2019 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Frank E. Curtis

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Library utility command
AR = ar rv

# Library
LIB = libAMPLProblem.a

# Executable
EXE = solveAMPLProblem

# Libraries
NonOptLIB = "$(NONOPTDIR)"/NonOpt/src/libNonOpt.a
AMPLLIB   = "$(AMPLDIR)"/amplsolver.a

# Includes
INCLUDES = -I "$(NONOPTDIR)"/NonOpt/src -I "$(AMPLDIR)"

# Rule for all
all: $(LIB) $(EXE)

# Create library
$(LIB): AMPLProblem.o
	$(AR) $@ $<

# Rule for object
AMPLProblem.o: AMPLProblem.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Dependence
AMPLProblem.cpp: AMPLProblem.hpp

# Create executable
$(EXE): solveAMPLProblem.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(NonOptLIB) $(AMPLLIB) $(LIB) -L "$(LAPACKDIR)" -ldl -lblas -llapack

# Rule for object
solveAMPLProblem.o: solveAMPLProblem.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Dependencies
$(LIB): $(NonOptLIB)
$(EXE): $(NonOptLIB)

# Rules
$(NonOptLIB): $(wildcard "$(NONOPTDIR)"/NonOpt/src/*.hpp) $(wildcard "$(NONOPTDIR)"/NonOpt/src/*.cpp)
	$(MAKE) --directory="$(NONOPTDIR)"/NonOpt/src

# Clean
clean:
	rm -f AMPLProblem.o solveAMPLProblem.o

# Very clean
veryclean: clean
	rm -f $(LIB) $(EXE)
