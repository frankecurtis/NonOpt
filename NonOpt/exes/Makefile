# Copyright (C) 2019 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Frank E. Curtis

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -g -Wall -std=c++11

# Executable(s)
EXES = solveAMPLProblem solveNonsmoothProblem

# Sources, objects, etc.
sources     = $(EXES:=.cpp)
objects     = $(sources:.cpp=.o)
dependlists = $(sources:.cpp=.d)
dsyms       = $(sources:.cpp=.dSYM)

# Libraries
NonOptLIB         = $(NonOptDIR)/libNonOpt.a
NonOptProblemsLIB = $(NonOptProblemsDIR)/libNonOptProblems.a
AMPLLIB           = $(AMPLDIR)/amplsolver.a

# Includes
INCLUDES = -I $(NonOptDIR) -I $(NonOptProblemsDIR) -I $(AMPLDIR)

# Rule for all
all: $(EXES)

# Dependence for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(NonOptLIB) $(NonOptProblemsLIB) $(AMPLLIB) -L/Users/frankecurtis/Dropbox/software/lapack-3.8.0 -ldl -lblas -llapack

# Rule for .cpp.o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ -I/usr/local/include

# Clean
clean:
	rm -f $(objects) $(dependlists)
	rm -rf $(dsyms)

# Very clean
veryclean: clean
	rm -f $(EXES)
