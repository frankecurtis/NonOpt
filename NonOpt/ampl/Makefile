# Copyright (C) 2019 Frank E. Curtis
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
NonOptLIB = $(NONOPTDIR)/NonOpt/src/libNonOpt.a
AMPLLIB   = $(AMPLDIR)/amplsolver.a

# Includes
INCLUDES = -I $(NONOPTDIR)/NonOpt/src -I $(AMPLDIR) -I /usr/local/include

# Rule for all
all: $(EXES)

# Dependence for executable
$(EXES): % : %.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $< $(NonOptLIB) $(AMPLLIB) -L $(LAPACKDIR) -ldl -lblas -llapack

# Dependencies for executable
$(EXES): $(NonOptLIB)

# Rules for libraries
$(NonOptLIB):
	$(MAKE) --directory=$(NONOPTDIR)/NonOpt/src

# Dependencies for libraries
$(NonOptLIB): $(wildcard $(NONOPTDIR)/NonOpt/src/*.hpp) $(wildcard $(NONOPTDIR)/NonOpt/src/*.cpp)

# Rule for .cpp.o
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
