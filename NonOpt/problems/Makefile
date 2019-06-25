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

# Sources, objects, etc.
sources = \
ActiveFaces.cpp \
AMPLProblem.cpp \
BrownFunction_2.cpp \
ChainedCB3_1.cpp \
ChainedCB3_2.cpp \
ChainedCrescent_1.cpp \
ChainedCrescent_2.cpp \
ChainedLQ.cpp \
ChainedMifflin_2.cpp \
MaxQ.cpp \
MxHilb.cpp \
QuadPoly.cpp
objects = $(sources:.cpp=.o)
dependlists = $(sources:.cpp=.d)

# Library
LIB = libNonOptProblems.a

# Create library
$(LIB): $(objects)
	$(AR) $@ $(objects)

# Includes
INCLUDES = -I $(NonOptDIR) -I $(AMPLDIR)

# Rule for .cpp.o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -f $(objects) $(dependlists)

# Very clean
veryclean: clean
	rm -f $(LIB)
