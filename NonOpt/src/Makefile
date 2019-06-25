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

# Set sources and objects
sources = \
NonOptDirectionComputationCuttingPlane.cpp \
NonOptDirectionComputationGradient.cpp \
NonOptDirectionComputationGradientCombination.cpp \
NonOptInverseHessianUpdateBFGS.cpp \
NonOptLineSearchBacktracking.cpp \
NonOptLineSearchWeakWolfe.cpp \
NonOptOptions.cpp \
NonOptPoint.cpp \
NonOptPointSetUpdateProximity.cpp \
NonOptQuantities.cpp \
NonOptQPSolverActiveSet.cpp \
NonOptReporter.cpp \
NonOptStrategies.cpp \
NonOptSymmetricMatrixDense.cpp \
NonOptSymmetricMatrixLimitedMemory.cpp \
NonOptVector.cpp \
NonOptSolver.cpp
objects = $(sources:.cpp=.o)
dependlists = $(sources:.cpp=.d)

# Library
LIB = libNonOpt.a

# Create library
$(LIB): $(objects)
	$(AR) $@ $(objects)

# Rule for .cpp.o
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(objects) $(dependlists)

# Very clean
veryclean: clean
	rm -f $(LIB)
