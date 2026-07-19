# Copyright (C) 2025 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Frank E. Curtis

# C++ compiler
CXX = g++

# C++ compiler flags (optimized release build by default)
CXXFLAGS = -Wall -std=c++11 -O3 -DNDEBUG

# Debug build flags (unoptimized, with symbols; ASSERT_EXCEPTION is unaffected
# by -DNDEBUG either way since it is a custom macro, not <cassert>)
DEBUG_CXXFLAGS = -Wall -std=c++11 -g -O0

# Library utility command
AR = ar rv

# Set sources, etc.
headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)

# Library
LIB = libNonOpt.a

# Rule for all
all: $(LIB)

# Debug build (use: make debug)
debug: CXXFLAGS := $(DEBUG_CXXFLAGS)
debug: $(LIB)

# Create library
$(LIB): $(objects)
	$(AR) $@ $(objects)

# Rule for objects
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Objects depend on headers and sources
$(objects): $(headers)
$(objects): $(sources)

# Clean
clean:
	rm -f $(objects) $(depends)

# Very clean
veryclean: clean
	rm -f $(LIB)
