# Copyright (C) 2025 Frank E. Curtis
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

# Set sources, etc.
headers = $(wildcard *.hpp)
sources = $(wildcard *.cpp)
objects = $(sources:.cpp=.o)
depends = $(sources:.cpp=.d)

# Library
LIB = libNonOpt.a

# Rule for all
all: $(LIB)

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
