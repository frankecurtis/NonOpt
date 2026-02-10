# Copyright (C) 2026 Lara Zebiane, Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Lara Zebiane, Frank E. Curtis

# C++ compiler
CXX = g++

# C++ compiler flags
CXXFLAGS = -Wall -std=c++11 -shared -fPIC

# Python flags
PYINCLUDES = $(shell python3 -m pybind11 --includes)
PYSUFFIX   = $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

# Libraries
NonOptLIB = "$(NONOPTDIR)"/NonOpt/src/libNonOpt.a

# Includes
INCLUDES = -I "$(NONOPTDIR)"/NonOpt/src -I $(PYBINDDIR) $(PYINCLUDES)

# Target
EXE = nonoptpy$(PYSUFFIX)

# Rule for all
all: $(EXE)

# Dependency for executable
$(EXE): $(NonOptLIB)

# Rules for libraries
$(NonOptLIB):
	$(MAKE) --directory="$(NONOPTDIR)"/NonOpt/src

# Dependencies for libraries
$(NonOptLIB): $(wildcard "$(NONOPTDIR)"/NonOpt/src/*.hpp) $(wildcard "$(NONOPTDIR)"/NonOpt/src/*.cpp)

# Rule for executable
$(EXE): PythonProblem.cpp $(NonOptLIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) PythonProblem.cpp -o $(EXE) -undefined dynamic_lookup $(NonOptLIB) -ldl -lblas -llapack

# Clean
clean:
	rm -f nonoptpy*.so

# Very clean
veryclean: clean