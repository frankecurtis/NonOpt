// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cstdio>
#include <cstring>
#include <string>

#include "ActiveFaces.hpp"
#include "BrownFunction_2.hpp"
#include "ChainedCB3_1.hpp"
#include "ChainedCB3_2.hpp"
#include "ChainedCrescent_1.hpp"
#include "ChainedCrescent_2.hpp"
#include "ChainedLQ.hpp"
#include "ChainedMifflin_2.hpp"
#include "MaxQ.hpp"
#include "MxHilb.hpp"
#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"
#include "QuadPoly.hpp"
#include "Test29_11.hpp"
#include "Test29_13.hpp"
#include "Test29_17.hpp"
#include "Test29_19.hpp"
#include "Test29_2.hpp"
#include "Test29_20.hpp"
#include "Test29_22.hpp"
#include "Test29_24.hpp"
#include "Test29_5.hpp"
#include "Test29_6.hpp"

using namespace NonOpt;

// Main function
int main(int argc, char* argv[])
{

  // Set usage string
  std::string usage("Usage: ./solveProblem ProblemName\n"
                    "       where ProblemName is name of problem in problems subdirectory.\n");

  // Check number of input arguments
  if (argc < 2) {
    printf("Too few arguments. Quitting.\n");
    printf("%s", usage.c_str());
    return 1;
  }

  // Declare problem dimension
  int const dimension = 1000;

  // Declare problem
  std::shared_ptr<Problem> problem;
  if (strcmp(argv[1], "ActiveFaces") == 0) {
    problem = std::make_shared<ActiveFaces>(dimension);
  }
  else if (strcmp(argv[1], "BrownFunction_2") == 0) {
    problem = std::make_shared<BrownFunction_2>(dimension);
  }
  else if (strcmp(argv[1], "ChainedCB3_1") == 0) {
    problem = std::make_shared<ChainedCB3_1>(dimension);
  }
  else if (strcmp(argv[1], "ChainedCB3_2") == 0) {
    problem = std::make_shared<ChainedCB3_2>(dimension);
  }
  else if (strcmp(argv[1], "ChainedCrescent_1") == 0) {
    problem = std::make_shared<ChainedCrescent_1>(dimension);
  }
  else if (strcmp(argv[1], "ChainedCrescent_2") == 0) {
    problem = std::make_shared<ChainedCrescent_2>(dimension);
  }
  else if (strcmp(argv[1], "ChainedLQ") == 0) {
    problem = std::make_shared<ChainedLQ>(dimension);
  }
  else if (strcmp(argv[1], "ChainedMifflin_2") == 0) {
    problem = std::make_shared<ChainedMifflin_2>(dimension);
  }
  else if (strcmp(argv[1], "MaxQ") == 0) {
    problem = std::make_shared<MaxQ>(dimension);
  }
  else if (strcmp(argv[1], "MxHilb") == 0) {
    problem = std::make_shared<MxHilb>(dimension);
  }
  else if (strcmp(argv[1], "QuadPoly") == 0) {
    problem = std::make_shared<QuadPoly>(dimension, 10, 5, 10.0);
  }
  else if (strcmp(argv[1], "Test29_2") == 0) {
    problem = std::make_shared<Test29_2>(dimension);
  }
  else if (strcmp(argv[1], "Test29_5") == 0) {
    problem = std::make_shared<Test29_5>(dimension);
  }
  else if (strcmp(argv[1], "Test29_6") == 0) {
    problem = std::make_shared<Test29_6>(dimension);
  }
  else if (strcmp(argv[1], "Test29_11") == 0) {
    problem = std::make_shared<Test29_11>(dimension);
  }
  else if (strcmp(argv[1], "Test29_13") == 0) {
    problem = std::make_shared<Test29_13>(dimension);
  }
  else if (strcmp(argv[1], "Test29_17") == 0) {
    problem = std::make_shared<Test29_17>(dimension);
  }
  else if (strcmp(argv[1], "Test29_19") == 0) {
    problem = std::make_shared<Test29_19>(dimension);
  }
  else if (strcmp(argv[1], "Test29_20") == 0) {
    problem = std::make_shared<Test29_20>(dimension);
  }
  else if (strcmp(argv[1], "Test29_22") == 0) {
    problem = std::make_shared<Test29_22>(dimension);
  }
  else if (strcmp(argv[1], "Test29_24") == 0) {
    problem = std::make_shared<Test29_24>(dimension);
  }
  else {
    printf("Invalid problem name. Quitting.\n");
    printf("%s", usage.c_str());
    return 1;
  }

  // Declare solver object
  NonOptSolver nonopt;

  // Modify options from file
  nonopt.options()->modifyOptionsFromFile("nonopt.opt");

  // Optimize
  nonopt.optimize(problem);

  // Return
  return 0;

} // end main
