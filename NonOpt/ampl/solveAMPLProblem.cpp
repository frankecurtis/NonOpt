// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "AMPLProblem.hpp"
#include "NonOptSolver.hpp"

#include <cstdio>
#include <cstring>
#include <memory>

using namespace NonOpt;

int main(int argc, char* argv[])
{

  // Set usage string
  std::string usage("Usage: ./solveAMPLProblem ProblemName\n"
                    "       where ProblemName is name of nl file in ampl subdirectory.\n");

  // Check number of input arguments
  if (argc < 2) {
    printf("Too few arguments. Quitting.\n");
    printf("%s",usage.c_str());
    return 1;
  }

  // Open file
  char nl_file[100];
  strcat(nl_file, argv[1]);
  strcat(nl_file, (char*)".nl");

  // Declare pointer to problem
  std::shared_ptr<AMPLProblem> problem = std::make_shared<AMPLProblem>(nl_file);

  // Declare algorithm
  NonOptSolver nonopt;

  // Optimize
  nonopt.optimize(problem);

  // Return
  return 0;

}  // end main
