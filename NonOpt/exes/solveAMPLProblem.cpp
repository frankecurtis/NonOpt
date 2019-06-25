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

  // Open file
  char nl_file[100];
  strcpy(nl_file, (char*)"cute_unconstrained/");
  strcat(nl_file, argv[1]);
  strcat(nl_file, (char*)".nl");

  // Declare pointer to problem
  std::shared_ptr<AMPLProblem> problem = std::make_shared<AMPLProblem>(nl_file);

  // Declare algorithm
  NonOptSolver nonopt;

  // Modify direction computation option
  nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", argv[2]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "inverse_hessian_update", argv[3]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "line_search", argv[4]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "point_set_update", argv[5]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "qp_solver", argv[6]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "symmetric_matrix", argv[7]);

  // Declare file report
  std::shared_ptr<FileReport> r(new FileReport("f", R_NL, R_PER_INNER_ITERATION));

  // Open file
  char out_file[100];
  strcpy(out_file, (char*)"output/");
  strcat(out_file, argv[1]);
  strcat(out_file, (char*)"_");
  strcat(out_file, argv[2]);
  strcat(out_file, (char*)"_");
  strcat(out_file, argv[3]);
  strcat(out_file, (char*)"_");
  strcat(out_file, argv[4]);
  strcat(out_file, (char*)"_");
  strcat(out_file, argv[5]);
  strcat(out_file, (char*)"_");
  strcat(out_file, argv[6]);
  strcat(out_file, (char*)"_");
  strcat(out_file, argv[7]);
  strcat(out_file, (char*)".out");
  r->open(out_file);

  // Add to reporter
  nonopt.reporter()->addReport(r);

  // Optimize
  nonopt.optimize(problem);

  // Return
  return 0;

}  // end main
