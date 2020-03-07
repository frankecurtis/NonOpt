// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"

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

#include <cstdio>
#include <cstring>
#include <memory>

using namespace NonOpt;

// Main function
int main(int argc, char* argv[])
{
	int const p_size=200;
  // Declare problem
  std::shared_ptr<Problem> problem;
  if (strcmp(argv[1], "ActiveFaces") == 0) {
    problem = std::make_shared<ActiveFaces>(p_size);
  }
  else if (strcmp(argv[1], "BrownFunction_2") == 0) {
    problem = std::make_shared<BrownFunction_2>(p_size);
  }
  else if (strcmp(argv[1], "ChainedCB3_1") == 0) {
    problem = std::make_shared<ChainedCB3_1>(p_size);
  }
  else if (strcmp(argv[1], "ChainedCB3_2") == 0) {
    problem = std::make_shared<ChainedCB3_2>(p_size);
  }
  else if (strcmp(argv[1], "ChainedCrescent_1") == 0) {
    problem = std::make_shared<ChainedCrescent_1>(p_size);
  }
  else if (strcmp(argv[1], "ChainedCrescent_2") == 0) {
    problem = std::make_shared<ChainedCrescent_2>(p_size);
  }
  else if (strcmp(argv[1], "ChainedLQ") == 0) {
    problem = std::make_shared<ChainedLQ>(p_size);
  }
  else if (strcmp(argv[1], "ChainedMifflin_2") == 0) {
    problem = std::make_shared<ChainedMifflin_2>(p_size);
  }
  else if (strcmp(argv[1], "MaxQ") == 0) {
    problem = std::make_shared<MaxQ>(p_size);
  }
  else if (strcmp(argv[1], "MxHilb") == 0) {
    problem = std::make_shared<MxHilb>(p_size);
  }
  else if (strcmp(argv[1], "Test29_2") == 0) {
    problem = std::make_shared<Test29_2>(p_size);
  }
  else if (strcmp(argv[1], "Test29_5") == 0) {
    problem = std::make_shared<Test29_5>(p_size);
  }
  else if (strcmp(argv[1], "Test29_6") == 0) {
    problem = std::make_shared<Test29_6>(p_size);
  }
  else if (strcmp(argv[1], "Test29_11") == 0) {
    problem = std::make_shared<Test29_11>(p_size);
  }
  else if (strcmp(argv[1], "Test29_13") == 0) {
    problem = std::make_shared<Test29_13>(p_size);
  }
  else if (strcmp(argv[1], "Test29_17") == 0) {
    problem = std::make_shared<Test29_17>(p_size);
  }
  else if (strcmp(argv[1], "Test29_19") == 0) {
    problem = std::make_shared<Test29_19>(p_size);
  }
  else if (strcmp(argv[1], "Test29_20") == 0) {
    problem = std::make_shared<Test29_20>(p_size);
  }
  else if (strcmp(argv[1], "Test29_22") == 0) {
    problem = std::make_shared<Test29_22>(p_size);
  }
  else if (strcmp(argv[1], "Test29_24") == 0) {
    problem = std::make_shared<Test29_24>(p_size);
  }
  else {
    problem = std::make_shared<QuadPoly>(50, 10, 5, 10.0);
  }

  // Declare algorithm
  NonOptSolver nonopt;

  // Modify direction computation option
  nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", argv[2]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "inverse_hessian_update", argv[3]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "line_search", argv[4]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "point_set_update", argv[5]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "qp_solver", argv[6]);
  nonopt.options()->modifyStringValue(nonopt.reporter(), "symmetric_matrix", argv[7]);

//  // Declare file report
 // std::shared_ptr<FileReport> r(new FileReport("f", R_NL, R_PER_INNER_ITERATION));
//  std::shared_ptr<FileReport> r_NL(new FileReport("f", R_NL, R_PER_INNER_ITERATION));
//  std::shared_ptr<FileReport> r_QP(new FileReport("f", R_QP, R_PER_INNER_ITERATION));
//
//  // Open file
//  char out_file[100];
//  strcpy(out_file, (char*)"output/");
//  strcat(out_file, argv[1]);
//  strcat(out_file, (char*)"_");
//  strcat(out_file, argv[2]);
//  strcat(out_file, (char*)"_");
//  strcat(out_file, argv[3]);
//  strcat(out_file, (char*)"_");
//  strcat(out_file, argv[4]);
//  strcat(out_file, (char*)"_");
//  strcat(out_file, argv[5]);
//  strcat(out_file, (char*)"_");
//  strcat(out_file, argv[6]);
//  strcat(out_file, (char*)"_");
//  strcat(out_file, argv[7]);
//  strcat(out_file, (char*)".out");
//  r->open(out_file);
//
//  // Add to reporter
  //nonopt.reporter()->addReport(r);
  nonopt.reporter()->deleteReports();


  //nonopt.options()->modifyBoolValue(nonopt.reporter(), "check_derivatives", true);

  // Optimize
  nonopt.optimize(problem);


  // Return
  return 0;

}  // end main
