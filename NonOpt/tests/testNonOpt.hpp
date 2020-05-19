// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTNONOPT_HPP__
#define __TESTNONOPT_HPP__

#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"

#include "MaxQ.hpp"

using namespace NonOpt;

// Implementation of test
int testNonOptImplementation()
{

  // Declare NonOptSolver
  NonOptSolver nonopt;

  // Change print level
  std::shared_ptr<Report> report = nonopt.reporter()->report("default");
  report->setTypeAndLevel(R_NL, R_BASIC);

  // Declare problem pointer
  std::shared_ptr<Problem> problem = std::make_shared<MaxQ>(10);

  // Print message
  nonopt.reporter()->printf(R_NL, R_BASIC, "solving MaxQ(10):\n");

  // Optimize
  nonopt.optimize(problem);

  // Return
  return 0;

} // end testNonOptImplementation

#endif /* __TESTNONOPT_HPP__ */
