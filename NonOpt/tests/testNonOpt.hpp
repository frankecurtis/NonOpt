// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTNONOPT_HPP__
#define __TESTNONOPT_HPP__

#include "NonOptProblem.hpp"
#include "NonOptReporter.hpp"
#include "NonOptSolver.hpp"

#include "MaxQ.hpp"

using namespace NonOpt;

// Implementation of test
int testNonOptImplementation()
{

  // Declare NonOptSolver
  NonOptSolver nonopt;

  // Declare problem pointer
  std::shared_ptr<Problem> problem = std::make_shared<MaxQ>(10);

  // Print message
  Reporter r;
  std::shared_ptr<StreamReport> s(new StreamReport("temp", R_NL, R_BASIC));
  s->setStream(&std::cout);
  r.addReport(s);
  r.printf(R_NL, R_BASIC, "solving MaxQ(10):\n");

  // Optimize
  nonopt.optimize(problem);

  // Return
  return 0;

} // end testNonOptImplementation

#endif /* __TESTNONOPT_HPP__ */
