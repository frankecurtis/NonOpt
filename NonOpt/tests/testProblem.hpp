// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTPROBLEM_HPP__
#define __TESTPROBLEM_HPP__

#include <iostream>

#include "MaxQ.hpp"
#include "NonOptReporter.hpp"

using namespace NonOpt;

// Main function
int testProblemImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter reporter;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> s(new StreamReport("s", R_NL, R_BASIC));

    // Set stream report to standard output
    s->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(s);

  } // end if

  // Declare problem
  MaxQ problem(50);

  // Declare problem size
  int n;

  // Get number of variables
  problem.numberOfVariables(n);

  // Check length
  if (n != 50) {
    result = 1;
  }

  // Print number of variables
  reporter.printf(R_NL, R_BASIC, "Number of variables (should be 50): %d\n", n);

  // Declare quantities
  double* x = new double[n];
  double f;
  double* g = new double[n];

  // Get initial point
  problem.initialPoint(n, x);

  // Print initial point and check values
  reporter.printf(R_NL, R_BASIC, "Initial point... should be [1, ..., 25, -26, ..., -50]:\n");
  for (int i = 0; i < n; i++) {
    if (i < 25) {
      if (x[i] < i + 1 - 1e-12 || x[i] > i + 1 + 1e-12) {
        result = 1;
      }
    }
    else {
      if (x[i] < -i - 1 - 1e-12 || x[i] > -i - 1 + 1e-12) {
        result = 1;
      }
    }
    reporter.printf(R_NL, R_BASIC, "%+23.16e\n", x[i]);
  } // end for

  // Evaluate objective value
  problem.evaluateObjective(n, x, f);

  // Check objective value
  if (f < 2500.0 - 1e-12 || f > 2500.0 + 1e-12) {
    result = 1;
  }

  // Print objective value
  reporter.printf(R_NL, R_BASIC, "Objective value... should be 2500: %+23.16e\n", f);

  // Evaluate gradient value
  problem.evaluateGradient(n, x, g);

  // Print gradient and check values
  reporter.printf(R_NL, R_BASIC, "Gradient... should be [0, ..., 0, -100]:\n");
  for (int i = 0; i < n; i++) {
    if (i < n - 1) {
      if (g[i] < 0.0 - 1e-12 || g[i] > 0.0 + 1e-12) {
        result = 1;
      }
    }
    else {
      if (g[i] < -100.0 - 1e-12 || g[i] > -100.0 + 1e-12) {
        result = 1;
      }
    }
    reporter.printf(R_NL, R_BASIC, "%+23.16e\n", g[i]);
  } // end for

  // Delete objects
  delete[] x;
  delete[] g;

  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      reporter.printf(R_NL, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    }
    else {
      reporter.printf(R_NL, R_BASIC, "TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testProblemImplementation

#endif /* __TESTPROBLEM_HPP__ */
