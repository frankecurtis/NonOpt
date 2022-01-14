// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cstdio>

#include "testNonOpt.hpp"
#include "testOptions.hpp"
#include "testPoint.hpp"
#include "testProblem.hpp"
#include "testQPSolver.hpp"
#include "testQuantities.hpp"
#include "testReporter.hpp"
#include "testSymmetricMatrix.hpp"
#include "testVector.hpp"

// Main function
int main()
{

  // Initialize result
  int result = 0;

  // Run tests
  printf("testing Quantities................. ");
  if (!testQuantitiesImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testQuantities for details)\n");
  }
  printf("testing Options.................... ");
  if (!testOptionsImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testOptions for details)\n");
  }
  printf("testing Point...................... ");
  if (!testPointImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testPoint for details)\n");
  }
  printf("testing Problem.................... ");
  if (!testProblemImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testProblem for details)\n");
  }
  printf("testing QPSolver................... ");
  if (!testQPSolverImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testQPSolver for details)\n");
  }
  printf("testing Reporter................... ");
  if (!testReporterImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testReporter for details)\n");
  }
  printf("testing SymmetricMatrix............ ");
  if (!testSymmetricMatrixImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testSymmetricMatrix for details)\n");
  }
  printf("testing Vector..................... ");
  if (!testVectorImplementation(0)) {
    printf("success.\n");
  }
  else {
    result = 1;
    printf("failure! (run testVector for details)\n");
  }

  // Test algorithm if results are good
  if (result == 0) {
    printf("testing NonOpt..................... ");
    testNonOptImplementation();
  }
  else {
    printf("testing NonOpt..................... skipped (due to failure above).\n");
  }

  // Return
  return 0;

} // end main
