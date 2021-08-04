// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>
#include <cstdio>

#include "Test29_24.hpp"

// Constructor
Test29_24::Test29_24(int n)
  : number_of_variables_(n) {}

// Destructor
Test29_24::~Test29_24() {}

// Number of variables
bool Test29_24::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool Test29_24::initialPoint(int n,
                             double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = 1.0;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool Test29_24::evaluateObjective(int n,
                                  const double* x,
                                  double& f)
{

  // Evaluate maximum value
  f = fabs(2 * x[0] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[0]) - x[1]);
  for (int i = 1; i < n - 1; i++) {
    f = fmax(f, fabs(2 * x[i] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[i]) - x[i - 1] - x[i + 1]));
  }
  f = fmax(f, fabs(2 * x[n - 1] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[n - 1]) - x[n - 2] - 1.0));

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool Test29_24::evaluateObjectiveAndGradient(int n,
                                             const double* x,
                                             double& f,
                                             double* g)
{

  // Declare success
  bool success = true;

  // Initialize gradient and evaluate maximum value
  int index = 0;
  double term = 2 * x[0] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[0]) - x[1];
  double maximum = term;
  f = fabs(term);
  g[0] = 0.0;
  for (int i = 1; i < n - 1; i++) {
    term = 2 * x[i] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[i]) - x[i - 1] - x[i + 1];
    if (fabs(term) > f) {
      index = i;
      maximum = term;
      f = fabs(term);
    } // end if
    g[i] = 0.0;
  } // end for
  term = 2 * x[n - 1] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[n - 1]) - x[n - 2] - 1.0;
  if (fabs(term) > f) {
    index = n - 1;
    maximum = term;
    f = fabs(term);
  } // end if
  g[n - 1] = 0.0;

  // Evaluate gradient
  double sign = ((maximum >= 0.0) ? 1.0 : -1.0);
  g[index] = sign * (2.0 + (100.0 / ((double)(n * n + 2 * n + 1))) * cosh(10.0 * x[index]));
  if (isnan(g[index])) {
    success = false;
  }
  if (index > 0) {
    g[index - 1] = sign * (-1.0);
    if (isnan(g[index - 1])) {
      success = false;
    }
  } // end if
  if (index < n - 1) {
    g[index + 1] = sign * (-1.0);
    if (isnan(g[index + 1])) {
      success = false;
    }
  } // end if

  // Return
  return !isnan(f) && success;

} // end evaluateObjective

// Gradient value
bool Test29_24::evaluateGradient(int n,
                                 const double* x,
                                 double* g)
{

  // Declare success
  bool success = true;

  // Initialize gradient and evaluate maximum value
  int index = 0;
  double term = 2 * x[0] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[0]) - x[1];
  double maximum = term;
  double f = fabs(term);
  g[0] = 0.0;
  for (int i = 1; i < n - 1; i++) {
    term = 2 * x[i] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[i]) - x[i - 1] - x[i + 1];
    if (fabs(term) > f) {
      index = i;
      maximum = term;
      f = fabs(term);
    } // end if
    g[i] = 0.0;
  } // end for
  term = 2 * x[n - 1] + (10.0 / ((double)(n * n + 2 * n + 1))) * sinh(10.0 * x[n - 1]) - x[n - 2] - 1.0;
  if (fabs(term) > f) {
    index = n - 1;
    maximum = term;
    f = fabs(term);
  } // end if
  g[n - 1] = 0.0;

  // Evaluate gradient
  double sign = ((maximum >= 0.0) ? 1.0 : -1.0);
  g[index] = sign * (2.0 + (100.0 / ((double)(n * n + 2 * n + 1))) * cosh(10.0 * x[index]));
  if (isnan(g[index])) {
    success = false;
  }
  if (index > 0) {
    g[index - 1] = sign * (-1.0);
    if (isnan(g[index - 1])) {
      success = false;
    }
  } // end if
  if (index < n - 1) {
    g[index + 1] = sign * (-1.0);
    if (isnan(g[index + 1])) {
      success = false;
    }
  } // end if

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool Test29_24::finalizeSolution(int n,
                                 const double* x,
                                 double f,
                                 const double* g)
{
  return true;
}
