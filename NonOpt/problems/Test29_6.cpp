// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_6.hpp"

// Constructor
Test29_6::Test29_6(int n)
  : number_of_variables_(n) {}

// Destructor
Test29_6::~Test29_6() {}

// Number of variables
bool Test29_6::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool Test29_6::initialPoint(int n,
                            double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = -1.0;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool Test29_6::evaluateObjective(int n,
                                 const double* x,
                                 double& f)
{

  // Evaluate maximum value
  f = fabs((3.0 - 2.0 * x[0]) * x[0] + 1.0 - x[1]);
  for (int i = 1; i < n - 1; i++) {
    f = fmax(f, fabs((3.0 - 2.0 * x[i]) * x[i] + 1.0 - x[i - 1] - x[i + 1]));
  } // end for
  f = fmax(f, fabs((3.0 - 2.0 * x[n - 1]) * x[n - 1] + 1.0 - x[n - 2]));

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool Test29_6::evaluateObjectiveAndGradient(int n,
                                            const double* x,
                                            double& f,
                                            double* g)
{

  // Declare success
  bool success = true;

  // Initialize gradient and evaluate maximum value
  int index = 0;
  double term = (3.0 - 2.0 * x[0]) * x[0] + 1.0 - x[1];
  double maximum = term;
  f = fabs(term);
  g[0] = 0.0;
  for (int i = 1; i < n - 1; i++) {
    term = (3.0 - 2.0 * x[i]) * x[i] + 1.0 - x[i - 1] - x[i + 1];
    if (fabs(term) > f) {
      index = i;
      maximum = term;
      f = fabs(term);
    } // end if
    g[i] = 0.0;
  } // end for
  term = (3.0 - 2.0 * x[n - 1]) * x[n - 1] + 1.0 - x[n - 2];
  if (fabs(term) > f) {
    index = n - 1;
    maximum = term;
    f = fabs(term);
  } // end if
  g[n - 1] = 0.0;

  // Evaluate gradient
  double sign = ((maximum >= 0.0) ? 1.0 : -1.0);
  g[index] = sign * (3.0 - 4.0 * x[index]);
  if (isnan(g[index])) {
    success = false;
  }
  if (index > 0) {
    g[index - 1] = sign * (-1.0);
    if (isnan(g[index - 1])) {
      success = false;
    }
  }
  if (index < n - 1) {
    g[index + 1] = sign * (-1.0);
    if (isnan(g[index + 1])) {
      success = false;
    }
  }

  // Return
  return !isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool Test29_6::evaluateGradient(int n,
                                const double* x,
                                double* g)
{

  // Declare success
  bool success = true;

  // Initialize gradient and evaluate maximum value
  int index = 0;
  double term = (3.0 - 2.0 * x[0]) * x[0] + 1.0 - x[1];
  double maximum = term;
  double f = fabs(term);
  g[0] = 0.0;
  for (int i = 1; i < n - 1; i++) {
    term = (3.0 - 2.0 * x[i]) * x[i] + 1.0 - x[i - 1] - x[i + 1];
    if (fabs(term) > f) {
      index = i;
      maximum = term;
      f = fabs(term);
    } // end if
    g[i] = 0.0;
  } // end for
  term = (3.0 - 2.0 * x[n - 1]) * x[n - 1] + 1.0 - x[n - 2];
  if (fabs(term) > f) {
    index = n - 1;
    maximum = term;
    f = fabs(term);
  } // end if
  g[n - 1] = 0.0;

  // Evaluate gradient
  double sign = ((maximum >= 0.0) ? 1.0 : -1.0);
  g[index] = sign * (3.0 - 4.0 * x[index]);
  if (isnan(g[index])) {
    success = false;
  }
  if (index > 0) {
    g[index - 1] = sign * (-1.0);
    if (isnan(g[index - 1])) {
      success = false;
    }
  }
  if (index < n - 1) {
    g[index + 1] = sign * (-1.0);
    if (isnan(g[index + 1])) {
      success = false;
    }
  }

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool Test29_6::finalizeSolution(int n,
                                const double* x,
                                double f,
                                const double* g)
{
  return true;
}
