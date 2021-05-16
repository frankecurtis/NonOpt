// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_17.hpp"

// Constructor
Test29_17::Test29_17(int n)
  : number_of_variables_(n) {}

// Destructor
Test29_17::~Test29_17() {}

// Number of variables
bool Test29_17::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool Test29_17::initialPoint(int n,
                             double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = 1.0 / (double)n;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool Test29_17::evaluateObjective(int n,
                                  const double* x,
                                  double& f)
{

  // Evaluate maximum value
  f = 0.0;
  for (int i = 0; i < n; i++) {
    int j = i / 5;
    f = fmax(f, fabs(5 - (double)(j + 1) * (1 - cos(x[i])) - sin(x[i]) - cos(x[5 * j]) - cos(x[5 * j + 1]) - cos(x[5 * j + 2]) - cos(x[5 * j + 3]) - cos(x[5 * j + 4])));
  } // end for

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool Test29_17::evaluateObjectiveAndGradient(int n,
                                             const double* x,
                                             double& f,
                                             double* g)
{

  // Evaluate gradient
  int index = 0;
  f = 0.0;
  double maximum = 0.0;
  double term;
  int j;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    j = i / 5;
    term = 5 - (double)(j + 1) * (1 - cos(x[i])) - sin(x[i]) - cos(x[5 * j]) - cos(x[5 * j + 1]) - cos(x[5 * j + 2]) - cos(x[5 * j + 3]) - cos(x[5 * j + 4]);
    if (fabs(term) > f) {
      index = i;
      f = fabs(term);
      maximum = term;
    } // end if
  }   // end for
  j = index / 5;
  double sign = ((maximum >= 0) ? 1.0 : -1.0);
  g[index] += sign * (-(double)(j + 1) * sin(x[index]) - cos(x[index]));
  g[5 * j] += sign * sin(x[5 * j]);
  g[5 * j + 1] += sign * sin(x[5 * j + 1]);
  g[5 * j + 2] += sign * sin(x[5 * j + 2]);
  g[5 * j + 3] += sign * sin(x[5 * j + 3]);
  g[5 * j + 4] += sign * sin(x[5 * j + 4]);

  // Return
  return !isnan(g[index]) && !isnan(g[5 * j]) && !isnan(g[5 * j + 1]) && !isnan(g[5 * j + 2]) && !isnan(g[5 * j + 3]) && !isnan(g[5 * j + 4]);

} // end evaluateObjectiveAndGradient

// Gradient value
bool Test29_17::evaluateGradient(int n,
                                 const double* x,
                                 double* g)
{

  // Evaluate gradient
  int index = 0;
  double f = 0.0;
  double maximum = 0.0;
  double term;
  int j;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    j = i / 5;
    term = 5 - (double)(j + 1) * (1 - cos(x[i])) - sin(x[i]) - cos(x[5 * j]) - cos(x[5 * j + 1]) - cos(x[5 * j + 2]) - cos(x[5 * j + 3]) - cos(x[5 * j + 4]);
    if (fabs(term) > f) {
      index = i;
      f = fabs(term);
      maximum = term;
    } // end if
  }   // end for
  j = index / 5;
  double sign = ((maximum >= 0) ? 1.0 : -1.0);
  g[index] += sign * (-(double)(j + 1) * sin(x[index]) - cos(x[index]));
  g[5 * j] += sign * sin(x[5 * j]);
  g[5 * j + 1] += sign * sin(x[5 * j + 1]);
  g[5 * j + 2] += sign * sin(x[5 * j + 2]);
  g[5 * j + 3] += sign * sin(x[5 * j + 3]);
  g[5 * j + 4] += sign * sin(x[5 * j + 4]);

  // Return
  return !isnan(g[index]) && !isnan(g[5 * j]) && !isnan(g[5 * j + 1]) && !isnan(g[5 * j + 2]) && !isnan(g[5 * j + 3]) && !isnan(g[5 * j + 4]);

} // end evaluateGradient

// Finalize solution
bool Test29_17::finalizeSolution(int n,
                                 const double* x,
                                 double f,
                                 const double* g)
{
  return true;
}
