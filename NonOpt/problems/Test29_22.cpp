// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_22.hpp"

// Constructor
Test29_22::Test29_22(int n)
  : number_of_variables_(n) {}

// Destructor
Test29_22::~Test29_22() {}

// Number of variables
bool Test29_22::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool Test29_22::initialPoint(int n,
                             double* x)
{

  // Set initial point
  for (int i = 0; i < n - 1; i++) {
    x[i] = ((double)(i + 1) / (double)n) * (((double)(i + 1) / (double)n) - 1.0);
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool Test29_22::evaluateObjective(int n,
                                  const double* x,
                                  double& f)
{

  // Evaluate maximum value
  f = fabs(2 * x[0] + (1.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[0] + (double)(1) / ((double)(n + 1)) + 1.0, 3.0) - x[1]);
  for (int i = 1; i < n - 1; i++) {
    f = fmax(f, fabs(2 * x[i] + (1.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[i] + (double)(i + 1) / ((double)(n + 1)) + 1.0, 3.0) - x[i - 1] - x[i + 1]));
  }
  f = fmax(f, fabs(2 * x[n - 1] + (1.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[n - 1] + (double)n / ((double)(n + 1)) + 1.0, 3.0) - x[n - 2]));

  // Return
  return true;

} // end evaluateObjective

// Gradient value
bool Test29_22::evaluateGradient(int n,
                                 const double* x,
                                 double* g)
{

  // Initialize gradient and evaluate maximum value
  int max_ind = 0;
  double term = 2 * x[0] + (1.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[0] + (double)(1) / ((double)(n + 1)) + 1.0, 3.0) - x[1];
  double max_term = term;
  double max_val = fabs(term);
  g[0] = 0.0;
  for (int i = 1; i < n - 1; i++) {
    term = 2 * x[i] + (1.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[i] + (double)(i + 1) / ((double)(n + 1)) + 1.0, 3.0) - x[i - 1] - x[i + 1];
    if (fabs(term) > max_val) {
      max_ind = i;
      max_term = term;
      max_val = fabs(term);
    } // end if
    g[i] = 0.0;
  } // end for
  term = 2 * x[n - 1] + (1.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[n - 1] + (double)n / ((double)(n + 1)) + 1.0, 3.0) - x[n - 2];
  if (fabs(term) > max_val) {
    max_ind = n - 1;
    max_term = term;
    max_val = fabs(term);
  } // end if
  g[n - 1] = 0.0;

  // Evaluate gradient
  double sign = ((max_term >= 0.0) ? 1.0 : -1.0);
  g[max_ind] = sign * (2.0 + (3.0 / (2.0 * (double)(n * n + 2 * n + 1))) * pow(x[max_ind] + (double)(max_ind + 1) / ((double)(n + 1)) + 1.0, 2.0));
  if (max_ind > 0) {
    g[max_ind - 1] = sign * (-1.0);
  }
  if (max_ind < n - 1) {
    g[max_ind + 1] = sign * (-1.0);
  }

  // Return
  return true;

} // end evaluateGradient

// Finalize solution
bool Test29_22::finalizeSolution(int n,
                                 const double* x,
                                 double f,
                                 const double* g)
{
  return true;
}
