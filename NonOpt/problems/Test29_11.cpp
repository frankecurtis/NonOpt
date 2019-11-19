// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_11.hpp"

// Constructor
Test29_11::Test29_11(int n)
    : number_of_variables_(n) {}

// Destructor
Test29_11::~Test29_11() {}

// Number of variables
bool Test29_11::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_11::initialPoint(int n,
                             double* x)
{

  // Set initial point
  for (int i = 0; i < n - 1; i++) {
    x[i] = 0.5;
  }
  x[n - 1] = -2.0;

  // Return
  return true;

}  // end initialPoint

// Objective value
bool Test29_11::evaluateObjective(int n,
                                  const double* x,
                                  double& f)
{

  // Evaluate sum of absolute values
  f = 0.0;
  for (int k = 1; k <= 2 * n - 2; k++) {
    int i = (k + 1) / 2;
    if (k % 2 == 1) {
      f += fabs(x[i - 1] + x[i] * ((5.0 - x[i]) * x[i] - 2.0) - 13.0);
    }
    else {
      f += fabs(x[i - 1] + x[i] * ((1.0 + x[i]) * x[i] - 14.0) - 29.0);
    }
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_11::evaluateGradient(int n,
                                 const double* x,
                                 double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }

  // Evaluate gradient
  for (int k = 1; k <= 2 * n - 2; k++) {
    int i = (k + 1) / 2;
    if (k % 2 == 1) {
      double term = x[i-1] + x[i] * ((5.0 - x[i]) * x[i] - 2.0) - 13.0;
      double sign = ((term >= 0.0) ? 1.0 : -1.0);
      g[i-1] += sign * (1.0);
      g[i] += sign * (-3.0 * x[i] * x[i] + 10.0 * x[i] - 2.0);
    }  // end if
    else {
      double term = x[i - 1] + x[i] * ((1.0 + x[i]) * x[i] - 14.0) - 29.0;
      double sign = ((term >= 0.0) ? 1.0 : -1.0);
      g[i-1] += sign * (1.0);
      g[i] += sign * (3.0 * x[i] * x[i] + 2.0 * x[i] - 14.0);
    }  // end else
  }    // end for

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_11::finalizeSolution(int n,
                                 const double* x,
                                 double f,
                                 const double* g)
{
  return true;
}
