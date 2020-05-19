// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedCrescent_1.hpp"

// Constructor
ChainedCrescent_1::ChainedCrescent_1(int n)
  : number_of_variables_(n) {}

// Destructor
ChainedCrescent_1::~ChainedCrescent_1() {}

// Number of variables
bool ChainedCrescent_1::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool ChainedCrescent_1::initialPoint(int n,
                                     double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    if (i % 2 == 0) {
      x[i] = -1.5;
    }
    else {
      x[i] = 2.0;
    }
  } // end for

  // Return
  return true;

} // end initialPoint

// Objective value
bool ChainedCrescent_1::evaluateObjective(int n,
                                          const double* x,
                                          double& f)
{

  // Evaluate sums
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (int i = 0; i < n - 1; i++) {
    sum1 = sum1 + pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0;
    sum2 = sum2 - pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0;
  } // end for

  // Evaluate objective
  f = fmax(sum1, sum2);

  // Return
  return true;

} // end evaluateObjective

// Gradient value
bool ChainedCrescent_1::evaluateGradient(int n,
                                         const double* x,
                                         double* g)
{

  // Initialize gradient and evaluate sums
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (int i = 0; i < n - 1; i++) {
    g[i] = 0.0;
    sum1 = sum1 + pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0;
    sum2 = sum2 - pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0;
  } // end for
  g[n - 1] = 0.0;

  // Evaluate gradient
  if (sum1 >= sum2) {
    for (int i = 0; i < n - 1; i++) {
      g[i] = g[i] + 2.0 * x[i];
      g[i + 1] = g[i + 1] + 2.0 * (x[i + 1] - 1.0) + 1.0;
    } // end for
  }   // end if
  else {
    for (int i = 0; i < n - 1; i++) {
      g[i] = g[i] - 2.0 * x[i];
      g[i + 1] = g[i + 1] - 2.0 * (x[i + 1] - 1.0) + 1.0;
    } // end for
  }   // end else

  // Return
  return true;

} // end evaluateGradient

// Finalize solution
bool ChainedCrescent_1::finalizeSolution(int n,
                                         const double* x,
                                         double f,
                                         const double* g)
{
  return true;
}
