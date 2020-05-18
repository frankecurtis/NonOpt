// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedCrescent_2.hpp"

// Constructor
ChainedCrescent_2::ChainedCrescent_2(int n)
    : number_of_variables_(n) {}

// Destructor
ChainedCrescent_2::~ChainedCrescent_2() {}

// Number of variables
bool ChainedCrescent_2::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool ChainedCrescent_2::initialPoint(int n,
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
  }  // end for

  // Return
  return true;

}  // end initialPoint

// Objective value
bool ChainedCrescent_2::evaluateObjective(int n,
                                          const double* x,
                                          double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f = f + fmax(pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0,
                 -pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0);
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool ChainedCrescent_2::evaluateGradient(int n,
                                         const double* x,
                                         double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }

  // Evaluate gradient
  for (int i = 0; i < n - 1; i++) {
    if (pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0 >=
        -pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0) {
      g[i] = g[i] + 2.0 * x[i];
      g[i + 1] = g[i + 1] + 2.0 * (x[i + 1] - 1.0) + 1.0;
    }  // end if
    else {
      g[i] = g[i] - 2.0 * x[i];
      g[i + 1] = g[i + 1] - 2.0 * (x[i + 1] - 1.0) + 1.0;
    }  // end else
  }    // end for

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool ChainedCrescent_2::finalizeSolution(int n,
                                         const double* x,
                                         double f,
                                         const double* g)
{
  return true;
}
