// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedLQ.hpp"

// Constructor
ChainedLQ::ChainedLQ(int n)
    : number_of_variables_(n) {}

// Destructor
ChainedLQ::~ChainedLQ() {}

// Number of variables
bool ChainedLQ::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool ChainedLQ::initialPoint(int n,
                             double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = -0.5;
  }

  // Return
  return true;

}  // end initialPoint

// Objective value
bool ChainedLQ::evaluateObjective(int n,
                                  const double* x,
                                  double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f = f + fmax(-x[i] - x[i + 1],
                 -x[i] - x[i + 1] + (pow(x[i], 2) + pow(x[i + 1], 2) - 1.0));
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool ChainedLQ::evaluateGradient(int n,
                                 const double* x,
                                 double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }

  // Evaluate gradient
  for (int i = 0; i < n - 1; i++) {
    if (-x[i] - x[i + 1] >= -x[i] - x[i + 1] + (pow(x[i], 2) + pow(x[i + 1], 2) - 1.0)) {
      g[i] = g[i] - 1.0;
      g[i + 1] = g[i + 1] - 1.0;
    }  // end if
    else {
      g[i] = g[i] - 1.0 + 2.0 * x[i];
      g[i + 1] = g[i + 1] - 1.0 + 2.0 * x[i + 1];
    }  // end else
  }    // end for

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool ChainedLQ::finalizeSolution(int n,
                                 const double* x,
                                 double f,
                                 const double* g)
{
  return true;
}
