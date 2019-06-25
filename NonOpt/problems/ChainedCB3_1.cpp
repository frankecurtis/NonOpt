// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedCB3_1.hpp"

// Constructor
ChainedCB3_1::ChainedCB3_1() {}

// Destructor
ChainedCB3_1::~ChainedCB3_1() {}

// Number of variables
bool ChainedCB3_1::numberOfVariables(int& n)
{

  // Set number of variables
  n = 50;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool ChainedCB3_1::initialPoint(int n,
                                double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = 2.0;
  }

  // Return
  return true;

}  // end initialPoint

// Objective value
bool ChainedCB3_1::evaluateObjective(int n,
                                     const double* x,
                                     double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f = f + fmax(pow(x[i], 4) + pow(x[i + 1], 2),
                 fmax(pow(2.0 - x[i], 2) + pow(2.0 - x[i + 1], 2),
                      2.0 * exp(-x[i] + x[i + 1])));
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool ChainedCB3_1::evaluateGradient(int n,
                                    const double* x,
                                    double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }

  // Evaluate gradient
  for (int i = 0; i < n - 1; i++) {
    double term1 = pow(x[i], 4) + pow(x[i + 1], 2);
    double term2 = pow(2.0 - x[i], 2) + pow(2.0 - x[i + 1], 2);
    double term3 = 2.0 * exp(-x[i] + x[i + 1]);
    if (term1 >= term2 && term1 >= term3) {
      g[i] = g[i] + 4.0 * pow(x[i], 3);
      g[i + 1] = g[i + 1] + 2.0 * x[i + 1];
    }  // end if
    else if (term2 >= term1 && term2 >= term3) {
      g[i] = g[i] - 2.0 * (2.0 - x[i]);
      g[i + 1] = g[i + 1] - 2.0 * (2.0 - x[i + 1]);
    }  // end else if
    else {
      g[i] = g[i] - 2.0 * exp(-x[i] + x[i + 1]);
      g[i + 1] = g[i + 1] + 2.0 * exp(-x[i] + x[i + 1]);
    }  // end else
  }    // end for

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool ChainedCB3_1::finalizeSolution(int n,
                                    const double* x,
                                    double f,
                                    const double* g)
{
  return true;
}
