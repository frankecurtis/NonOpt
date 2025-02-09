// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedCB3_1.hpp"

// Constructor
ChainedCB3_1::ChainedCB3_1(int n)
  : number_of_variables_(n) {}

// Destructor
ChainedCB3_1::~ChainedCB3_1() {}

// Number of variables
bool ChainedCB3_1::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

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

} // end initialPoint

// Objective value
bool ChainedCB3_1::evaluateObjective(int n,
                                     const double* x,
                                     double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f += fmax(pow(x[i], 4) + pow(x[i + 1], 2),
              fmax(pow(2.0 - x[i], 2) + pow(2.0 - x[i + 1], 2),
                   2.0 * exp(-x[i] + x[i + 1])));
  } // end for

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool ChainedCB3_1::evaluateObjectiveAndGradient(int n,
                                                const double* x,
                                                double& f,
                                                double* g)
{

  // Declare success
  bool success = true;

  // Evaluate objective and gradient
  f = 0.0;
  g[0] = 0.0;
  double term1, term2, term3;
  for (int i = 0; i < n - 1; i++) {
    term1 = pow(x[i], 4) + pow(x[i + 1], 2);
    term2 = pow(2.0 - x[i], 2) + pow(2.0 - x[i + 1], 2);
    term3 = 2.0 * exp(-x[i] + x[i + 1]);
    f += fmax(term1, fmax(term2, term3));
    g[i + 1] = 0.0;
    if (term1 >= term2 && term1 >= term3) {
      g[i] += 4.0 * pow(x[i], 3);
      g[i + 1] += 2.0 * x[i + 1];
    } // end if
    else if (term2 >= term1 && term2 >= term3) {
      g[i] += -2.0 * (2.0 - x[i]);
      g[i + 1] += -2.0 * (2.0 - x[i + 1]);
    } // end else if
    else {
      g[i] += -2.0 * exp(-x[i] + x[i + 1]);
      g[i + 1] += 2.0 * exp(-x[i] + x[i + 1]);
    } // end else
    if (isnan(g[i]) || isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return !isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool ChainedCB3_1::evaluateGradient(int n,
                                    const double* x,
                                    double* g)
{

  // Declare success
  bool success = true;

  // Evaluate gradient
  g[0] = 0.0;
  double term1, term2, term3;
  for (int i = 0; i < n - 1; i++) {
    g[i + 1] = 0.0;
    term1 = pow(x[i], 4) + pow(x[i + 1], 2);
    term2 = pow(2.0 - x[i], 2) + pow(2.0 - x[i + 1], 2);
    term3 = 2.0 * exp(-x[i] + x[i + 1]);
    if (term1 >= term2 && term1 >= term3) {
      g[i] += 4.0 * pow(x[i], 3);
      g[i + 1] += 2.0 * x[i + 1];
    } // end if
    else if (term2 >= term1 && term2 >= term3) {
      g[i] += -2.0 * (2.0 - x[i]);
      g[i + 1] += -2.0 * (2.0 - x[i + 1]);
    } // end else if
    else {
      g[i] += -2.0 * exp(-x[i] + x[i + 1]);
      g[i + 1] += 2.0 * exp(-x[i] + x[i + 1]);
    } // end else
    if (isnan(g[i]) || isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool ChainedCB3_1::finalizeSolution(int n,
                                    const double* x,
                                    double f,
                                    const double* g)
{
  return true;
}
