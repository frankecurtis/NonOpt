// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedMifflin_2.hpp"

// Constructor
ChainedMifflin_2::ChainedMifflin_2(int n)
  : number_of_variables_(n) {}

// Destructor
ChainedMifflin_2::~ChainedMifflin_2() {}

// Number of variables
bool ChainedMifflin_2::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool ChainedMifflin_2::initialPoint(int n,
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
bool ChainedMifflin_2::evaluateObjective(int n,
                                         const double* x,
                                         double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f = f - x[i] + 2.0 * (pow(x[i], 2) + pow(x[i + 1], 2) - 1.0) + 1.75 * fabs(pow(x[i], 2) + pow(x[i + 1], 2) - 1.0);
  }

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool ChainedMifflin_2::evaluateObjectiveAndGradient(int n,
                                                    const double* x,
                                                    double& f,
                                                    double* g)
{

  // Declare success
  bool success = true;

  // Evaluate objective and gradient
  f = 0.0;
  g[0] = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f = f - x[i] + 2.0 * (pow(x[i], 2) + pow(x[i + 1], 2) - 1.0) + 1.75 * fabs(pow(x[i], 2) + pow(x[i + 1], 2) - 1.0);
    g[i + 1] = 0.0;
    if (pow(x[i], 2) + pow(x[i + 1], 2) - 1.0 >= 0.0) {
      g[i] += -1.0 + 7.5 * x[i];
      g[i + 1] += 7.5 * x[i + 1];
    } // end if
    else {
      g[i] += -1.0 + 0.5 * x[i];
      g[i + 1] += 0.5 * x[i + 1];
    } // end else
    if (isnan(g[i]) || isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return !isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool ChainedMifflin_2::evaluateGradient(int n,
                                        const double* x,
                                        double* g)
{

  // Declare success
  bool success = true;

  // Evaluate gradient
  g[0] = 0.0;
  for (int i = 0; i < n - 1; i++) {
    g[i + 1] = 0.0;
    if (pow(x[i], 2) + pow(x[i + 1], 2) - 1.0 >= 0.0) {
      g[i] += -1.0 + 7.5 * x[i];
      g[i + 1] += 7.5 * x[i + 1];
    } // end if
    else {
      g[i] += -1.0 + 0.5 * x[i];
      g[i + 1] += 0.5 * x[i + 1];
    } // end else
    if (isnan(g[i]) || isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool ChainedMifflin_2::finalizeSolution(int n,
                                        const double* x,
                                        double f,
                                        const double* g)
{
  return true;
}
