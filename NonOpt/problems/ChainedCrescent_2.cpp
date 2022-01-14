// Copyright (C) 2022 Frank E. Curtis
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

} // end numberOfVariables

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
  } // end for

  // Return
  return true;

} // end initialPoint

// Objective value
bool ChainedCrescent_2::evaluateObjective(int n,
                                          const double* x,
                                          double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f += fmax(pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0,
              -pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0);
  } // end for

  // Return
  return !std::isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool ChainedCrescent_2::evaluateObjectiveAndGradient(int n,
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
    f += fmax(pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0,
              -pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0);
    g[i + 1] = 0.0;
    if (pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0 >=
        -pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0) {
      g[i] += 2.0 * x[i];
      g[i + 1] += 2.0 * (x[i + 1] - 1.0) + 1.0;
    } // end if
    else {
      g[i] += -2.0 * x[i];
      g[i + 1] += -2.0 * (x[i + 1] - 1.0) + 1.0;
    } // end else
    if (std::isnan(g[i]) || std::isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return !std::isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool ChainedCrescent_2::evaluateGradient(int n,
                                         const double* x,
                                         double* g)
{

  // Declare success
  bool success = true;

  // Evaluate gradient
  g[0] = 0.0;
  for (int i = 0; i < n - 1; i++) {
    g[i + 1] = 0.0;
    if (pow(x[i], 2) + pow(x[i + 1] - 1.0, 2) + x[i + 1] - 1.0 >=
        -pow(x[i], 2) - pow(x[i + 1] - 1.0, 2) + x[i + 1] + 1.0) {
      g[i] += 2.0 * x[i];
      g[i + 1] += 2.0 * (x[i + 1] - 1.0) + 1.0;
    } // end if
    else {
      g[i] += -2.0 * x[i];
      g[i + 1] += -2.0 * (x[i + 1] - 1.0) + 1.0;
    } // end else
    if (std::isnan(g[i]) || std::isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool ChainedCrescent_2::finalizeSolution(int n,
                                         const double* x,
                                         double f,
                                         const double* g)
{
  return true;
}
