// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "BrownFunction_2.hpp"

// Constructor
BrownFunction_2::BrownFunction_2(int n)
  : number_of_variables_(n) {}

// Destructor
BrownFunction_2::~BrownFunction_2() {}

// Number of variables
bool BrownFunction_2::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool BrownFunction_2::initialPoint(int n,
                                   double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = pow(-1.0, (i + 1));
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool BrownFunction_2::evaluateObjective(int n,
                                        const double* x,
                                        double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n - 1; i++) {
    f = f + pow(fabs(x[i]), (x[i + 1] * x[i + 1] + 1.0)) + pow(fabs(x[i + 1]), (x[i] * x[i] + 1.0));
  }

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool BrownFunction_2::evaluateObjectiveAndGradient(int n,
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
    f = f + pow(fabs(x[i]), (x[i + 1] * x[i + 1] + 1.0)) + pow(fabs(x[i + 1]), (x[i] * x[i] + 1.0));
    g[i + 1] = 0.0;
    if (x[i] >= 0.0 && x[i + 1] >= 0.0) {
      g[i] += (pow(x[i + 1], 2) + 1.0) * pow(x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(x[i + 1]) * pow(x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(x[i]) * pow(x[i], pow(x[i + 1], 2) + 1.0) + (pow(x[i], 2) + 1.0) * pow(x[i + 1], pow(x[i], 2));
    } // end if
    else if (x[i] >= 0.0 && x[i + 1] < 0.0) {
      g[i] += (pow(x[i + 1], 2) + 1.0) * pow(x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(-x[i + 1]) * pow(-x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(x[i]) * pow(x[i], pow(x[i + 1], 2) + 1.0) - (pow(x[i], 2) + 1.0) * pow(-x[i + 1], pow(x[i], 2));
    } // end else if
    else if (x[i] < 0.0 && x[i + 1] >= 0.0) {
      g[i] += -(pow(x[i + 1], 2) + 1.0) * pow(-x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(x[i + 1]) * pow(x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(-x[i]) * pow(-x[i], pow(x[i + 1], 2) + 1.0) + (pow(x[i], 2) + 1.0) * pow(x[i + 1], pow(x[i], 2));
    } // end else if
    else {
      g[i] += -(pow(x[i + 1], 2) + 1.0) * pow(-x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(-x[i + 1]) * pow(-x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(-x[i]) * pow(-x[i], pow(x[i + 1], 2) + 1.0) - (pow(x[i], 2) + 1.0) * pow(-x[i + 1], pow(x[i], 2));
    } // end else
    if (isnan(g[i]) || isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return !isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool BrownFunction_2::evaluateGradient(int n,
                                       const double* x,
                                       double* g)
{

  // Declare success
  bool success = true;

  // Evaluate gradient
  g[0] = 0.0;
  for (int i = 0; i < n - 1; i++) {
    g[i + 1] = 0.0;
    if (x[i] >= 0.0 && x[i + 1] >= 0.0) {
      g[i] += (pow(x[i + 1], 2) + 1.0) * pow(x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(x[i + 1]) * pow(x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(x[i]) * pow(x[i], pow(x[i + 1], 2) + 1.0) + (pow(x[i], 2) + 1.0) * pow(x[i + 1], pow(x[i], 2));
    } // end if
    else if (x[i] >= 0.0 && x[i + 1] < 0.0) {
      g[i] += (pow(x[i + 1], 2) + 1.0) * pow(x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(-x[i + 1]) * pow(-x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(x[i]) * pow(x[i], pow(x[i + 1], 2) + 1.0) - (pow(x[i], 2) + 1.0) * pow(-x[i + 1], pow(x[i], 2));
    } // end else if
    else if (x[i] < 0.0 && x[i + 1] >= 0.0) {
      g[i] += -(pow(x[i + 1], 2) + 1.0) * pow(-x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(x[i + 1]) * pow(x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(-x[i]) * pow(-x[i], pow(x[i + 1], 2) + 1.0) + (pow(x[i], 2) + 1.0) * pow(x[i + 1], pow(x[i], 2));
    } // end else if
    else {
      g[i] += -(pow(x[i + 1], 2) + 1.0) * pow(-x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(-x[i + 1]) * pow(-x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] += 2.0 * x[i + 1] * log(-x[i]) * pow(-x[i], pow(x[i + 1], 2) + 1.0) - (pow(x[i], 2) + 1.0) * pow(-x[i + 1], pow(x[i], 2));
    } // end else
    if (isnan(g[i]) || isnan(g[i + 1])) {
      success = false;
    }
  } // end for

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool BrownFunction_2::finalizeSolution(int n,
                                       const double* x,
                                       double f,
                                       const double* g)
{
  return true;
}
