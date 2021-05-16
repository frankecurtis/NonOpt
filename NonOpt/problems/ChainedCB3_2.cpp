// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ChainedCB3_2.hpp"

// Constructor
ChainedCB3_2::ChainedCB3_2(int n)
  : number_of_variables_(n) {}

// Destructor
ChainedCB3_2::~ChainedCB3_2() {}

// Number of variables
bool ChainedCB3_2::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool ChainedCB3_2::initialPoint(int n,
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
bool ChainedCB3_2::evaluateObjective(int n,
                                     const double* x,
                                     double& f)
{

  // Evaluate sums
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  for (int i = 0; i < n - 1; i++) {
    sum1 += pow(x[i], 4) + pow(x[i + 1], 2);
    sum2 += pow(2 - x[i], 2) + pow(2 - x[i + 1], 2);
    sum3 += 2.0 * exp(-x[i] + x[i + 1]);
  } // end for

  // Evaluate objective
  f = fmax(sum1, fmax(sum2, sum3));

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool ChainedCB3_2::evaluateObjectiveAndGradient(int n,
                                                const double* x,
                                                double& f,
                                                double* g)
{

  // Declare success
  bool success = true;

  // Initialize gradient and evaluate sums
  g[0] = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  for (int i = 0; i < n - 1; i++) {
    g[i+1] = 0.0;
    sum1 += pow(x[i], 4) + pow(x[i + 1], 2);
    sum2 += pow(2 - x[i], 2) + pow(2 - x[i + 1], 2);
    sum3 += 2.0 * exp(-x[i] + x[i + 1]);
  } // end for

  // Evaluate objective
  f = fmax(sum1, fmax(sum2, sum3));

  // Determine index of maximum
  int index = 1;
  if (sum1 >= sum2 && sum1 >= sum3) {
    index = 1;
  }
  else if (sum2 >= sum1 && sum2 >= sum3) {
    index = 2;
  }
  else {
    index = 3;
  }

  // Evaluate gradient
  if (index == 1) {
    for (int i = 0; i < n - 1; i++) {
      g[i] += 4.0 * pow(x[i], 3);
      g[i + 1] += 2.0 * x[i + 1];
      if (isnan(g[i]) || isnan(g[i+1])) {
        success = false;
      }
    } // end for
  }   // end if
  else if (index == 2) {
    for (int i = 0; i < n - 1; i++) {
      g[i] += -2.0 * (2.0 - x[i]);
      g[i + 1] += -2.0 * (2.0 - x[i + 1]);
      if (isnan(g[i]) || isnan(g[i+1])) {
        success = false;
      }
    } // end for
  }   // end else if
  else {
    for (int i = 0; i < n - 1; i++) {
      g[i] += -2.0 * exp(-x[i] + x[i + 1]);
      g[i + 1] += 2.0 * exp(-x[i] + x[i + 1]);
      if (isnan(g[i]) || isnan(g[i+1])) {
        success = false;
      }
    } // end for
  }   // end else

  // Return
  return !isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool ChainedCB3_2::evaluateGradient(int n,
                                    const double* x,
                                    double* g)
{

  // Declare success
  bool success = true;

  // Initialize gradient and evaluate sums
  g[0] = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  for (int i = 0; i < n - 1; i++) {
    g[i+1] = 0.0;
    sum1 += pow(x[i], 4) + pow(x[i + 1], 2);
    sum2 += pow(2 - x[i], 2) + pow(2 - x[i + 1], 2);
    sum3 += 2.0 * exp(-x[i] + x[i + 1]);
  } // end for

  // Determine index of maximum
  int index = 1;
  if (sum1 >= sum2 && sum1 >= sum3) {
    index = 1;
  }
  else if (sum2 >= sum1 && sum2 >= sum3) {
    index = 2;
  }
  else {
    index = 3;
  }

  // Evaluate gradient
  if (index == 1) {
    for (int i = 0; i < n - 1; i++) {
      g[i] += 4.0 * pow(x[i], 3);
      g[i + 1] += 2.0 * x[i + 1];
      if (isnan(g[i]) || isnan(g[i+1])) {
        success = false;
      }
    } // end for
  }   // end if
  else if (index == 2) {
    for (int i = 0; i < n - 1; i++) {
      g[i] += -2.0 * (2.0 - x[i]);
      g[i + 1] += -2.0 * (2.0 - x[i + 1]);
      if (isnan(g[i]) || isnan(g[i+1])) {
        success = false;
      }
    } // end for
  }   // end else if
  else {
    for (int i = 0; i < n - 1; i++) {
      g[i] += -2.0 * exp(-x[i] + x[i + 1]);
      g[i + 1] += 2.0 * exp(-x[i] + x[i + 1]);
      if (isnan(g[i]) || isnan(g[i+1])) {
        success = false;
      }
    } // end for
  }   // end else

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool ChainedCB3_2::finalizeSolution(int n,
                                    const double* x,
                                    double f,
                                    const double* g)
{
  return true;
}
