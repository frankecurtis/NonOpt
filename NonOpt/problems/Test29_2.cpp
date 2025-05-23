// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_2.hpp"

// Constructor
Test29_2::Test29_2(int n)
  : number_of_variables_(n) {}

// Destructor
Test29_2::~Test29_2() {}

// Number of variables
bool Test29_2::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool Test29_2::initialPoint(int n,
                            double* x)
{

  // Set initial point
  for (int i = 0; i < n / 2; i++) {
    x[i] = double(i + 1) / (double)n;
  }
  for (int i = n / 2; i < n; i++) {
    x[i] = double(-i) / (double)n;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool Test29_2::evaluateObjective(int n,
                                 const double* x,
                                 double& f)
{

  // Evaluate maximum absolute value
  f = 0.0;
  for (int i = 0; i < n; i++) {
    f = fmax(f, fabs(x[i]));
  } // end for

  // Return
  return !std::isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool Test29_2::evaluateObjectiveAndGradient(int n,
                                            const double* x,
                                            double& f,
                                            double* g)
{

  // Initialize gradient and evaluate maximum absolute value
  int index = 0;
  f = 0.0;
  for (int i = 0; i < n; i++) {
    if (fabs(x[i]) > f) {
      index = i;
      f = fabs(x[i]);
    }
    g[i] = 0.0;
  }
  g[index] = ((x[index] > 0.0) ? 1.0 : ((x[index] < 0.0) ? -1.0 : 0.0));

  // Return
  return !std::isnan(f) && !std::isnan(g[index]);

} // end evaluateObjectiveAndGradient

// Gradient value
bool Test29_2::evaluateGradient(int n,
                                const double* x,
                                double* g)
{

  // Initialize gradient and evaluate maximum absolute value
  int index = 0;
  double f = 0.0;
  for (int i = 0; i < n; i++) {
    if (fabs(x[i]) > f) {
      index = i;
      f = fabs(x[i]);
    }
    g[i] = 0.0;
  }
  g[index] = ((x[index] > 0.0) ? 1.0 : ((x[index] < 0.0) ? -1.0 : 0.0));

  // Return
  return !std::isnan(g[index]);

} // end evaluateGradient

// Finalize solution
bool Test29_2::finalizeSolution(int n,
                                const double* x,
                                double f,
                                const double* g)
{
  return true;
}
