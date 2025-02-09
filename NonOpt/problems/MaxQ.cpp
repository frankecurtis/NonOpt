// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "MaxQ.hpp"

// Constructor
MaxQ::MaxQ(int n)
  : number_of_variables_(n) {}

// Destructor
MaxQ::~MaxQ() {}

// Number of variables
bool MaxQ::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool MaxQ::initialPoint(int n,
                        double* x)
{

  // Set initial point
  for (int i = 0; i < n / 2; i++) {
    x[i] = i + 1;
  }
  for (int i = n / 2; i < n; i++) {
    x[i] = -i - 1;
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool MaxQ::evaluateObjective(int n,
                             const double* x,
                             double& f)
{

  // Evaluate maximum of squares
  f = 0.0;
  for (int i = 0; i < n; i++) {
    f = fmax(f, pow(x[i], 2.0));
  }

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool MaxQ::evaluateObjectiveAndGradient(int n,
                                        const double* x,
                                        double& f,
                                        double* g)
{

  // Initialize objective
  f = 0.0;

  // Initalize gradient and evaluate max of squares
  double f_temp = 0.0;
  int index = 0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    f_temp = pow(x[i], 2.0);
    if (f_temp > f) {
      f = f_temp;
      index = i;
    } // end if
  }   // end for

  // Evaluate gradient
  g[index] = 2 * x[index];

  // Return
  return !isnan(f) && !isnan(g[index]);

} // end evaluateObjectiveAndGradient

// Gradient value
bool MaxQ::evaluateGradient(int n,
                            const double* x,
                            double* g)
{

  // Initialize gradient and evaluate index of maximum of squares
  double f = 0.0;
  double f_temp = 0.0;
  int index = 0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    f_temp = pow(x[i], 2.0);
    if (f_temp > f) {
      f = f_temp;
      index = i;
    } // end if
  }   // end for

  // Evaluate gradient
  g[index] = 2 * x[index];

  // Return
  return !isnan(g[index]);

} // end evaluateGradient

// Finalize solution
bool MaxQ::finalizeSolution(int n,
                            const double* x,
                            double f,
                            const double* g)
{
  return true;
}
