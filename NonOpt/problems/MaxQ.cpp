// Copyright (C) 2019 Frank E. Curtis
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
  return true;

} // end evaluateObjective

// Gradient value
bool MaxQ::evaluateGradient(int n,
                            const double* x,
                            double* g)
{

  // Initialize gradient and evaluate index of maximum of squares
  double max_val = 0.0;
  int max_ind = 0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    double temp = pow(x[i], 2.0);
    if (temp > max_val) {
      max_val = temp;
      max_ind = i;
    }
  } // end for

  // Evaluate gradient
  g[max_ind] = 2 * x[max_ind];

  // Return
  return true;

} // end evaluateGradient

// Finalize solution
bool MaxQ::finalizeSolution(int n,
                            const double* x,
                            double f,
                            const double* g)
{
  return true;
}
