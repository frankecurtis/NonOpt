// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "MxHilb.hpp"

// Constructor
MxHilb::MxHilb(int n)
    : number_of_variables_(n) {}

// Destructor
MxHilb::~MxHilb() {}

// Number of variables
bool MxHilb::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool MxHilb::initialPoint(int n,
                          double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = 1.0;
  }

  // Return
  return true;

}  // end initialPoint

// Objective value
bool MxHilb::evaluateObjective(int n,
                               const double* x,
                               double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum = sum + x[j] / ((double)i + (double)j + 1.0);
    }
    f = fmax(f, fabs(sum));
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool MxHilb::evaluateGradient(int n,
                              const double* x,
                              double* g)
{

  // Evaluate sums
  double* sum = new double[n];
  for (int i = 0; i < n; i++) {
    sum[i] = 0.0;
    for (int j = 0; j < n; j++) {
      sum[i] = sum[i] + x[j] / ((double)i + (double)j + 1.0);
    }
  }  // end for

  // Evaluate maximum of absolute values
  double max_val = 0.0;
  int max_ind = 0;
  for (int i = 0; i < n; i++) {
    if (fabs(sum[i]) > max_val) {
      max_val = fabs(sum[i]);
      max_ind = i;
    }
  }  // end for

  // Evaluate gradient
  if (sum[max_ind] >= 0.0) {
    for (int j = 0; j < n; j++) {
      g[j] = 1.0 / ((double)max_ind + (double)j + 1.0);
    }
  }  // end if
  else {
    for (int j = 0; j < n; j++) {
      g[j] = -1.0 / ((double)max_ind + (double)j + 1.0);
    }
  }  // end else

  // Delete sums
  delete[] sum;

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool MxHilb::finalizeSolution(int n,
                              const double* x,
                              double f,
                              const double* g)
{
  return true;
}
