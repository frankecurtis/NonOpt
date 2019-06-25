// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ActiveFaces.hpp"

// Constructor
ActiveFaces::ActiveFaces() {}

// Destructor
ActiveFaces::~ActiveFaces() {}

// Number of variables
bool ActiveFaces::numberOfVariables(int& n)
{

  // Set number of variables
  n = 50;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool ActiveFaces::initialPoint(int n,
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
bool ActiveFaces::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

  // Evaluate maximum value
  f = 0.0;
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    f = fmax(f, log(fabs(x[i]) + 1.0));
    sum = sum + x[i];
  }  // end for
  f = fmax(f, log(fabs(sum) + 1.0));

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool ActiveFaces::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
  double max_val = 0.0;
  int max_ind = 0;
  double sum = 0.0;
  double temp = 0.0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    temp = log(fabs(x[i]) + 1.0);
    if (temp > max_val) {
      max_val = temp;
      max_ind = i;
    }  // end if
    sum = sum + x[i];
  }  // end for
  temp = log(fabs(sum) + 1.0);
  if (temp > max_val) {
    max_val = temp;
    max_ind = -1;
  }  // end if

  // Evaluate gradient
  if (max_ind >= 0) {
    if (x[max_ind] >= 0) {
      g[max_ind] = 1 / (x[max_ind] + 1);
    }
    else {
      g[max_ind] = -1 / (-x[max_ind] + 1);
    }
  }  // end if
  else {
    if (-sum >= 0) {
      for (int i = 0; i < n; i++) {
        g[i] = -1 / (-sum + 1);
      }
    }
    else {
      for (int i = 0; i < n; i++) {
        g[i] = 1 / (sum + 1);
      }
    }
  }  // end else

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool ActiveFaces::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
