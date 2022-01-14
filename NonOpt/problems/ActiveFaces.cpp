// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ActiveFaces.hpp"

// Constructor
ActiveFaces::ActiveFaces(int n)
  : number_of_variables_(n) {}

// Destructor
ActiveFaces::~ActiveFaces() {}

// Number of variables
bool ActiveFaces::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

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

} // end initialPoint

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
  } // end for
  f = fmax(f, log(fabs(sum) + 1.0));

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool ActiveFaces::evaluateObjectiveAndGradient(int n,
                                               const double* x,
                                               double& f,
                                               double* g)
{

  // Initialize gradient and evaluate maximum value
  f = 0.0;
  int index = 0;
  double sum = 0.0;
  double temp = 0.0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    temp = log(fabs(x[i]) + 1.0);
    if (temp > f) {
      f = temp;
      index = i;
    } // end if
    sum = sum + x[i];
  } // end for
  temp = log(fabs(sum) + 1.0);
  if (temp > f) {
    f = temp;
    index = -1;
  } // end if

  // Declare success
  bool success = !isnan(f);

  // Evaluate gradient
  if (index >= 0) {
    if (x[index] >= 0) {
      g[index] = 1 / (x[index] + 1);
    }
    else {
      g[index] = -1 / (-x[index] + 1);
    }
    success = !isnan(g[index]);
  } // end if
  else {
    if (-sum >= 0) {
      for (int i = 0; i < n; i++) {
        g[i] = -1 / (-sum + 1);
      }
      success = !isnan(-1 / (-sum + 1));
    }
    else {
      for (int i = 0; i < n; i++) {
        g[i] = 1 / (sum + 1);
      }
      success = !isnan(1 / (sum + 1));
    }
  } // end else

  // Return
  return success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool ActiveFaces::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
  double f = 0.0;
  int index = 0;
  double sum = 0.0;
  double temp = 0.0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
    temp = log(fabs(x[i]) + 1.0);
    if (temp > f) {
      f = temp;
      index = i;
    } // end if
    sum = sum + x[i];
  } // end for
  temp = log(fabs(sum) + 1.0);
  if (temp > f) {
    f = temp;
    index = -1;
  } // end if

  // Declare success
  bool success = true;

  // Evaluate gradient
  if (index >= 0) {
    if (x[index] >= 0) {
      g[index] = 1 / (x[index] + 1);
    }
    else {
      g[index] = -1 / (-x[index] + 1);
    }
    success = !isnan(g[index]);
  } // end if
  else {
    if (-sum >= 0) {
      for (int i = 0; i < n; i++) {
        g[i] = -1 / (-sum + 1);
      }
      success = !isnan(-1 / (-sum + 1));
    }
    else {
      for (int i = 0; i < n; i++) {
        g[i] = 1 / (sum + 1);
      }
      success = !isnan(1 / (sum + 1));
    }
  } // end else

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool ActiveFaces::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
