// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_5.hpp"

// Constructor
Test29_5::Test29_5(int n)
  : number_of_variables_(n) {}

// Destructor
Test29_5::~Test29_5() {}

// Number of variables
bool Test29_5::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool Test29_5::initialPoint(int n,
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
bool Test29_5::evaluateObjective(int n,
                                 const double* x,
                                 double& f)
{

  // Evaluate objective
  f = 0.0;
  for (int i = 0; i < n; i++) {
    double term = 0.0;
    for (int j = 0; j < n; j++) {
      term += x[j] / double(i + j + 1);
    }
    f += fabs(term);
  } // end for

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool Test29_5::evaluateObjectiveAndGradient(int n,
                                            const double* x,
                                            double& f,
                                            double* g)
{

  // Declare success
  bool success = true;

  // Evaluate objective and gradient
  f = 0.0;
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }
  for (int i = 0; i < n; i++) {
    double term = 0.0;
    for (int j = 0; j < n; j++) {
      term += x[j] / double(i + j + 1);
    }
    f += fabs(term);
    double sign = ((term > 0.0) ? 1.0 : ((term < 0.0) ? -1.0 : 0.0));
    for (int j = 0; j < n; j++) {
      g[j] += sign / double(i + j + 1);
      if (isnan(g[j])) {
        success = false;
      }
    } // end for
  } // end for

  // Return
  return !isnan(f) && success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool Test29_5::evaluateGradient(int n,
                                const double* x,
                                double* g)
{

  // Declare success
  bool success = true;

  // Evaluate gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }
  for (int i = 0; i < n; i++) {
    double term = 0.0;
    for (int j = 0; j < n; j++) {
      term += x[j] / double(i + j + 1);
    }
    double sign = ((term > 0.0) ? 1.0 : ((term < 0.0) ? -1.0 : 0.0));
    for (int j = 0; j < n; j++) {
      g[j] += sign / double(i + j + 1);
      if (isnan(g[j])) {
        success = false;
      }
    }
  } // end for

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool Test29_5::finalizeSolution(int n,
                                const double* x,
                                double f,
                                const double* g)
{
  return true;
}
