// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "MxHilb.hpp"

// Constructor
MxHilb::MxHilb(int n)
  : number_of_variables_(n)
{

  // Declare array for sums
  sum_ = new double[n];

} // end constructor

// Destructor
MxHilb::~MxHilb()
{

  // Delete array
  if (sum_ != nullptr) {
    delete[] sum_;
    sum_ = nullptr;
  } // end if

} // end destructor

// Number of variables
bool MxHilb::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

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

} // end initialPoint

// Objective value
bool MxHilb::evaluateObjective(int n,
                               const double* x,
                               double& f)
{

  // Evaluate objective
  f = 0.0;
  double sum;
  for (int i = 0; i < n; i++) {
    sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum = sum + x[j] / ((double)i + (double)j + 1.0);
    }
    f = fmax(f, fabs(sum));
  } // end for

  // Return
  return !std::isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool MxHilb::evaluateObjectiveAndGradient(int n,
                                          const double* x,
                                          double& f,
                                          double* g)
{

  // Evaluate sums
  f = 0.0;
  int index = 0;
  for (int i = 0; i < n; i++) {
    sum_[i] = 0.0;
    for (int j = 0; j < n; j++) {
      sum_[i] = sum_[i] + x[j] / ((double)i + (double)j + 1.0);
    }
    if (fabs(sum_[i]) > f) {
      f = fabs(sum_[i]);
      index = i;
    } // end if
  }   // end for

  // Declare success
  bool success = !std::isnan(f);

  // Evaluate gradient
  if (sum_[index] >= 0.0) {
    for (int j = 0; j < n; j++) {
      g[j] = 1.0 / ((double)index + (double)j + 1.0);
      if (std::isnan(g[j])) {
        success = false;
      }
    } // end for
  }   // end if
  else {
    for (int j = 0; j < n; j++) {
      g[j] = -1.0 / ((double)index + (double)j + 1.0);
      if (std::isnan(g[j])) {
        success = false;
      }
    } // end for
  }   // end else

  // Return
  return success;

} // end evaluateObjectiveAndGradient

// Gradient value
bool MxHilb::evaluateGradient(int n,
                              const double* x,
                              double* g)
{

  // Evaluate sums
  double f = 0.0;
  int index = 0;
  for (int i = 0; i < n; i++) {
    sum_[i] = 0.0;
    for (int j = 0; j < n; j++) {
      sum_[i] = sum_[i] + x[j] / ((double)i + (double)j + 1.0);
    }
    if (fabs(sum_[i]) > f) {
      f = fabs(sum_[i]);
      index = i;
    } // end if
  }   // end for

  // Declare success
  bool success = true;

  // Evaluate gradient
  if (sum_[index] >= 0.0) {
    for (int j = 0; j < n; j++) {
      g[j] = 1.0 / ((double)index + (double)j + 1.0);
      if (std::isnan(g[j])) {
        success = false;
      }
    } // end for
  }   // end if
  else {
    for (int j = 0; j < n; j++) {
      g[j] = -1.0 / ((double)index + (double)j + 1.0);
      if (std::isnan(g[j])) {
        success = false;
      }
    } // end for
  }   // end else

  // Return
  return success;

} // end evaluateGradient

// Finalize solution
bool MxHilb::finalizeSolution(int n,
                              const double* x,
                              double f,
                              const double* g)
{
  return true;
}
