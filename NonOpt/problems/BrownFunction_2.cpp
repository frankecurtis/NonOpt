// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>
#include "setDim.hpp"
#include "BrownFunction_2.hpp"

// Constructor
BrownFunction_2::BrownFunction_2() {}

// Destructor
BrownFunction_2::~BrownFunction_2() {}

// Number of variables
bool BrownFunction_2::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

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

}  // end initialPoint

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
  return true;

}  // end evaluateObjective

// Gradient value
bool BrownFunction_2::evaluateGradient(int n,
                                       const double* x,
                                       double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0;
  }

  // Evaluate gradient
  for (int i = 0; i < n - 1; i++) {
    if (x[i] >= 0.0 && x[i + 1] >= 0.0) {
      g[i] = g[i] + (pow(x[i + 1], 2) + 1.0) * pow(x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(x[i + 1]) * pow(x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] = g[i + 1] + 2.0 * x[i + 1] * log(x[i]) * pow(x[i], pow(x[i + 1], 2) + 1.0) + (pow(x[i], 2) + 1.0) * pow(x[i + 1], pow(x[i], 2));
    }  // end if
    else if (x[i] >= 0.0 && x[i + 1] < 0.0) {
      g[i] = g[i] + (pow(x[i + 1], 2) + 1.0) * pow(x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(-x[i + 1]) * pow(-x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] = g[i + 1] + 2.0 * x[i + 1] * log(x[i]) * pow(x[i], pow(x[i + 1], 2) + 1.0) - (pow(x[i], 2) + 1.0) * pow(-x[i + 1], pow(x[i], 2));
    }  // end else if
    else if (x[i] < 0.0 && x[i + 1] >= 0.0) {
      g[i] = g[i] - (pow(x[i + 1], 2) + 1.0) * pow(-x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(x[i + 1]) * pow(x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] = g[i + 1] + 2.0 * x[i + 1] * log(-x[i]) * pow(-x[i], pow(x[i + 1], 2) + 1.0) + (pow(x[i], 2) + 1.0) * pow(x[i + 1], pow(x[i], 2));
    }  // end else if
    else {
      g[i] = g[i] - (pow(x[i + 1], 2) + 1.0) * pow(-x[i], pow(x[i + 1], 2)) + 2.0 * x[i] * log(-x[i + 1]) * pow(-x[i + 1], pow(x[i], 2) + 1.0);
      g[i + 1] = g[i + 1] + 2.0 * x[i + 1] * log(-x[i]) * pow(-x[i], pow(x[i + 1], 2) + 1.0) - (pow(x[i], 2) + 1.0) * pow(-x[i + 1], pow(x[i], 2));
    }  // end else
  }    // end for

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool BrownFunction_2::finalizeSolution(int n,
                                       const double* x,
                                       double f,
                                       const double* g)
{
  return true;
}
