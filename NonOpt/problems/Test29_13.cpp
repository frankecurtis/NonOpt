// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>

#include "Test29_13.hpp"

// Constructor
Test29_13::Test29_13(int n)
    : number_of_variables_(n) {}

// Destructor
Test29_13::~Test29_13() {}

// Number of variables
bool Test29_13::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_13::initialPoint(int n,
                             double* x)
{

  // Set initial point
  for (int i = 1; i <= n; i++) {
    if (i % 4 == 0) {
      x[i - 1] = 0.8;
    }
    else if (i % 4 == 1) {
      x[i - 1] = -0.8;
    }
    else if (i % 4 == 2) {
      x[i - 1] = 1.2;
    }
    else {
      x[i - 1] = -1.2;
    }
  }  // end for

  // Return
  return true;

}  // end initialPoint

// Objective value
bool Test29_13::evaluateObjective(int n,
                                  const double* x,
                                  double& f)
{

  // Declare parameters
  double y[4] = {-14.4, -6.8, -4.2, -3.2};

  // Evaluate sum of absolute values
  f = 0.0;
  for (int k = 1; k <= 2 * n - 4; k++) {
    int i = 2 * ((k + 3) / 4) - 2;
    int l = (k - 1) % 4 + 1;
    double term1 = y[l - 1];
    for (int h = 1; h <= 3; h++) {
      double term2 = 0.0;
      for (int j = 1; j <= 4; j++) {
        if (x[i + j - 1] == 0.0) {
          term2 *= 0.0;
        }
        else {
          term2 *= (x[i + j - 1] / fabs(x[i + j - 1])) * fabs(pow(x[i + j - 1], (double)j / (double)(h * l)));
        }
      }
      term1 += ((double)(h * h) / (double)l) * term2;
    }
    f += fabs(term1);
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_13::evaluateGradient(int n,
                                 const double* x,
                                 double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 0.0;
  }

  // Declare parameters
  double y[4] = {-14.4, -6.8, -4.2, -3.2};

  // Evaluate gradient
  for (int k = 1; k <= 2 * n - 4; k++) {
    int i = 2 * ((k + 3) / 4) - 2;
    int l = (k - 1) % 4 + 1;
    double term1 = y[l - 1];
    double p[3][4];
    for (int h = 1; h <= 3; h++) {
      double term2 = 0.0;
      for (int j = 1; j <= 4; j++) {
        if (x[i + j - 1] == 0.0) {
          term2 *= 0.0;
        }
        else {
          term2 *= (x[i + j - 1] / fabs(x[i + j - 1])) * fabs(pow(x[i + j - 1], (double)j / (double)(h * l)));
        }
        p[h - 1][j - 1] = term2;
      }  // end for
      term1 += ((double)(h * h) / (double)l) * term2;
    }  // end for
    double sign = ((term1 >= 0.0) ? 1.0 : -1.0);
    for (int h = 1; h <= 3; h++) {
      double product = sign * ((double)(h * h) / (double)l);
      for (int j = 1; j <= 4; j++) {
        for (int m = 1; m <= 4; m++) {
          if (j != m) {
            product *= p[h - 1][m - 1];
          }
        }  // end for
        double sign_x = ((x[i + j - 1] >= 0.0) ? 1.0 : -1.0);
        g[i + j - 1] += product * sign_x * ((double)j / double(h * l)) * pow(fabs(x[i + j - 1]), ((double)j / (double)(h * l) - 1.0));
      }  // end for
    }    // end for
  }      // end for

  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_13::finalizeSolution(int n,
                                 const double* x,
                                 double f,
                                 const double* g)
{
  return true;
}
