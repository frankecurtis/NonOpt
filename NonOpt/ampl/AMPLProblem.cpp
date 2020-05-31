// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cassert>
#include <cmath>
#include <cstring>
#include <memory>

#include "AMPLProblem.hpp"
#include "asl.h"

#define asl cur_ASL

// Constructor
AMPLProblem::AMPLProblem(char* stub)
{

  // Assert that stub has been set
  assert(stub != nullptr);

  // Allocate space for stub
  stub_ = (char*)malloc(strlen(stub) + 1);

  // Set stub_
  strcpy(stub_, stub);

  // Declare nl file
  FILE* nl;

  // Allocate AMPL
  ASL_alloc(ASL_read_fg);

  // Set interface file
  nl = jac0dim(stub, (fint)strlen(stub));

  // Allocate initial point
  X0 = (real*)Malloc(n_var * sizeof(real));

  // Read from nl
  fg_read(nl, 0);

} // end AMPLProblem

// Destructor
AMPLProblem::~AMPLProblem() {
  delete [] stub_;
}

// Number of variables
bool AMPLProblem::numberOfVariables(int& n)
{

  // Set number of variables
  n = n_var;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool AMPLProblem::initialPoint(int n,
                               double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = X0[i];
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool AMPLProblem::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

  // Copy x
  double* x_non_const = new double[n];
  for (int i = 0; i < n; i++) {
    x_non_const[i] = x[i];
  }

  // Evaluate objective
  int nerror = 0;
  f = objval(0, x_non_const, &nerror);

  // Delete x
  delete [] x_non_const;

  // Determine evaluation success
  bool evaluation_success = true;
  if (nerror > 0) {
    evaluation_success = false;
  }

  // Return
  return evaluation_success;

} // end evaluateObjective

// Gradient value
bool AMPLProblem::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Copy x
  double* x_non_const = new double[n];
  for (int i = 0; i < n; i++) {
    x_non_const[i] = x[i];
  }

  // Evaluate gradient
  int nerror = 0;
  objgrd(0, x_non_const, g, &nerror);

  // Delete x
  delete [] x_non_const;

  // Determine evaluation success
  bool evaluation_success = true;
  if (nerror > 0) {
    evaluation_success = false;
  }

  // Return
  return evaluation_success;

} // end evaluateGradient

// Finalize solution
bool AMPLProblem::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
