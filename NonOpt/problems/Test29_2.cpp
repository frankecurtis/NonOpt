// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include "setDim.hpp"
#include "Test29_2.hpp"

// Constructor
Test29_2::Test29_2() {}

// Destructor
Test29_2::~Test29_2() {}

// Number of variables
bool Test29_2::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_2::initialPoint(int n,
                               double* x)
{

  // Set initial point
  for (int i = 1; i <=n/2; i++) {
    x[i-1] = i/(double)n;
  }
  for (int i = n/2+1; i <=n; i++) {
    x[i-1] = -i/(double)n;
  }

  // Return
  return true;

}  // end initialPoint

// Objective value
bool Test29_2::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

  // Evaluate maximum value
  f = -1.0;
  for (int i = 0; i < n; i++) {
	  if(fabs(x[i])>f){
		  f=fabs(x[i]);
	  }
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_2::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
  int max_ind = 0;
  double max_val=-1.0;
  for (int i = 0; i < n; i++) {
	  if(fabs(x[i])>max_val){
		  max_val=fabs(x[i]);
		  max_ind=i;
	  }
	  g[i] = 0.0;
  }
  if(x[max_ind]>=0){
	  g[max_ind]=1.0;
  }
  else{
	  g[max_ind]=-1.0;
  }






  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_2::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
