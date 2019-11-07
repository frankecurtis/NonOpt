// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include <math.h>
#include "setDim.hpp"
#include "Test29_24.hpp"
#include <vector>

// Constructor
Test29_24::Test29_24() {}

// Destructor
Test29_24::~Test29_24() {}

// Number of variables
bool Test29_24::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_24::initialPoint(int n,
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
bool Test29_24::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

	  // Evaluate maximum value

	std::vector<double> term(n,0.0);
	  for (int i = 2; i <= n-1; i++) {
		  term[i-1]=2*x[i-1]+10.0/pow(n+1,2)*sinh(10.0*x[i-1])-x[i-2]-x[i];
	  }
	  term[0]=2*x[0]+10.0/pow(n+1,2)*sinh(10.0*x[0])-x[1];
	  term[n-1]=2*x[n-1]+10.0/pow(n+1,2)*sinh(10.0*x[n-1])-x[n-2];

	  f = -1.0;
	  for (int i = 0; i < n; i++) {
		  f=fmax(f,fabs(term[i]));
	  }

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_24::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
	std::vector<double> term(n,0.0);
	  for (int i = 2; i <= n-1; i++) {
		  term[i-1]=2*x[i-1]+10.0/pow(n+1,2)*sinh(10.0*x[i-1])-x[i-2]-x[i];
	  }
	  term[0]=2*x[0]+10.0/pow(n+1,2)*sinh(10.0*x[0])-x[1];
	  term[n-1]=2*x[n-1]+10.0/pow(n+1,2)*sinh(10.0*x[n-1])-x[n-2];

  int max_ind = 0;
  double max_val=-1.0;

  for (int i = 0; i < n; i++) {
	  g[i] = 0.0;
	  if(fabs(term[i])>max_val){
		  max_val=fabs(term[i]);
		  max_ind=i;
	  }
  }

  if(max_ind==0){
	  if(term[max_ind]>=0){
		  g[max_ind]=2+100.0/pow(n+1,2)*cosh(10.0*x[max_ind]);
		  g[max_ind+1]=-1.0;
	  }
	  else{
		  g[max_ind]=-2-100.0/pow(n+1,2)*cosh(10.0*x[max_ind]);
		  g[max_ind+1]=1.0;
	  }
  }
  else if(max_ind==n-1){
	  if(term[max_ind]>=0){
		  g[max_ind]=2+100.0/pow(n+1,2)*cosh(10.0*x[max_ind]);
		  g[max_ind-1]=-1.0;
	  }
	  else{
		  g[max_ind]=-2-100.0/pow(n+1,2)*cosh(10.0*x[max_ind]);
		  g[max_ind-1]=1.0;
	  }
  }
  else {
	  if(term[max_ind]>=0){
		  g[max_ind]=2+100.0/pow(n+1,2)*cosh(10.0*x[max_ind]);
		  g[max_ind+1]=-1.0;
		  g[max_ind-1]=-1.0;
	  }
	  else{
		  g[max_ind]=-2-100.0/pow(n+1,2)*cosh(10*x[max_ind]);
		  g[max_ind+1]=1.0;
		  g[max_ind-1]=1.0;
	  }
  }







  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_24::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
