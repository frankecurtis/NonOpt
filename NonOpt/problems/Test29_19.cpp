// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include "setDim.hpp"
#include "Test29_19.hpp"
#include <vector>

// Constructor
Test29_19::Test29_19() {}

// Destructor
Test29_19::~Test29_19() {}

// Number of variables
bool Test29_19::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_19::initialPoint(int n,
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
bool Test29_19::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

	  // Evaluate maximum value

	  std::vector<double> term(n,0.0);
	  for (int i = 1; i < n-1; i++) {
		  term[i]=pow(((3-2*x[i])*x[i]+1-x[i-1]-2*x[i+1]),2);
	  }
	  term[0]=(3-2*x[0])*x[0]+1-2*pow(x[1],2);
	  term[n-1]=(3-2*x[n-1])*x[n-1]+1-pow(x[n-2],2);

	  f = x[0];
	  for (int i = 1; i < n; i++) {
		  f=fmax(f,term[i]);
	  }

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_19::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
	  std::vector<double> term(n,0.0);
	  for (int i = 1; i < n-1; i++) {
		  term[i]=pow(((3-2*x[i])*x[i]+1-x[i-1]-2*x[i+1]),2);
	  }
	  term[0]=(3-2*x[0])*x[0]+1-2*pow(x[1],2);
	  term[n-1]=(3-2*x[n-1])*x[n-1]+1-pow(x[n-2],2);

  int max_ind = 0;
  int max_val=term[0];
  g[0]=0.0;
  for (int i = 1; i < n; i++) {
	  g[i] = 0.0;
	  if(term[i]>max_val){
		  max_ind=i;
	  }
  }

  if(max_ind==0){
	  g[max_ind]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind+1])*(3-4*x[max_ind]);
	  g[max_ind+1]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind+1])*(-2);
  }
  else if(max_ind==n-1){
	  g[max_ind]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind+1])*(3-4*x[max_ind]);
	  g[max_ind-1]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind-1])*(-1);
  }
  else{
	  g[max_ind]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind+1])*(3-4*x[max_ind]);
	  g[max_ind+1]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind+1])*(-1);
	  g[max_ind-1]=2*((3-2*x[max_ind])*x[max_ind]+1-2*x[max_ind-1])*(-1);
  }







  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_19::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
