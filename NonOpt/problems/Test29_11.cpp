// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include "setDim.hpp"
#include <vector>
#include "Test29_11.hpp"

// Constructor
Test29_11::Test29_11() {}

// Destructor
Test29_11::~Test29_11() {}

// Number of variables
bool Test29_11::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_11::initialPoint(int n,
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
bool Test29_11::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

  // Evaluate maximum value
  f = 0.0;
  std::vector<double> term(2*n-2,0.0);
  for (int k = 1; k < 2*n-1; k++) {
	  int i=(int) (k+1)/2;
	  if(k%2==1){
		  term[k-1]=x[i-1]+x[i]*((5-x[i])*x[i]-2)-13;
	  }
	  else{
		  term[k-1]=x[i-1]+x[i]*((1+x[i])*x[i]-14)-29;
	  }
	  f+=fabs(term[k-1]);
  }  // end for

  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_11::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
  std::vector<double> term(2*n-2,0.0);

  for (int k = 1; k < 2*n-1; k++) {
	  int i=(int) (k+1)/2;
	  if(k%2==1){
		  term[k-1]=x[i-1]+x[i]*((5-x[i])*x[i]-2)-13;
	  }
	  else{
		  term[k-1]=x[i-1]+x[i]*((1+x[i])*x[i]-14)-29;
	  }

	  if(term[k-1]>=0){
		  if(k%2==0){
			  g[i-1]+=1.0;
			  g[i]+=(1+x[i])*x[i]-14+x[i]*(2*x[i]+1);
		  }
		  else{
			  g[i-1]+=1.0;
			  g[i]+=(5-x[i])*x[i]-2+x[i]*(-2*x[i]+5);
		  }
	  }
	  else{
		  if(k%2==0){
			  g[i-1]-=1.0;
			  g[i]+=-(1+x[i])*x[i]-14+x[i]*(2*x[i]+1);
		  }
		  else{
			  g[i-1]-=1.0;
			  g[i]+=-(5-x[i])*x[i]-2+x[i]*(-2*x[i]+5);
		  }
	  }
  }  // end for









  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_11::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
