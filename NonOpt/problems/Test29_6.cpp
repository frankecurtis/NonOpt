// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include "setDim.hpp"
#include "Test29_6.hpp"
#include <vector>

// Constructor
Test29_6::Test29_6() {}

// Destructor
Test29_6::~Test29_6() {}

// Number of variables
bool Test29_6::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_6::initialPoint(int n,
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
bool Test29_6::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

  // Evaluate maximum value
  f = 0.0;
  std::vector<double> term(n,0.0);
  for (int i = 1; i < n-1; i++) {
	  term[i]=(3-2*x[i])*x[i]+1-x[i-1]-x[i+1];
  }  // end for
  term[0]=(3-2*x[0])*x[0]+1-x[1];
  term[n-1]=(3-2*x[n-1])*x[n-1]+1-x[n-2];

  for (int i = 0; i < n; i++) {
	  f=fmax(f,fabs(term[i]));
  }
  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_6::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
  int max_ind = 0;
  double max_val=0.0;
  std::vector<double> term(n,0.0);

  for (int i = 1; i < n-1; i++) {
	  term[i]=(3-2*x[i])*x[i]+1-x[i-1]-x[i+1];
  }  // end for
  term[0]=(3-2*x[0])*x[0]+1-x[1];
  term[n-1]=(3-2*x[n-1])*x[n-1]+1-x[n-2];
  for(int i=0;i<n;i++){
	  g[i]=0.0;
	  if(fabs(term[i])>max_val){
		  max_ind=i;
	  }
  }

  if(max_ind==0){
	  if(term[max_ind]>=0){
		  g[max_ind]=3.0-4*x[max_ind];
		  g[max_ind+1]=-1.0;
	  }
	  else{
		  g[max_ind]=-3.0+4*x[max_ind];
		  g[max_ind]=1.0;
	  }
  }
  else if(max_ind==n-1){
	  if(term[max_ind]>=0){
		  g[max_ind]=3.0-4*x[max_ind];
		  g[max_ind-1]=-1.0;
	  }
	  else{
		  g[max_ind]=-3.0+4*x[max_ind];
		  g[max_ind-1]=1.0;
	  }
  }
  else {
	  if(term[max_ind]>=0){
		  g[max_ind]=3.0-4*x[max_ind];
		  g[max_ind+1]=-1.0;
		  g[max_ind-1]=-1.0;
	  }
	  else{
		  g[max_ind]=-3.0+4*x[max_ind];
		  g[max_ind-1]=1.0;
		  g[max_ind+1]=1.0;
	  }
  }









  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_6::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
