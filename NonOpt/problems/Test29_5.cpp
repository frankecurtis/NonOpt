// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include "setDim.hpp"
#include "Test29_5.hpp"
#include <vector>

// Constructor
Test29_5::Test29_5() {}

// Destructor
Test29_5::~Test29_5() {}

// Number of variables
bool Test29_5::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

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

}  // end initialPoint

// Objective value
bool Test29_5::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

  // Evaluate maximum value
  f = 0.0;
  std::vector<double> term(n,0.0);
  for (int i = 0; i < n; i++) {
	  for(int j=0;j<n;j++){
		  term[i]+=x[j]/(i+j+1);
	  }
	  f+=fabs(term[i]);
  }


  // Return
  return true;

}  // end evaluateObjective

// Gradient value
bool Test29_5::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

  // Initialize gradient and evaluate maximum value
	std::vector<double> term(n,0.0);
  for (int i = 0; i < n; i++) {
	  g[i]=0.0;
  }
  for (int i = 0; i < n; i++) {
	 for(int j=0;j<n;j++){
		term[i]+=x[j]/(i+j+1);
	 }
	 if(term[i]>=0){
		 for(int j=0;j<n;j++){
			 g[j]+=1.0/(i+j+1);
		 }
	 }
	 else{
		 for(int j=0;j<n;j++){
			 g[j]+=-1.0/(i+j+1);
		 }
	 }
   }






  // Return
  return true;

}  // end evaluateGradient

// Finalize solution
bool Test29_5::finalizeSolution(int n,
                                   const double* x,
                                   double f,
                                   const double* g)
{
  return true;
}
