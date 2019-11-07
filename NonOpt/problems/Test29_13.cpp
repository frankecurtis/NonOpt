// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Minhan Li

#include <cmath>
#include <vector>
#include "NonOptBLAS.hpp"
#include "setDim.hpp"
#include "Test29_13.hpp"

// Constructor
Test29_13::Test29_13() {}

// Destructor
Test29_13::~Test29_13() {}

int sgn(double val) {
   return ((0.0 < val) - (val < 0.0));
}
// Number of variables
bool Test29_13::numberOfVariables(int& n)
{

  // Set number of variables
	setDim di;
  n = di.getDim();


  // Return
  return true;

}  // end numberOfVariables

// Initial point
bool Test29_13::initialPoint(int n,
                               double* x)
{

  // Set initial point
//  for (int i = 0; i < n; i++) {
//    x[i] = 1.0;
//  }
  for(int i=1;i<=n;i++){
	  if(i%4==0){x[i-1]=0.8;}
	  else if(i%4==1){x[i-1]=-0.8;}
	  else if(i%4==2){x[i-1]=1.2;}
	  else if(i%4==3){x[i-1]=-1.2;}
  }

  // Return
  return true;

}  // end initialPoint

// Objective value
bool Test29_13::evaluateObjective(int n,
                                    const double* x,
                                    double& f)
{

 double* xx =new double[n];

 for(int j=0;j<n;j++){
	 xx[j]=x[j];
 }

  // Evaluate maximum value
  f = 0.0;
 double y[4]={-14.4,-6.8,-4.2,-3.2};

  for (int ka = 1; ka <=2*n-4; ka++) {
	  int i=2*(int)(ka+3)/4-2;
	  int l=(ka-1)%4+1;
	  double fa=y[l-1];
	  for(int k=1;k<4;k++){
		  double a=pow(k,2)/(double)l;
		  for(int j=1;j<5;j++){
			  if(xx[i+j-1]==0){
				  //a*=sgn(x[i+j-1])*fabs(pow(1e-16,j/(k*l)));
				  xx[i+j-1]=1e-16;
			  }
			  else{
			  a*=sgn(xx[i+j-1])*fabs(pow(xx[i+j-1],j/(k*l)));
			  }
		  }
		  fa+=a;
	  }
	  f+=fabs(fa);
  }  // end for

  delete []xx;
  // Return
  return true;

}  // end evaluateObjective


// Gradient value
bool Test29_13::evaluateGradient(int n,
                                   const double* x,
                                   double* g)
{

	 double* xx =new double[n];

	 for(int j=0;j<n;j++){
		 xx[j]=x[j];
	 }
  // Initialize gradient and evaluate maximum value

  for (int i = 0; i < n; i++) {
	  g[i] = 0.0;
  }
  double y [4]={-14.4,-6.8,-4.2,-3.2};


  for (int ka = 1; ka <=2*n-4; ka++) {
	  int i=2*((ka+3)/4)-2;
	  int l=(ka-1)%4+1;
	  double fa=y[l-1];
	  std::vector<double> term(3,0.0);
	  for(int k=1;k<4;k++){
		  double a=pow(k,2)/(double)l;
		  for(int j=1;j<5;j++){
			  if(xx[i+j-1]==0){
				  //a*=sgn(xx[i+j-1])*fabs(pow(1e-16,j/(k*l)));
				  xx[i+j-1]=1e-16;
			  }
			  else{
			  a*=sgn(xx[i+j-1])*fabs(pow(xx[i+j-1],j/(k*l)));
			  }
		  }
		  term[k-1]=a;

		  fa+=a;
	  }
	  int ai=sgn(fa);
	  for(int k=1;k<4;k++){
//		  double aa=pow(k,2)/(double)l;
//		  for(int j=1;j<5;j++){
//			  if(xx[i+j-1]==0){
//				  //aa*=sgn(x[i+j-1])*fabs(pow(1e-16,j/(k*l)));
//				  xx[i+j-1]=1e-16;
//			  }
//			  else{
//			  aa*=sgn(xx[i+j-1])*fabs(pow(xx[i+j-1],j/(k*l)));
//			  }
//		  }
		  for(int j=1;j<5;j++){
			  g[i+j-1]+=(j/(k*l))*term[k-1]/(xx[i+j-1])*ai;
		  }

	  }
  }  // end for


  delete []xx;




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
