// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "ImageDenoising.hpp"

// Constructor
ImageDenoising::ImageDenoising(int rows,
                               int cols,
                               double* image,
                               double regularization,
                               double smoothing)
  : number_of_variables_(rows*cols),
    rows_(rows),
    cols_(cols),
    regularization_(regularization),
    smoothing_(smoothing)
{

  // Initialize image array
  image_ = new double[rows*cols];

  // Evaluate data-fitting term
  for (int i = 0; i < rows*cols; i++) {
    image_[i] = image[i];
  }

} // end ImageDenoising

// Destructor
ImageDenoising::~ImageDenoising()
{

  // Delete array
  if (image_ != nullptr) {
    delete[] image_;
    image_ = nullptr;
  } // end if

} // end ~ImageDenoising

// Number of variables
bool ImageDenoising::numberOfVariables(int& n)
{

  // Set number of variables
  n = number_of_variables_;

  // Return
  return true;

} // end numberOfVariables

// Initial point
bool ImageDenoising::initialPoint(int n,
                                  double* x)
{

  // Set initial point
  for (int i = 0; i < n; i++) {
    x[i] = image_[i];
  }

  // Return
  return true;

} // end initialPoint

// Objective value
bool ImageDenoising::evaluateObjective(int n,
                                       const double* x,
                                       double& f)
{

  // Initialize objective value
  f = 0.0;

  // Evaluate data-fitting term
  for (int i = 0; i < n; i++) {
    f += pow(x[i] - image_[i], 2.0);
  }

  // Evaluate total variation term
  for (int r = 0; r < rows_ - 1; r++) {
    for (int c = 0; c < cols_ - 1; c++) {
      f += regularization_ * sqrt(smoothing_ + pow(x[r*cols_ + c + 1] - x[r*cols_ + c], 2.0) + pow(x[(r+1)*cols_ + c] - x[r*cols_ + c], 2.0));
    }
  }

  // Return
  return !isnan(f);

} // end evaluateObjective

// Objective and gradient value
bool ImageDenoising::evaluateObjectiveAndGradient(int n,
                                                  const double* x,
                                                  double& f,
                                                  double* g)
{

  // Initialize objective value
  f = 0.0;

  // Evaluate data-fitting term
  for (int i = 0; i < n; i++) {
    f    += pow(x[i] - image_[i], 2.0);
    g[i] = 2 * (x[i] - image_[i]);
  }

  // Evaluate total variation term
  for (int r = 0; r < rows_ - 1; r++) {
    for (int c = 0; c < cols_ - 1; c++) {
      double sum = smoothing_ + pow(x[r*cols_ + c + 1] - x[r*cols_ + c], 2.0) + pow(x[(r+1)*cols_ + c] - x[r*cols_ + c], 2.0);
      f += regularization_ * sqrt(sum);
      g[r*cols_ + c    ] += regularization_ * (2.0 * x[r*cols_ + c] - x[r*cols_ + c + 1] - x[(r+1)*cols_ + c]) / sqrt(sum);
      g[r*cols_ + c + 1] += regularization_ * (x[r*cols_ + c + 1] - x[r*cols_ + c]) / sqrt(sum);
      g[(r+1)*cols_ + c] += regularization_ * (x[(r+1)*cols_ + c] - x[r*cols_ + c]) / sqrt(sum);
    }
  }

  // Return
  return true;

} // end evaluateObjectiveAndGradient

// Gradient value
bool ImageDenoising::evaluateGradient(int n,
                                      const double* x,
                                      double* g)
{

  // Initialize gradient
  for (int i = 0; i < n; i++) {
    g[i] = 2 * (x[i] - image_[i]);
  }

  // Update based on total variation
  for (int r = 0; r < rows_ - 1; r++) {
    for (int c = 0; c < cols_ - 1; c++) {
      double sum = smoothing_ + pow(x[r*cols_ + c + 1] - x[r*cols_ + c], 2.0) + pow(x[(r+1)*cols_ + c] - x[r*cols_ + c], 2.0);
      g[r*cols_ + c    ] += regularization_ * (2.0 * x[r*cols_ + c] - x[r*cols_ + c + 1] - x[(r+1)*cols_ + c]) / sqrt(sum);
      g[r*cols_ + c + 1] += regularization_ * (x[r*cols_ + c + 1] - x[r*cols_ + c]) / sqrt(sum);
      g[(r+1)*cols_ + c] += regularization_ * (x[(r+1)*cols_ + c] - x[r*cols_ + c]) / sqrt(sum);
    }
  }

  // Return
  return true;

} // end evaluateGradient

// Finalize solution
bool ImageDenoising::finalizeSolution(int n,
                                      const double* x,
                                      double f,
                                      const double* g)
{
  return true;
}
