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
                               int regularizer,
                               double regularization,
                               double smoothing)
  : number_of_variables_(rows*cols),
    rows_(rows),
    cols_(cols),
    regularizer_(regularizer),
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
      f += regularization_ * (termFunction(x[r*cols_ + c+1] - x[r*cols_ + c]) + termFunction(x[(r+1)*cols_ + c] - x[r*cols_ + c]));
    }
  }
  for (int r = 0; r < rows_ - 1; r++) {
    f += regularization_ * termFunction(x[r*cols_ + cols_-1] - x[(r+1)*cols_ + cols_-1]);
  }
  for (int c = 0; c < cols_ - 1; c++) {
    f += regularization_ * termFunction(x[(rows_-1)*cols_ + c+1] - x[(rows_-1)*cols_ + c]);
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

  // Update based on total variation
  for (int r = 0; r < rows_ - 1; r++) {
    for (int c = 0; c < cols_ - 1; c++) {
      f += regularization_ * (termFunction(x[r*cols_ + c+1] - x[r*cols_ + c]) + termFunction(x[(r+1)*cols_ + c] - x[r*cols_ + c]));
      if (x[r*cols_ + c+1] > x[r*cols_ + c]) {
        g[r*cols_ + c+1] += regularization_ * termDerivative(x[r*cols_ + c+1] - x[r*cols_ + c]);
        g[r*cols_ + c  ] -= regularization_ * termDerivative(x[r*cols_ + c+1] - x[r*cols_ + c]);
      }
      else if (x[r*cols_ + c+1] < x[r*cols_ + c]) {
        g[r*cols_ + c+1] -= regularization_ * termDerivative(x[r*cols_ + c] - x[r*cols_ + c+1]);
        g[r*cols_ + c  ] += regularization_ * termDerivative(x[r*cols_ + c] - x[r*cols_ + c+1]);
      }
      if (x[(r+1)*cols_ + c] > x[r*cols_ + c]) {
        g[(r+1)*cols_ + c] += regularization_ * termDerivative(x[(r+1)*cols_ + c] - x[r*cols_ + c]);
        g[ r   *cols_ + c] -= regularization_ * termDerivative(x[(r+1)*cols_ + c] - x[r*cols_ + c]);
      }
      else if (x[(r+1)*cols_ + c] < x[r*cols_ + c]) {
        g[(r+1)*cols_ + c] -= regularization_ * termDerivative(x[r*cols_ + c] - x[(r+1)*cols_ + c]);
        g[ r   *cols_ + c] += regularization_ * termDerivative(x[r*cols_ + c] - x[(r+1)*cols_ + c]);
      }
    }
  }
  for (int r = 0; r < rows_ - 1; r++) {
    f += regularization_ * termFunction(x[r*cols_ + cols_ - 1] - x[(r+1)*cols_ + cols_ - 1]);
    if (x[r*cols_ + cols_-1] > x[(r+1)*cols_ + cols_-1]) {
      g[ r   *cols_ + cols_-1] += regularization_ * termDerivative(x[r*cols_ + cols_-1] - x[(r+1)*cols_ + cols_-1]);
      g[(r+1)*cols_ + cols_-1] -= regularization_ * termDerivative(x[r*cols_ + cols_-1] - x[(r+1)*cols_ + cols_-1]);
    }
    else if (x[r*cols_ + cols_-1] < x[(r+1)*cols_ + cols_-1]) {
      g[ r   *cols_ + cols_-1] -= regularization_ * termDerivative(x[(r+1)*cols_ + cols_-1] - x[r*cols_ + cols_-1]);
      g[(r+1)*cols_ + cols_-1] += regularization_ * termDerivative(x[(r+1)*cols_ + cols_-1] - x[r*cols_ + cols_-1]);
    }
  }
  for (int c = 0; c < cols_ - 1; c++) {
    f += regularization_ * termFunction(x[(rows_-1)*cols_ + c + 1] - x[(rows_-1)*cols_ + c]);
    if (x[(rows_-1)*cols_ + c+1] > x[(rows_-1)*cols_ + c]) {
      g[(rows_-1)*cols_ + c+1] += regularization_ * termDerivative(x[(rows_-1)*cols_ + c+1] - x[(rows_-1)*cols_ + c]);
      g[(rows_-1)*cols_ + c  ] -= regularization_ * termDerivative(x[(rows_-1)*cols_ + c+1] - x[(rows_-1)*cols_ + c]);
    }
    else if (x[(rows_-1)*cols_ + c+1] < x[(rows_-1)*cols_ + c]) {
      g[(rows_-1)*cols_ + c+1] -= regularization_ * termDerivative(x[(rows_-1)*cols_ + c] - x[(rows_-1)*cols_ + c+1]);
      g[(rows_-1)*cols_ + c  ] += regularization_ * termDerivative(x[(rows_-1)*cols_ + c] - x[(rows_-1)*cols_ + c+1]);
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
      if (x[r*cols_ + c+1] > x[r*cols_ + c]) {
        g[r*cols_ + c+1] += regularization_ * termDerivative(x[r*cols_ + c+1] - x[r*cols_ + c]);
        g[r*cols_ + c  ] -= regularization_ * termDerivative(x[r*cols_ + c+1] - x[r*cols_ + c]);
      }
      else if (x[r*cols_ + c+1] < x[r*cols_ + c]) {
        g[r*cols_ + c+1] -= regularization_ * termDerivative(x[r*cols_ + c] - x[r*cols_ + c+1]);
        g[r*cols_ + c  ] += regularization_ * termDerivative(x[r*cols_ + c] - x[r*cols_ + c+1]);
      }
      if (x[(r+1)*cols_ + c] > x[r*cols_ + c]) {
        g[(r+1)*cols_ + c] += regularization_ * termDerivative(x[(r+1)*cols_ + c] - x[r*cols_ + c]);
        g[ r   *cols_ + c] -= regularization_ * termDerivative(x[(r+1)*cols_ + c] - x[r*cols_ + c]);
      }
      else if (x[(r+1)*cols_ + c] < x[r*cols_ + c]) {
        g[(r+1)*cols_ + c] -= regularization_ * termDerivative(x[r*cols_ + c] - x[(r+1)*cols_ + c]);
        g[ r   *cols_ + c] += regularization_ * termDerivative(x[r*cols_ + c] - x[(r+1)*cols_ + c]);
      }
    }
  }
  for (int r = 0; r < rows_ - 1; r++) {
    if (x[r*cols_ + cols_-1] > x[(r+1)*cols_ + cols_-1]) {
      g[ r   *cols_ + cols_-1] += regularization_ * termDerivative(x[r*cols_ + cols_-1] - x[(r+1)*cols_ + cols_-1]);
      g[(r+1)*cols_ + cols_-1] -= regularization_ * termDerivative(x[r*cols_ + cols_-1] - x[(r+1)*cols_ + cols_-1]);
    }
    else if (x[r*cols_ + cols_-1] < x[(r+1)*cols_ + cols_-1]) {
      g[ r   *cols_ + cols_-1] -= regularization_ * termDerivative(x[(r+1)*cols_ + cols_-1] - x[r*cols_ + cols_-1]);
      g[(r+1)*cols_ + cols_-1] += regularization_ * termDerivative(x[(r+1)*cols_ + cols_-1] - x[r*cols_ + cols_-1]);
    }
  }
  for (int c = 0; c < cols_ - 1; c++) {
    if (x[(rows_-1)*cols_ + c+1] > x[(rows_-1)*cols_ + c]) {
      g[(rows_-1)*cols_ + c+1] += regularization_ * termDerivative(x[(rows_-1)*cols_ + c+1] - x[(rows_-1)*cols_ + c]);
      g[(rows_-1)*cols_ + c  ] -= regularization_ * termDerivative(x[(rows_-1)*cols_ + c+1] - x[(rows_-1)*cols_ + c]);
    }
    else if (x[(rows_-1)*cols_ + c+1] < x[(rows_-1)*cols_ + c]) {
      g[(rows_-1)*cols_ + c+1] -= regularization_ * termDerivative(x[(rows_-1)*cols_ + c] - x[(rows_-1)*cols_ + c+1]);
      g[(rows_-1)*cols_ + c  ] += regularization_ * termDerivative(x[(rows_-1)*cols_ + c] - x[(rows_-1)*cols_ + c+1]);
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
} // end finalizeSolution

// Term function
double ImageDenoising::termFunction(const double t)
{

  // Switch on regularizer
  switch(regularizer_) {
    case 0:
      return smoothing_ * fabs(t);
      break;
    case 1:
      return smoothing_ * log(1.0 + smoothing_ * fabs(t));
      break;
    case 2:
      return smoothing_ * fabs(t) / (1.0 + smoothing_ * fabs(t));
      break;
    case 3:
      return 0.5 * (smoothing_ - pow(fmax(0.0, smoothing_ - fabs(t)), 2.0) / smoothing_);
      break;
    default:
      return fabs(t);
  } // end switch

} // end termFunction

// Term derivative
double ImageDenoising::termDerivative(const double t)
{

  // Switch on regularizer (t presumed to be nonnegative)
  switch(regularizer_) {
    case 0:
      return smoothing_;
      break;
    case 1:
      return pow(smoothing_, 2.0) / (1.0 + smoothing_ * t);
      break;
    case 2:
      return smoothing_ / pow(1.0 + smoothing_ * t, 2.0);
      break;
    case 3:
      if (t >= smoothing_) {
        return 0.0;
      }
      else {
        return (smoothing_ - t) / smoothing_;
      }
      break;
    default:
      return 1.0;
  } // end switch

} // end termFunction