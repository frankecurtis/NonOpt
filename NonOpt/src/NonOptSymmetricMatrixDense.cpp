// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptBLAS.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptSymmetricMatrixDense.hpp"

namespace NonOpt
{

// Destructor
SymmetricMatrixDense::~SymmetricMatrixDense()
{

  // Delete array
  if (hessian_inverse_values_ != nullptr) {
    delete[] hessian_inverse_values_;
  }

}  // end destructor

// Add options
void SymmetricMatrixDense::addOptions(Options* options,
                                      const Reporter* reporter)
{

  // Add double options

  // Add integer options

}  // end addOptions

// Set options
void SymmetricMatrixDense::setOptions(const Options* options,
                                      const Reporter* reporter)
{

  // Read double options

  // Read integer options

}  // end setOptions

// Initialize
void SymmetricMatrixDense::initialize(const Options* options,
                                      Quantities* quantities,
                                      const Reporter* reporter)
{

  // Set as identity
  setAsDiagonal(quantities->numberOfVariables(), 1.0);

}  // end initialize

// Column
void const SymmetricMatrixDense::columnHessian(int column_index,
                                               Vector& column)
{

  // Asserts
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");
  ASSERT_EXCEPTION(size_ == column.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column has incorrect length.");

  // Set inputs for blas
  int increment1 = size_;
  int increment2 = 1;

  // Copy elements
  dcopy_(&size_, &hessian_values_[column_index], &increment1, column.valuesModifiable(), &increment2);

}  // end columnHessian

// Column
void const SymmetricMatrixDense::columnHessianInverse(int column_index,
                                                      Vector& column)
{

  // Asserts
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");
  ASSERT_EXCEPTION(size_ == column.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column has incorrect length.");

  // Set inputs for blas
  int increment1 = size_;
  int increment2 = 1;

  // Copy elements
  dcopy_(&size_, &hessian_inverse_values_[column_index], &increment1, column.valuesModifiable(), &increment2);

}  // end columnHessianInverse

// Element of Hessian
double const SymmetricMatrixDense::elementHessian(int row_index,
                                                  int column_index)
{

  // Asserts
  ASSERT_EXCEPTION(row_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is negative.");
  ASSERT_EXCEPTION(row_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is too large.");
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");

  // Return element
  return hessian_values_[row_index * size_ + column_index];

}  // end elementHessian

// Element of Inverse Hessian
double const SymmetricMatrixDense::elementHessianInverse(int row_index,
                                                         int column_index)
{

  // Asserts
  ASSERT_EXCEPTION(row_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is negative.");
  ASSERT_EXCEPTION(row_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is too large.");
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");

  // Return element
  return hessian_inverse_values_[row_index * size_ + column_index];

}  // end elementHessianInverse

// Inner product
double SymmetricMatrixDense::innerProductHessian(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Set inputs for blas
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, hessian_values_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

  // Return product
  return ddot_(&size_, product.values(), &increment, vector.values(), &increment);

}  // end innerProductHessian

// Inner product
double SymmetricMatrixDense::innerProductHessianInverse(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Set inputs for blas
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, hessian_inverse_values_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

  // Return product
  return ddot_(&size_, product.values(), &increment, vector.values(), &increment);

}  // end innerProductHessianInverse

void SymmetricMatrixDense::matrixVectorProductHessian(const Vector& vector,
                                                      Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Set inputs for blas
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, hessian_values_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

}  // end matrixVectorProductHessian

// Matrix-vector product
void SymmetricMatrixDense::matrixVectorProductHessianInverse(const Vector& vector,
                                                             Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Set inputs for blas
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, hessian_inverse_values_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

}  // end matrixVectorProductHessianInverse

// Set as diagonal matrix
void SymmetricMatrixDense::setAsDiagonal(int size,
                                         double value)
{

  // Assert
  ASSERT_EXCEPTION(value > 0.0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Value is nonpositive.");

  // Check current size
  if (size_ != size) {

    // Delete previous array, if exists
    if (hessian_values_!= nullptr) {
      delete[] hessian_values_;
    }
    if (hessian_inverse_values_!= nullptr) {
      delete[] hessian_inverse_values_;
    }

    // Store size
    size_ = size;

    // Store length
    length_ = size * size;

    // Allocate array
    hessian_values_ = new double[length_];
    hessian_inverse_values_ = new double[length_];

  }  // end if

  // Set inputs for blas
  double zero_value = 0.0;
  int increment1 = 0;
  int increment2 = 1;

  // Initialize values
  dcopy_(&length_, &zero_value, &increment1, hessian_values_, &increment2);
  dcopy_(&length_, &zero_value, &increment1, hessian_inverse_values_, &increment2);

  // Set diagonal entries
  for (int i = 0; i < length_; i = i + size_ + 1) {
    hessian_values_[i] = value;
   hessian_inverse_values_[i] = 1.0 / value;
  }

}  // end setAsDiagonal

// Symmetric update
void SymmetricMatrixDense::updateBFGS(const Vector& s,
                                      const Vector& y)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == s.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector s has incorrect length.");
  ASSERT_EXCEPTION(size_ == y.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector y has incorrect length.");

  // Declare temporary vector
  double* Hs = new double[size_];
  double* Wy = new double[size_];

  // Set inputs for blas
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, hessian_values_, &size_, s.values(), &increment, &scale2, Hs, &increment);
  dsymv_(&upper_lower, &size_, &scale1, hessian_inverse_values_, &size_, y.values(), &increment, &scale2, Wy, &increment);

  // Declare scalars
  double sHs = ddot_(&size_, s.values(), &increment, Hs, &increment);
  double yWy = ddot_(&size_, y.values(), &increment, Wy, &increment);
  double sy = ddot_(&size_, y.values(), &increment, s.values(), &increment);

  // Set scale
  double scale = -1.0 / sHs;

  // Perform symmetric rank-1 update (to add -H*s*s'*H/(s'*H*s))
  dsyr_(&upper_lower, &size_, &scale, Hs, &increment, hessian_values_, &size_);

  // Set scale
  scale = 1.0 / sy;

  // Perform symmetric rank-1 update (to add y*y'/(s'*y))
  dsyr_(&upper_lower, &size_, &scale, y.values(), &increment, hessian_values_, &size_);

  // Set scale
  scale = (1.0 + yWy / sy) / sy;

  // Perform symmetric rank-1 update (to add (1+yMy/sy)*s*s')
  dsyr_(&upper_lower, &size_, &scale, s.values(), &increment, hessian_inverse_values_, &size_);

  // Set input for lapack
  scale = -(1.0 / sy);

  // Perform symmetric rank-2 update (to add -(1/sy)*s*My'-(1/sy)*My*s')
  dsyr2_(&upper_lower, &size_, &scale, s.values(), &increment, Wy, &increment, hessian_inverse_values_, &size_);

  // Complete matrix
  for (int i = 1; i < size_; i++) {
    for (int j = 0; j < i; j++) {
      hessian_values_[i * size_ + j] = hessian_values_[j * size_ + i];
      hessian_inverse_values_[i * size_ + j] = hessian_inverse_values_[j * size_ + i];
    }
  }  // end for

  // Delete intermediate vector
  if (Hs != nullptr) {
    delete[] Hs;
  }
  if (Wy != nullptr) {
    delete[] Wy;
  }

}  // end updateBFGS

// Print
void SymmetricMatrixDense::print(const Reporter* reporter,
                                 std::string name) const
{

  // Print elements of Hessian
  for (int i = 0; i < length_; i++) {
    reporter->printf(R_NL, R_BASIC, "Hessian %s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), hessian_values_[i]);
    reporter->printf(R_QP, R_BASIC, "Hessian %s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), hessian_values_[i]);
  }  // end for

  // Print elements of inverse Hessian
  for (int i = 0; i < length_; i++) {
    reporter->printf(R_NL, R_BASIC, "Hessian inverse %s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), hessian_inverse_values_[i]);
    reporter->printf(R_QP, R_BASIC, "Hessian inverse %s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), hessian_inverse_values_[i]);
  }  // end for

}  // end print

}  // namespace NonOpt
