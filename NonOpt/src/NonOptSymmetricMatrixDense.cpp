// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptSymmetricMatrixDense.hpp"
#include "NonOptBLASLAPACK.hpp"
#include "NonOptDeclarations.hpp"

namespace NonOpt
{

// Destructor
SymmetricMatrixDense::~SymmetricMatrixDense()
{

  // Delete arrays
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if
  if (values_of_inverse_ != nullptr) {
    delete[] values_of_inverse_;
    values_of_inverse_ = nullptr;
  } // end if

} // end destructor

// Add options
void SymmetricMatrixDense::addOptions(Options* options)
{

  // Add double options

  // Add integer options

} // end addOptions

// Set options
void SymmetricMatrixDense::setOptions(Options* options)
{

  // Read double options

  // Read integer options

  // Read string options
  options->valueAsString("approximate_hessian_update", type_);

} // end setOptions

// Initialize
void SymmetricMatrixDense::initialize(const Options* options,
                                      Quantities* quantities,
                                      const Reporter* reporter)
{

  // Set as identity
  setAsDiagonal(quantities->numberOfVariables(), 1.0);

} // end initialize

// Column
void const SymmetricMatrixDense::column(int column_index,
                                        Vector& column)
{

  // Asserts
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");
  ASSERT_EXCEPTION(size_ == column.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column has incorrect length.");

  // Set inputs for BLASLAPACK
  int length = column_index + 1;
  int increment1 = size_;
  int increment2 = 1;

  // Copy first column_index+1 elements
  dcopy_(&length, &values_[column_index], &increment1, column.valuesModifiable(), &increment2);

  // Set input for BLASLAPACK
  length = size_ - (column_index + 1);

  // Copy remaining size_-(column_index+1) elements
  if (length > 0) {
    dcopy_(&length, &values_[column_index * (size_ + 1) + 1], &increment2, &column.valuesModifiable()[column_index + 1], &increment2);
  }

} // end column

// Column of inverse
void const SymmetricMatrixDense::columnOfInverse(int column_index,
                                                 Vector& column)
{

  // Asserts
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");
  ASSERT_EXCEPTION(size_ == column.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column has incorrect length.");

  // Set inputs for BLASLAPACK
  int length = column_index + 1;
  int increment1 = size_;
  int increment2 = 1;

  // Copy elements
  dcopy_(&length, &values_of_inverse_[column_index], &increment1, column.valuesModifiable(), &increment2);

  // Set input for BLASLAPACK
  length = size_ - (column_index + 1);

  // Copy remaining size_-(column_index+1) elements
  if (length > 0) {
    dcopy_(&length, &values_of_inverse_[column_index * (size_ + 1) + 1], &increment2, &column.valuesModifiable()[column_index + 1], &increment2);
  }

} // end columnOfInverse

// Element
double const SymmetricMatrixDense::element(int row_index,
                                           int column_index)
{

  // Asserts
  ASSERT_EXCEPTION(row_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is negative.");
  ASSERT_EXCEPTION(row_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is too large.");
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");

  // Return element
  if (row_index > column_index) {
    int temp = row_index;
    row_index = column_index;
    column_index = temp;
  } // end if
  return values_[row_index * size_ + column_index];

} // end element

// Element of inverse
double const SymmetricMatrixDense::elementOfInverse(int row_index,
                                                    int column_index)
{

  // Asserts
  ASSERT_EXCEPTION(row_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is negative.");
  ASSERT_EXCEPTION(row_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is too large.");
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");

  // Return element
  if (row_index > column_index) {
    int temp = row_index;
    row_index = column_index;
    column_index = temp;
  } // end if
  return values_of_inverse_[row_index * size_ + column_index];

} // end elementOfInverse

// Inner product
double SymmetricMatrixDense::innerProduct(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, values_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

  // Return product
  return ddot_(&size_, product.values(), &increment, vector.values(), &increment);

} // end innerProduct

// Inner product
double SymmetricMatrixDense::innerProductOfInverse(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, values_of_inverse_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

  // Return product
  return ddot_(&size_, product.values(), &increment, vector.values(), &increment);

} // end innerProductOfInverse

// Matrix-vector product
void SymmetricMatrixDense::matrixVectorProduct(const Vector& vector,
                                               Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, values_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

} // end matrixVectorProduct

// Matrix-vector product of inverse
void SymmetricMatrixDense::matrixVectorProductOfInverse(const Vector& vector,
                                                        Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, values_of_inverse_, &size_, vector.values(), &increment, &scale2, product.valuesModifiable(), &increment);

} // end matrixVectorProductOfInverse

// Set as diagonal matrix
void SymmetricMatrixDense::setAsDiagonal(int size,
                                         double value)
{

  // Assert
  ASSERT_EXCEPTION(value >= 0.0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Value is nonpositive.");

  // Check current size
  if (size_ != size) {

    // Delete previous array, if exists
    if (values_ != nullptr) {
      delete[] values_;
      values_ = nullptr;
    } // end if
    if (values_of_inverse_ != nullptr) {
      delete[] values_of_inverse_;
      values_of_inverse_ = nullptr;
    } // end if

    // Set size
    size_ = size;

    // Set length
    length_ = size * size;

    // Allocate arrays
    values_ = new double[length_];
    values_of_inverse_ = new double[length_];

  } // end if

  // Set inputs for BLASLAPACK
  double zero_value = 0.0;
  int increment1 = 0;
  int increment2 = 1;

  // Initialize values
  dcopy_(&length_, &zero_value, &increment1, values_, &increment2);
  dcopy_(&length_, &zero_value, &increment1, values_of_inverse_, &increment2);

  // Set diagonal entries
  for (int i = 0; i < length_; i = i + size_ + 1) {
    values_[i] = value;
    values_of_inverse_[i] = 1.0 / value;
  } // end for

} // end setAsDiagonal

// Set element to value
void SymmetricMatrixDense::setElement(int row_index,
                                      int column_index,
                                      double value)
{
  // Asserts
  ASSERT_EXCEPTION(row_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is negative.");
  ASSERT_EXCEPTION(row_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is too large.");
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");
  // Adjust indices, if necessary
  if (row_index > column_index) {
    int temp = row_index;
    row_index = column_index;
    column_index = temp;
  } // end if
  // Set element
  values_[row_index * size_ + column_index] = value;
} // end setElement

// Symmetric update
void SymmetricMatrixDense::update(const Vector& s,
                                  const Vector& y)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == s.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector s has incorrect length.");
  ASSERT_EXCEPTION(size_ == y.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector y has incorrect length.");

  // Call appropriate update method
  if (type_.compare("BFGS") == 0) {
    updateBFGS(s, y);
  }
  else if (type_.compare("DFP") == 0) {
    updateDFP(s, y);
  }

} // end update

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

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, values_, &size_, s.values(), &increment, &scale2, Hs, &increment);
  dsymv_(&upper_lower, &size_, &scale1, values_of_inverse_, &size_, y.values(), &increment, &scale2, Wy, &increment);

  // Declare scalars
  double sHs = ddot_(&size_, s.values(), &increment, Hs, &increment);
  double yWy = ddot_(&size_, y.values(), &increment, Wy, &increment);
  double sy = ddot_(&size_, y.values(), &increment, s.values(), &increment);

  // Set scale
  double scale = -1.0 / sHs;

  // Perform symmetric rank-1 update (to add -H*s*s'*H/(s'*H*s))
  dsyr_(&upper_lower, &size_, &scale, Hs, &increment, values_, &size_);

  // Set scale
  scale = 1.0 / sy;

  // Perform symmetric rank-1 update (to add y*y'/(s'*y))
  dsyr_(&upper_lower, &size_, &scale, y.values(), &increment, values_, &size_);

  // Set scale
  scale = (1.0 + yWy / sy) / sy;

  // Perform symmetric rank-1 update (to add (1+yWy/sy)/sy*s*s')
  dsyr_(&upper_lower, &size_, &scale, s.values(), &increment, values_of_inverse_, &size_);

  // Set input for BLASLAPACK
  scale = -(1.0 / sy);

  // Perform symmetric rank-2 update (to add -(1/sy)*s*Wy'-(1/sy)*Wy*s')
  dsyr2_(&upper_lower, &size_, &scale, s.values(), &increment, Wy, &increment, values_of_inverse_, &size_);

  // Delete intermediate vector
  if (Hs != nullptr) {
    delete[] Hs;
    Hs = nullptr;
  } // end if
  if (Wy != nullptr) {
    delete[] Wy;
    Wy = nullptr;
  } // end if

} // end updateBFGS

// Symmetric update
void SymmetricMatrixDense::updateDFP(const Vector& s,
                                     const Vector& y)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == s.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector s has incorrect length.");
  ASSERT_EXCEPTION(size_ == y.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector y has incorrect length.");

  // Declare temporary vector
  double* Hs = new double[size_];
  double* Wy = new double[size_];

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  double scale1 = 1.0;
  int increment = 1;
  double scale2 = 0.0;

  // Compute matrix-vector product
  dsymv_(&upper_lower, &size_, &scale1, values_, &size_, s.values(), &increment, &scale2, Hs, &increment);
  dsymv_(&upper_lower, &size_, &scale1, values_of_inverse_, &size_, y.values(), &increment, &scale2, Wy, &increment);

  // Declare scalars
  double sHs = ddot_(&size_, s.values(), &increment, Hs, &increment);
  double yWy = ddot_(&size_, y.values(), &increment, Wy, &increment);
  double sy = ddot_(&size_, y.values(), &increment, s.values(), &increment);

  // Set scale
  double scale = -1.0 / yWy;

  // Perform symmetric rank-1 update (to add -W*y*y'*W/(y'*W*y))
  dsyr_(&upper_lower, &size_, &scale, Wy, &increment, values_of_inverse_, &size_);

  // Set scale
  scale = 1.0 / sy;

  // Perform symmetric rank-1 update (to add s*s'/(s'*y))
  dsyr_(&upper_lower, &size_, &scale, s.values(), &increment, values_of_inverse_, &size_);

  // Set scale
  scale = (1.0 + sHs / sy) / sy;

  // Perform symmetric rank-1 update (to add (1+sHs/sy)/sy*s*s')
  dsyr_(&upper_lower, &size_, &scale, y.values(), &increment, values_, &size_);

  // Set input for BLASLAPACK
  scale = -(1.0 / sy);

  // Perform symmetric rank-2 update (to add -(1/sy)*y*Hs'-(1/sy)*Hs*y')
  dsyr2_(&upper_lower, &size_, &scale, y.values(), &increment, Hs, &increment, values_, &size_);

  // Delete intermediate vector
  if (Hs != nullptr) {
    delete[] Hs;
    Hs = nullptr;
  } // end if
  if (Wy != nullptr) {
    delete[] Wy;
    Wy = nullptr;
  } // end if

} // end updateDFP

// Print
void SymmetricMatrixDense::print(const Reporter* reporter,
                                 std::string name) const
{

  // Print elements of symmetric matrix
  reporter->printf(R_NL, R_BASIC, "Matrix:\n");
  reporter->printf(R_QP, R_BASIC, "Matrix:\n");
  for (int i = 0; i < length_; i++) {
    if (row_(i) <= col_(i)) {
      reporter->printf(R_NL, R_BASIC, "%s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), values_[i]);
      reporter->printf(R_QP, R_BASIC, "%s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), values_[i]);
    } // end if
  }   // end for

  // Print elements of symmetric matrix inverse
  reporter->printf(R_NL, R_BASIC, "Matrix inverse:\n");
  reporter->printf(R_QP, R_BASIC, "Matrix inverse:\n");
  for (int i = 0; i < length_; i++) {
    if (row_(i) <= col_(i)) {
      reporter->printf(R_NL, R_BASIC, "%s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), values_of_inverse_[i]);
      reporter->printf(R_QP, R_BASIC, "%s[%6d][%6d]=%+23.16e\n", name.c_str(), row_(i), col_(i), values_of_inverse_[i]);
    } // end if
  }   // end for

} // end print

} // namespace NonOpt
