// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptBLAS.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptSymmetricMatrixLimitedMemory.hpp"

namespace NonOpt
{

// Destructor
SymmetricMatrixLimitedMemory::~SymmetricMatrixLimitedMemory() {}

// Add options
void SymmetricMatrixLimitedMemory::addOptions(Options* options,
                                              const Reporter* reporter)
{
  // Add double options

  // Add integer options
  options->addIntegerOption(reporter,
                            "SMLM_history",
                            100,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limited-memory history length.\n"
                            "Default value: 100.");

}  // end addOptions

// Set options
void SymmetricMatrixLimitedMemory::setOptions(const Options* options,
                                              const Reporter* reporter)
{

  // Read double options

  // Read integer options
  options->valueAsInteger(reporter, "SMLM_history", history_);

}  // end setOptions

// Initialize
void SymmetricMatrixLimitedMemory::initialize(const Options* options,
                                              Quantities* quantities,
                                              const Reporter* reporter)
{

  // Set as identity
  setAsDiagonal(quantities->numberOfVariables(), 1.0);

}  // end initialize

// Column
void const SymmetricMatrixLimitedMemory::columnHessian(int column_index,
                                                       Vector& column)
{

  // Asserts
  ASSERT_EXCEPTION(false, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Getting column of limited-memory Hessian not yet implemented.");

}  // end columnHessian

// Column
void const SymmetricMatrixLimitedMemory::columnHessianInverse(int column_index,
                                                              Vector& column)
{

  // Asserts
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");
  ASSERT_EXCEPTION(size_ == column.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column has incorrect length.");

  // Initialize bool indicating column computed
  bool column_computed = false;
  int column_computed_index = -1;

  // Loop over computed columns
  for (int i = 0; i < (int)computed_column_indices_.size(); i++) {
    if (computed_column_indices_[i] == column_index) {
      column_computed = true;
      column_computed_index = i;
      break;
    }  // end if
  }    // end for

  // Commpute column
  if (column_computed) {
    column.copy(*computed_columns_[column_computed_index]);
  }
  else {

    // Create intermediate value vectors
    double* a = new double[(int)s_.size()];
    double* b = new double[(int)s_.size()];

    // Initialize column
    column.scale(0.0);
    column.set(column_index, 1.0);

    // Set parameters for daxpy
    double value = 0.0;
    int increment = 1;

    // Two-loop
    for (int i = (int)s_.size() - 1; i >= 0; i--) {
      a[i] = rho_[i] * s_[i]->innerProduct(column);
      value = -a[i];
      daxpy_(&size_, &value, y_[i]->values(), &increment, column.valuesModifiable(), &increment);
    }  // end for
    column.scale(initial_diagonal_value_);
    for (int i = 0; i < (int)s_.size(); i++) {
      b[i] = rho_[i] * y_[i]->innerProduct(column);
      value = a[i] - b[i];
      daxpy_(&size_, &value, s_[i]->values(), &increment, column.valuesModifiable(), &increment);
    }  // end for

    // Save column
    std::shared_ptr<Vector> computed_column = column.makeNewCopy();

    // Push to set
    computed_columns_.push_back(computed_column);
    computed_column_indices_.push_back(column_index);

    // Delete intermediate value vectors
    delete[] a;
    delete[] b;

  }  // end else

}  // end columnHessianInverse

// Element
double const SymmetricMatrixLimitedMemory::elementHessian(int row_index,
                                                          int column_index)
{

  // Asserts
  ASSERT_EXCEPTION(false, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Getting element of limited-memory Hessian not yet implemented.");

}  // end elementHessian

// Element
double const SymmetricMatrixLimitedMemory::elementHessianInverse(int row_index,
                                                                 int column_index)
{

  // Asserts
  ASSERT_EXCEPTION(row_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is negative.");
  ASSERT_EXCEPTION(row_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Row index is too large.");
  ASSERT_EXCEPTION(column_index >= 0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is negative.");
  ASSERT_EXCEPTION(column_index < size_, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Column index is too large.");

  // Initialize bool indicating element in a computed column
  bool element_computed = false;

  // Initialize element
  double element_value = 0.0;

  // Loop over computed columns
  for (int i = 0; i < (int)computed_column_indices_.size(); i++) {
    if (computed_column_indices_[i] == row_index) {
      element_computed = true;
      element_value = computed_columns_[i]->values()[column_index];
      break;
    }  // end if
    if (computed_column_indices_[i] == column_index) {
      element_computed = true;
      element_value = computed_columns_[i]->values()[row_index];
      break;
    }  // end if
  }    // end for

  // Check of element already computed
  if (!element_computed) {
    Vector column_vector(size_);
    columnHessianInverse(column_index, column_vector);
    element_value = column_vector.values()[row_index];
  }  // end if

  // Return element
  return element_value;

}  // end elementHessianInverse

// Inner product
double SymmetricMatrixLimitedMemory::innerProductHessian(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Compute matrix-vector product
  matrixVectorProductHessian(vector, product);

  // Return inner product
  return vector.innerProduct(product);

}  // end innerProductHessian

// Inner product
double SymmetricMatrixLimitedMemory::innerProductHessianInverse(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Compute matrix-vector product
  matrixVectorProductHessianInverse(vector, product);

  // Return inner product
  return vector.innerProduct(product);

}  // end innerProductHessianInverse

// Matrix-vector product
void SymmetricMatrixLimitedMemory::matrixVectorProductHessian(const Vector& vector,
                                                              Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(false, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Matrix-vector product with limited-memory Hessian not yet implemented.");

}  // end matrixVectorProductHessian

// Matrix-vector product
void SymmetricMatrixLimitedMemory::matrixVectorProductHessianInverse(const Vector& vector,
                                                                     Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Create intermediate value vectors
  double* a = new double[(int)s_.size()];
  double* b = new double[(int)s_.size()];

  // Initialize product
  product.copy(vector);

  // Set parameters for daxpy
  double value = 0.0;
  int increment = 1;

  // Two-loop
  for (int i = (int)s_.size() - 1; i >= 0; i--) {
    a[i] = rho_[i] * s_[i]->innerProduct(product);
    value = -a[i];
    daxpy_(&size_, &value, y_[i]->values(), &increment, product.valuesModifiable(), &increment);
  }  // end for
  product.scale(initial_diagonal_value_);
  for (int i = 0; i < (int)s_.size(); i++) {
    b[i] = rho_[i] * y_[i]->innerProduct(product);
    value = a[i] - b[i];
    daxpy_(&size_, &value, s_[i]->values(), &increment, product.valuesModifiable(), &increment);
  }  // end for

  // Delete intermediate value vectors
  delete[] a;
  delete[] b;

}  // end matrixVectorProductHessianInverse

// Set as diagonal matrix
void SymmetricMatrixLimitedMemory::setAsDiagonal(int size,
                                                 double value)
{

  // Assert
  ASSERT_EXCEPTION(value > 0.0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Value is nonpositive.");

  // Set size
  size_ = size;

  // Set initial diagonal value
  initial_diagonal_value_ = value;

  // Clear s, y, and rho sets
  s_.clear();
  y_.clear();
  rho_.clear();

  // Clear computed columns set
  computed_columns_.clear();
  computed_column_indices_.clear();

}  // end setAsDiagonal

// Symmetric update
void SymmetricMatrixLimitedMemory::updateBFGS(const Vector& s,
                                              const Vector& y)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == s.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector s has incorrect length.");
  ASSERT_EXCEPTION(size_ == y.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector y has incorrect length.");

  // Create pointers to new vectors, set as copy of (s,y) pair
  std::shared_ptr<Vector> s_new = s.makeNewCopy();
  std::shared_ptr<Vector> y_new = y.makeNewCopy();

  // Compute rho value
  double rho_new = 1.0 / (s_new->innerProduct(*y_new));

  // Add pair
  s_.push_back(s_new);
  y_.push_back(y_new);
  rho_.push_back(rho_new);

  // Remove old pair
  while ((int)s_.size() >= history_) {
    s_.erase(s_.begin());
    y_.erase(y_.begin());
    rho_.erase(rho_.begin());
  }  // end while

  // Clear computed columns set
  computed_columns_.clear();
  computed_column_indices_.clear();

}  // end updateBFGS

// Print
void SymmetricMatrixLimitedMemory::print(const Reporter* reporter,
                                         std::string name) const
{

  // Print elements
  reporter->printf(R_NL, R_BASIC, "%s initial_diagonal_value=%+23.16e\n", name.c_str(), initial_diagonal_value_);
  reporter->printf(R_QP, R_BASIC, "%s initial_diagonal_value=%+23.16e\n", name.c_str(), initial_diagonal_value_);
  for (int j = 0; j < (int)s_.size(); j++) {
    for (int i = 0; i < size_; i++) {
      reporter->printf(R_NL, R_BASIC, "%s %6d-th pair: s[%6d]=%+23.16e, y[%6d]=%+23.16e\n", name.c_str(), j, i, s_[j]->values()[i], i, y_[j]->values()[i]);
      reporter->printf(R_QP, R_BASIC, "%s %6d-th pair: s[%6d]=%+23.16e, y[%6d]=%+23.16e\n", name.c_str(), j, i, s_[j]->values()[i], i, y_[j]->values()[i]);
    }  // end for
  }    // end for
  for (int j = 0; j < (int)s_.size(); j++) {
    reporter->printf(R_NL, R_BASIC, "%s rho[%6d]=%+23.16e\n", name.c_str(), j, rho_[j]);
    reporter->printf(R_QP, R_BASIC, "%s rho[%6d]=%+23.16e\n", name.c_str(), j, rho_[j]);
  }  // end for

}  // end print

}  // namespace NonOpt
