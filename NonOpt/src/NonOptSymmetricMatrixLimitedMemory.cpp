// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptSymmetricMatrixLimitedMemory.hpp"
#include "NonOptBLASLAPACK.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"

namespace NonOpt
{

// Destructor
SymmetricMatrixLimitedMemory::~SymmetricMatrixLimitedMemory()
{

  // Delete arrays
  if (compact_form_diagonal_ != nullptr) {
    delete[] compact_form_diagonal_;
    compact_form_diagonal_ = nullptr;
  } // end if
  if (compact_form_factorization_ != nullptr) {
    delete[] compact_form_factorization_;
    compact_form_factorization_ = nullptr;
  } // end if
  if (compact_form_inner_product_ != nullptr) {
    delete[] compact_form_inner_product_;
    compact_form_inner_product_ = nullptr;
  } // end if
  if (compact_form_lower_triangular_ != nullptr) {
    delete[] compact_form_lower_triangular_;
    compact_form_lower_triangular_ = nullptr;
  } // end if

} // end destructor

// Add options
void SymmetricMatrixLimitedMemory::addOptions(Options* options)
{
  // Add double options

  // Add integer options
  options->addIntegerOption("SMLM_history",
                            20,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limited-memory history length.\n"
                            "Default     : 20.");

} // end addOptions

// Set options
void SymmetricMatrixLimitedMemory::setOptions(Options* options)
{

  // Read double options

  // Read integer options
  options->valueAsInteger("SMLM_history", history_);

  // Read string options
  options->valueAsString("approximate_hessian_update", type_);

} // end setOptions

// Initialize
void SymmetricMatrixLimitedMemory::initialize(const Options* options,
                                              Quantities* quantities,
                                              const Reporter* reporter)
{

  // Reduce history to at most number of variables
  history_ = fmin(history_, quantities->numberOfVariables());

  // Set as identity
  setAsDiagonal(quantities->numberOfVariables(), 1.0);

  // Delete array, if it exists
  if (compact_form_diagonal_ != nullptr) {
    delete[] compact_form_diagonal_;
    compact_form_diagonal_ = nullptr;
  } // end if
  if (compact_form_inner_product_ != nullptr) {
    delete[] compact_form_inner_product_;
    compact_form_inner_product_ = nullptr;
  } // end if
  if (compact_form_lower_triangular_ != nullptr) {
    delete[] compact_form_lower_triangular_;
    compact_form_lower_triangular_ = nullptr;
  } // end if

  // Allocate memory for compact form matrices
  compact_form_diagonal_ = new double[history_];
  compact_form_inner_product_ = new double[history_ * history_];
  compact_form_lower_triangular_ = new double[history_ * history_];

} // end initialize

// Column
void const SymmetricMatrixLimitedMemory::column(int column_index,
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
    } // end if
  }   // end for

  // Compute column
  if (column_computed) {
    column.copy(*computed_columns_[column_computed_index]);
  }
  else {

    // Create unit vector
    Vector unit_vector(size_, 0.0);
    unit_vector.set(column_index, 1.0);

    // Perform matrix-vector product
    matrixVectorProduct(unit_vector, column);

    // Save column
    std::shared_ptr<Vector> computed_column = column.makeNewCopy();

    // Push to set
    computed_columns_.push_back(computed_column);
    computed_column_indices_.push_back(column_index);

  } // end else

} // end column

// Column of inverse
void const SymmetricMatrixLimitedMemory::columnOfInverse(int column_index,
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
  for (int i = 0; i < (int)computed_column_indices_of_inverse_.size(); i++) {
    if (computed_column_indices_of_inverse_[i] == column_index) {
      column_computed = true;
      column_computed_index = i;
      break;
    } // end if
  }   // end for

  // Compute column
  if (column_computed) {
    column.copy(*computed_columns_of_inverse_[column_computed_index]);
  }
  else {

    // Create unit vector
    Vector unit_vector(size_, 0.0);
    unit_vector.set(column_index, 1.0);

    // Perform matrix-vector product
    matrixVectorProductOfInverse(unit_vector, column);

    // Save column
    std::shared_ptr<Vector> computed_column = column.makeNewCopy();

    // Push to set
    computed_columns_of_inverse_.push_back(computed_column);
    computed_column_indices_of_inverse_.push_back(column_index);

  } // end else

} // end columnOfInverse

// Element
double const SymmetricMatrixLimitedMemory::element(int row_index,
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
    } // end if
    if (computed_column_indices_[i] == column_index) {
      element_computed = true;
      element_value = computed_columns_[i]->values()[row_index];
      break;
    } // end if
  }   // end for

  // Check if element already computed
  if (!element_computed) {
    Vector column_vector(size_);
    column(column_index, column_vector);
    element_value = column_vector.values()[row_index];
  } // end if

  // Return element
  return element_value;

} // end element

// Element of inverse
double const SymmetricMatrixLimitedMemory::elementOfInverse(int row_index,
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
  for (int i = 0; i < (int)computed_column_indices_of_inverse_.size(); i++) {
    if (computed_column_indices_of_inverse_[i] == row_index) {
      element_computed = true;
      element_value = computed_columns_of_inverse_[i]->values()[column_index];
      break;
    } // end if
    if (computed_column_indices_of_inverse_[i] == column_index) {
      element_computed = true;
      element_value = computed_columns_of_inverse_[i]->values()[row_index];
      break;
    } // end if
  }   // end for

  // Check of element already computed
  if (!element_computed) {
    Vector column_vector(size_);
    columnOfInverse(column_index, column_vector);
    element_value = column_vector.values()[row_index];
  } // end if

  // Return element
  return element_value;

} // end elementOfInverse

// Inner product
double SymmetricMatrixLimitedMemory::innerProduct(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Compute matrix-vector product
  matrixVectorProduct(vector, product);

  // Return inner product
  return vector.innerProduct(product);

} // end innerProduct

// Inner product of inverse
double SymmetricMatrixLimitedMemory::innerProductOfInverse(const Vector& vector)
{

  // Assert
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");

  // Create new vector
  Vector product(size_);

  // Compute matrix-vector product
  matrixVectorProductOfInverse(vector, product);

  // Return inner product
  return vector.innerProduct(product);

} // end innerProductOfInverse

// Matrix-vector product
void SymmetricMatrixLimitedMemory::matrixVectorProduct(const Vector& vector,
                                                       Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Call appropriate update method
  if (type_.compare("BFGS") == 0) {
    matrixVectorProductBFGS(vector, product);
  }
  else if (type_.compare("DFP") == 0) {
    matrixVectorProductDFP(vector, product);
  }

} // end matrixVectorProduct

// Matrix-vector product BFGS
void SymmetricMatrixLimitedMemory::matrixVectorProductBFGS(const Vector& vector,
                                                           Vector& product)
{

  // Check if no pairs in storage
  if ((int)s_.size() == 0) {

    // Compute product as multiple of vector
    product.copy(vector);
    product.scale(initial_diagonal_value_);

  } // end if
  else {

    // Declare factorization size
    int factor_size = 2 * (int)s_.size();
    int factor_length = factor_size * factor_size;

    // Set inputs for BLASLAPACK
    char u = 'U';
    int* v = new int[factor_size];
    int flag = 0;

    // Delete factorization array, if exists
    if (compact_form_factorization_ != nullptr) {
      delete[] compact_form_factorization_;
      compact_form_factorization_ = nullptr;
    } // end if

    // Declare array (matrix size = (2*number of pairs)^2)
    compact_form_factorization_ = new double[factor_length];

    // Fill (1,1) block
    for (int i = 0; i < (int)s_.size(); i++) {
      for (int j = i; j < (int)s_.size(); j++) {
        compact_form_factorization_[i + j * factor_size] = compact_form_inner_product_[i * history_ + j];
      }
    } // end for

    // Fill (1,2) block
    int shift_over = 2 * (int)s_.size() * (int)s_.size();
    for (int i = 0; i < (int)s_.size(); i++) {
      for (int j = 0; j < (int)s_.size(); j++) {
        if (i > j) {
          compact_form_factorization_[shift_over + i + j * factor_size] = compact_form_lower_triangular_[i * history_ + j];
        }
        else {
          compact_form_factorization_[shift_over + i + j * factor_size] = 0.0;
        }
      } // end for
    }   // end for

    // Fill (2,2) block
    int shift_down = (int)s_.size();
    for (int i = 0; i < (int)s_.size(); i++) {
      for (int j = 0; j < (int)s_.size(); j++) {
        if (i == j) {
          compact_form_factorization_[shift_over + shift_down + i + j * factor_size] = -compact_form_diagonal_[i];
        }
        else {
          compact_form_factorization_[shift_over + shift_down + i + j * factor_size] = 0.0;
        }
      } // end for
    }   // end for

    // Set inputs for BLASLAPACK
    double* w = new double[factor_length];
    int l = factor_length;

    // Compute factorization
    dsytrf_(&u, &factor_size, compact_form_factorization_, &factor_size, v, w, &l, &flag);

    // Indicate factorization done
    compact_form_factorized_ = true;

    // Delete temporary array
    if (w != nullptr) {
      delete[] w;
      w = nullptr;
    } // end if

    // Compute product with initial matrix
    product.copy(vector);
    product.scale(initial_diagonal_value_);

    // Compute right product
    Vector right_product(factor_size);
    for (int i = 0; i < factor_size; i++) {
      if (i < factor_size / 2) {
        right_product.set(i, s_[i]->innerProduct(product));
      }
      else {
        right_product.set(i, y_[i - factor_size / 2]->innerProduct(vector));
      }
    } // end for

    // Set inputs for BLASLAPACK
    int rhs = 1;

    // Compute solve with factorization
    dsytrs_(&u, &factor_size, &rhs, compact_form_factorization_, &factor_size, v, right_product.valuesModifiable(), &factor_size, &flag);

    // Complete outer product
    Vector outer_product(size_, 0.0);
    for (int i = 0; i < factor_size; i++) {
      if (i < factor_size / 2) {
        outer_product.addScaledVector(initial_diagonal_value_ * right_product.values()[i], *s_[i]);
      }
      else {
        outer_product.addScaledVector(right_product.values()[i], *y_[i - factor_size / 2]);
      }
    } // end for

    // Compute matrix-vector product
    product.addScaledVector(-1.0, outer_product);

    // Delete temporary array
    if (v != nullptr) {
      delete[] v;
      v = nullptr;
    } // end if

  } // end if

} // end matrixVectorProductBFGS

// Matrix-vector product DFP
void SymmetricMatrixLimitedMemory::matrixVectorProductDFP(const Vector& vector,
                                                          Vector& product)
{

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
    a[i] = rho_[i] * y_[i]->innerProduct(product);
    value = -a[i];
    daxpy_(&size_, &value, s_[i]->values(), &increment, product.valuesModifiable(), &increment);
  } // end for
  product.scale(initial_diagonal_value_);
  for (int i = 0; i < (int)s_.size(); i++) {
    b[i] = rho_[i] * s_[i]->innerProduct(product);
    value = a[i] - b[i];
    daxpy_(&size_, &value, y_[i]->values(), &increment, product.valuesModifiable(), &increment);
  } // end for

  // Delete intermediate value vectors
  if (a != nullptr) {
    delete[] a;
    a = nullptr;
  } // end if
  if (b != nullptr) {
    delete[] b;
    b = nullptr;
  } // end if

} // end matrixVectorProductDFP

// Matrix-vector product
void SymmetricMatrixLimitedMemory::matrixVectorProductOfInverse(const Vector& vector,
                                                                Vector& product)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == vector.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector has incorrect length.");
  ASSERT_EXCEPTION(size_ == product.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Product has incorrect length.");

  // Call appropriate update method
  if (type_.compare("BFGS") == 0) {
    matrixVectorProductOfInverseBFGS(vector, product);
  }
  else if (type_.compare("DFP") == 0) {
    matrixVectorProductOfInverseDFP(vector, product);
  }

} // end matrixVectorProductOfInverse

// Matrix-vector product BFGS
void SymmetricMatrixLimitedMemory::matrixVectorProductOfInverseBFGS(const Vector& vector,
                                                                    Vector& product)
{

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
  } // end for
  product.scale(1.0 / initial_diagonal_value_);
  for (int i = 0; i < (int)s_.size(); i++) {
    b[i] = rho_[i] * y_[i]->innerProduct(product);
    value = a[i] - b[i];
    daxpy_(&size_, &value, s_[i]->values(), &increment, product.valuesModifiable(), &increment);
  } // end for

  // Delete intermediate value vectors
  if (a != nullptr) {
    delete[] a;
    a = nullptr;
  } // end if
  if (b != nullptr) {
    delete[] b;
    b = nullptr;
  } // end if

} // end matrixVectorProductOfInverseBFGS

// Matrix-vector product DFP
void SymmetricMatrixLimitedMemory::matrixVectorProductOfInverseDFP(const Vector& vector,
                                                                   Vector& product)
{

  // Check if no pairs in storage
  if ((int)s_.size() == 0) {

    // Compute product as multiple of vector
    product.copy(vector);
    product.scale(1.0 / initial_diagonal_value_);

  } // end if
  else {

    // Declare factorization size
    int factor_size = 2 * (int)s_.size();
    int factor_length = factor_size * factor_size;

    // Set inputs for BLASLAPACK
    char u = 'U';
    int* v = new int[factor_size];
    int flag = 0;

    // Delete factorization array, if exists
    if (compact_form_factorization_ != nullptr) {
      delete[] compact_form_factorization_;
      compact_form_factorization_ = nullptr;
    } // end if

    // Declare array (matrix size = (2*number of pairs)^2)
    compact_form_factorization_ = new double[factor_length];

    // Fill (1,1) block
    for (int i = 0; i < (int)s_.size(); i++) {
      for (int j = i; j < (int)s_.size(); j++) {
        compact_form_factorization_[i + j * factor_size] = compact_form_inner_product_[i * history_ + j];
      }
    } // end for

    // Fill (1,2) block
    int shift_over = 2 * (int)s_.size() * (int)s_.size();
    for (int i = 0; i < (int)s_.size(); i++) {
      for (int j = 0; j < (int)s_.size(); j++) {
        if (i > j) {
          compact_form_factorization_[shift_over + i + j * factor_size] = compact_form_lower_triangular_[i * history_ + j];
        }
        else {
          compact_form_factorization_[shift_over + i + j * factor_size] = 0.0;
        }
      } // end for
    }   // end for

    // Fill (2,2) block
    int shift_down = (int)s_.size();
    for (int i = 0; i < (int)s_.size(); i++) {
      for (int j = 0; j < (int)s_.size(); j++) {
        if (i == j) {
          compact_form_factorization_[shift_over + shift_down + i + j * factor_size] = -compact_form_diagonal_[i];
        }
        else {
          compact_form_factorization_[shift_over + shift_down + i + j * factor_size] = 0.0;
        }
      } // end for
    }   // end for

    // Set inputs for BLASLAPACK
    double* w = new double[factor_length];
    int l = factor_length;

    // Compute factorization
    dsytrf_(&u, &factor_size, compact_form_factorization_, &factor_size, v, w, &l, &flag);

    // Indicate factorization done
    compact_form_factorized_ = true;

    // Delete temporary array
    if (w != nullptr) {
      delete[] w;
      w = nullptr;
    } // end if

    // Compute product with initial matrix
    product.copy(vector);
    product.scale(1.0 / initial_diagonal_value_);

    // Compute right product
    Vector right_product(factor_size);
    for (int i = 0; i < factor_size; i++) {
      if (i < factor_size / 2) {
        right_product.set(i, y_[i]->innerProduct(product));
      }
      else {
        right_product.set(i, s_[i - factor_size / 2]->innerProduct(vector));
      }
    } // end for

    // Set inputs for BLASLAPACK
    int rhs = 1;

    // Compute solve with factorization
    dsytrs_(&u, &factor_size, &rhs, compact_form_factorization_, &factor_size, v, right_product.valuesModifiable(), &factor_size, &flag);

    // Complete outer product
    Vector outer_product(size_, 0.0);
    for (int i = 0; i < factor_size; i++) {
      if (i < factor_size / 2) {
        outer_product.addScaledVector((1.0 / initial_diagonal_value_) * right_product.values()[i], *y_[i]);
      }
      else {
        outer_product.addScaledVector(right_product.values()[i], *s_[i - factor_size / 2]);
      }
    } // end for

    // Compute matrix-vector product
    product.addScaledVector(-1.0, outer_product);

    // Delete temporary array
    if (v != nullptr) {
      delete[] v;
      v = nullptr;
    } // end if

  } // end if

} // end matrixVectorProductOfInverseDFP

// Set as diagonal matrix
void SymmetricMatrixLimitedMemory::setAsDiagonal(int size,
                                                 double value)
{

  // Assert
  ASSERT_EXCEPTION(value > 0.0, NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Value is nonpositive.");

  // Clear s, y, and rho sets
  s_.clear();
  y_.clear();
  rho_.clear();

  // Clear computed columns sets
  computed_columns_.clear();
  computed_columns_of_inverse_.clear();
  computed_column_indices_.clear();
  computed_column_indices_of_inverse_.clear();

  // Reset factorization indicator
  compact_form_factorized_ = false;

  // Delete array, if it exists
  if (compact_form_diagonal_ != nullptr) {
    delete[] compact_form_diagonal_;
    compact_form_diagonal_ = nullptr;
  } // end if
  if (compact_form_factorization_ != nullptr) {
    delete[] compact_form_factorization_;
    compact_form_factorization_ = nullptr;
  } // end if
  if (compact_form_inner_product_ != nullptr) {
    delete[] compact_form_inner_product_;
    compact_form_inner_product_ = nullptr;
  } // end if
  if (compact_form_lower_triangular_ != nullptr) {
    delete[] compact_form_lower_triangular_;
    compact_form_lower_triangular_ = nullptr;
  } // end if

  // Set size
  size_ = size;

  // Set initial diagonal value
  initial_diagonal_value_ = value;

  // Allocate memory for compact form matrices
  compact_form_diagonal_ = new double[history_];
  compact_form_inner_product_ = new double[history_ * history_];
  compact_form_lower_triangular_ = new double[history_ * history_];

} // end setAsDiagonal

// Symmetric update
void SymmetricMatrixLimitedMemory::update(const Vector& s,
                                          const Vector& y)
{

  // Asserts
  ASSERT_EXCEPTION(size_ == s.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector s has incorrect length.");
  ASSERT_EXCEPTION(size_ == y.length(), NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION, "Symmetric matrix assert failed.  Vector y has incorrect length.");

  // Create pointers to new vectors, set as copy of (s,y) pair
  std::shared_ptr<Vector> s_new = s.makeNewCopy();
  std::shared_ptr<Vector> y_new = y.makeNewCopy();

  // Compute sy product
  double sy_new = s_new->innerProduct(*y_new);

  // Compute rho value
  double rho_new = 1.0 / sy_new;

  // Check whether to remove old information
  bool remove = ((int)s_.size() >= history_);

  // Add pair
  s_.push_back(s_new);
  y_.push_back(y_new);
  rho_.push_back(rho_new);

  // Remove old pair
  if (remove) {
    s_.erase(s_.begin());
    y_.erase(y_.begin());
    rho_.erase(rho_.begin());
  } // end while

  // Update diagonal matrix
  if (remove) {
    for (int i = 0; i < history_ - 1; i++) {
      compact_form_diagonal_[i] = compact_form_diagonal_[i + 1];
    } // end for
  }   // end if
  compact_form_diagonal_[(int)s_.size() - 1] = sy_new;

  // Update inner product matrix
  if (remove) {
    for (int i = 0; i < history_ - 1; i++) {
      for (int j = i; j < history_ - 1; j++) {
        compact_form_inner_product_[i * history_ + j] = compact_form_inner_product_[(i + 1) * history_ + (j + 1)];
      } // end for
    }   // end for
  }     // end if
  if (type_.compare("BFGS") == 0) {
    for (int i = 0; i < (int)s_.size() - 1; i++) {
      compact_form_inner_product_[i * history_ + (int)s_.size() - 1] = initial_diagonal_value_ * s_[i]->innerProduct(*s_new);
    } // end for
    compact_form_inner_product_[((int)s_.size() - 1) * history_ + (int)s_.size() - 1] = initial_diagonal_value_ * pow(s_new->norm2(), 2.0);
  } // end if
  else if (type_.compare("DFP") == 0) {
    for (int i = 0; i < (int)s_.size() - 1; i++) {
      compact_form_inner_product_[i * history_ + (int)s_.size() - 1] = (1.0 / initial_diagonal_value_) * y_[i]->innerProduct(*y_new);
    } // end for
    compact_form_inner_product_[((int)s_.size() - 1) * history_ + (int)s_.size() - 1] = (1.0 / initial_diagonal_value_) * pow(y_new->norm2(), 2.0);
  } // end else if

  // Update lower triangular matrix
  if (remove) {
    for (int i = 1; i < history_ - 1; i++) {
      for (int j = 0; j <= i - 1; j++) {
        compact_form_lower_triangular_[i * history_ + j] = compact_form_lower_triangular_[(i + 1) * history_ + (j + 1)];
      } // end for
    }   // end for
  }     // end if
  if (type_.compare("BFGS") == 0) {
    for (int j = 0; j < (int)s_.size() - 1; j++) {
      compact_form_lower_triangular_[((int)s_.size() - 1) * history_ + j] = s_new->innerProduct(*y_[j]);
    } // end for
  }   // end if
  else if (type_.compare("DFP") == 0) {
    for (int j = 0; j < (int)s_.size() - 1; j++) {
      compact_form_lower_triangular_[((int)s_.size() - 1) * history_ + j] = y_new->innerProduct(*s_[j]);
    } // end for
  }   // end else if

  // Clear computed columns sets
  computed_columns_.clear();
  computed_columns_of_inverse_.clear();
  computed_column_indices_.clear();
  computed_column_indices_of_inverse_.clear();

  // Reset factorization indicator
  compact_form_factorized_ = false;

  // Delete factorization array, if exists
  if (compact_form_factorization_ != nullptr) {
    delete[] compact_form_factorization_;
    compact_form_factorization_ = nullptr;
  } // end if

} // end update

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
    } // end for
  }   // end for
  for (int j = 0; j < (int)s_.size(); j++) {
    reporter->printf(R_NL, R_BASIC, "%s rho[%6d]=%+23.16e\n", name.c_str(), j, rho_[j]);
    reporter->printf(R_QP, R_BASIC, "%s rho[%6d]=%+23.16e\n", name.c_str(), j, rho_[j]);
  } // end for

} // end print

} // namespace NonOpt
