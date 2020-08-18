// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptBLASLAPACK.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptVector.hpp"

namespace NonOpt
{

// Constructor with given length; values initialized to zero
Vector::Vector(int length)
  : length_(length),
    max_computed_(true),
    min_computed_(true),
    norm1_computed_(true),
    norm2_computed_(true),
    normInf_computed_(true),
    max_value_(0.0),
    min_value_(0.0),
    norm1_value_(0.0),
    norm2_value_(0.0),
    normInf_value_(0.0)
{

  // Allocate array
  values_ = new double[length];

  // Set inputs for BLASLAPACK
  double value = 0.0;
  int increment1 = 0;
  int increment2 = 1;

  // Initialize values
  dcopy_(&length, &value, &increment1, values_, &increment2);

} // end constructor

// Constructor with given length; values initialized to given value
Vector::Vector(int length,
               double value)
  : length_(length),
    max_computed_(true),
    min_computed_(true),
    norm1_computed_(true),
    norm2_computed_(true),
    normInf_computed_(true),
    max_value_(value),
    min_value_(value)
{

  // Allocate array
  values_ = new double[length];

  // Set inputs for BLASLAPACK
  int increment1 = 0;
  int increment2 = 1;

  // Initialize values
  dcopy_(&length, &value, &increment1, values_, &increment2);

  // Compute norms
  norm1_value_ = (double)length * fabs(value);
  norm2_value_ = sqrt((double)length * pow(value,2.0));
  normInf_value_ = fabs(value);

} // end constructor

// Destructor; values array deleted
Vector::~Vector()
{

  // Delete array
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if

} // end destructor

// Print array with given name
void Vector::print(const Reporter* reporter,
                   std::string name) const
{

  // Print elements
  for (int i = 0; i < length_; i++) {
    reporter->printf(R_NL, R_BASIC, "%s[%6d]=%+23.16e\n", name.c_str(), i, values_[i]);
  }

} // end print

// Make new Vector as a copy
std::shared_ptr<Vector> Vector::makeNewCopy() const
{

  // Create new vector
  std::shared_ptr<Vector> vector(new Vector(length_));

  // Copy elements
  vector->copy(*this);

  // Return
  return vector;

} // end makeNewCopy

// Make new Vector by adding "scalar1" times this Vector to "scalar2" times other_vector
std::shared_ptr<Vector> Vector::makeNewLinearCombination(double scalar1,
                                                         double scalar2,
                                                         const Vector& other_vector) const
{

  // Create new vector
  std::shared_ptr<Vector> vector(new Vector(length_));

  // Copy + add elements
  vector->linearCombination(scalar1, *this, scalar2, other_vector);

  // Return
  return vector;

} // end makeNewLinearCombination

// Set length and initialize values to zero
void Vector::setLength(int length)
{

  // Store length
  length_ = length;

  // Delete previous array, if exists
  if (values_ != nullptr) {
    delete[] values_;
    values_ = nullptr;
  } // end if

  // Allocate array
  values_ = new double[length];

  // Set inputs for BLASLAPACK
  double value = 0.0;
  int increment1 = 0;
  int increment2 = 1;

  // Initialize values
  dcopy_(&length, &value, &increment1, values_, &increment2);

  // Compute scalar values
  max_computed_ = true;
  min_computed_ = true;
  norm1_computed_ = true;
  norm2_computed_ = true;
  normInf_computed_ = true;
  max_value_ = 0.0;
  min_value_ = 0.0;
  norm1_value_ = 0.0;
  norm2_value_ = 0.0;
  normInf_value_ = 0.0;

} // end setLength

// Set element with given index to given value
void Vector::set(int index,
                 double value)
{

  // Asserts
  ASSERT_EXCEPTION(index >= 0, NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Index is negative.");
  ASSERT_EXCEPTION(index < length_, NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Index is too large.");

  // Set value
  values_[index] = value;

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end set

// Copy elements of other_vector
void Vector::copy(const Vector& other_vector)
{

  // Assert
  ASSERT_EXCEPTION(length_ == other_vector.length(), NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Vector length is incorrect.");

  // Set inputs for BLASLAPACK
  int length = length_;
  int increment = 1;

  // Copy elements
  dcopy_(&length, other_vector.values(), &increment, values_, &increment);

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end copy

// Copy elements of double array
void Vector::copyArray(double* array)
{

  // Set inputs for BLASLAPACK
  int length = length_;
  int increment = 1;

  // Copy elements
  dcopy_(&length, array, &increment, values_, &increment);

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end copyArray

// Scale elements by given scalar
void Vector::scale(double scalar)
{

  // Set inputs for BLASLAPACK
  int length = length_;
  int increment = 1;

  // Scale elements
  if (scalar != 0.0) {
    dscal_(&length, &scalar, values_, &increment);
  }
  else {
    int incrementzero = 0;
    dcopy_(&length, &scalar, &incrementzero, values_, &increment);
  } // end else

  // Compute scalar values
  if (max_computed_) {
    if (scalar >= 0.0) {
      max_value_ = scalar*max_value_;
    }
    else {
      if (min_computed_) {
        max_value_ = scalar*min_value_;
      }
      else {
        max_computed_ = false;
      }
    }
  } // end if
  if (min_computed_) {
    if (scalar >= 0.0) {
      min_value_ = scalar*min_value_;
    }
    else {
      if (max_computed_) {
        min_value_ = scalar*max_value_;
      }
      else {
        min_computed_ = false;
      }
    }
  } // end if
  if (norm1_computed_) {
    norm1_value_ = fabs(scalar)*norm1_value_;
  }
  if (norm2_computed_) {
    norm2_value_ = fabs(scalar)*norm2_value_;
  }
  if (normInf_computed_) {
    normInf_value_ = fabs(scalar)*normInf_value_;
  }

} // end scale

// Add to this Vector "scalar" times other_vector
void Vector::addScaledVector(double scalar,
                             const Vector& other_vector)
{

  // Assert
  ASSERT_EXCEPTION(length_ == other_vector.length(), NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Vector length is incorrect.");

  // Set inputs for BLASLAPACK
  int length = length_;
  int increment = 1;

  // Add a*vector
  daxpy_(&length, &scalar, other_vector.values(), &increment, values_, &increment);

  // Reset scalar value bools
  max_computed_ = false;
  min_computed_ = false;
  norm1_computed_ = false;
  norm2_computed_ = false;
  normInf_computed_ = false;

} // end addScaledVector

// Set values as linear combination (scalar1*vector1 + scalar2*vector2)
void Vector::linearCombination(double scalar1,
                               const Vector& vector1,
                               double scalar2,
                               const Vector& vector2)
{

  // Asserts
  ASSERT_EXCEPTION(length_ == vector1.length(), NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Vector length is incorrect.");
  ASSERT_EXCEPTION(length_ == vector2.length(), NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Vector length is incorrect.");

  // Set inputs for BLASLAPACK
  int length = length_;
  int increment = 1;

  // Check scalar1
  if (scalar1 != 0.0) {

    // Copy vector1 elements
    dcopy_(&length, vector1.values(), &increment, values_, &increment);

    // Scale if scalar1 not 1.0
    if (scalar1 != 1.0) {
      dscal_(&length, &scalar1, values_, &increment);
    }

    // Add scalar2*vector2
    if (scalar2 != 0.0) {
      daxpy_(&length, &scalar2, vector2.values(), &increment, values_, &increment);
    }

  } // end if

  // Check scalar2
  else if (scalar2 != 0.0) {

    // Copy vector2 elements
    dcopy_(&length, vector2.values(), &increment, values_, &increment);

    // Scale if scalar2 not 1.0
    if (scalar2 != 1.0) {
      dscal_(&length, &scalar2, values_, &increment);
    }

  } // end else if

  else {

    // Set vector to zeros
    int incrementzero = 0;
    dcopy_(&length, &scalar1, &incrementzero, values_, &increment);

  } // end else

  // Reset scalar value bools
  if (scalar1 == 0.0 && scalar2 == 0.0) {
    max_computed_ = true;
    min_computed_ = true;
    norm1_computed_ = true;
    norm2_computed_ = true;
    normInf_computed_ = true;
    max_value_ = 0.0;
    min_value_ = 0.0;
    norm1_value_ = 0.0;
    norm2_value_ = 0.0;
    normInf_value_ = 0.0;
  }
  else {
    max_computed_ = false;
    min_computed_ = false;
    norm1_computed_ = false;
    norm2_computed_ = false;
    normInf_computed_ = false;
  }

} // end linearCombination

// Inner product with other_vector
double Vector::innerProduct(const Vector& other_vector) const
{

  // Assert
  ASSERT_EXCEPTION(length_ == other_vector.length(), NONOPT_VECTOR_ASSERT_EXCEPTION, "Vector assert failed.  Vector length is incorrect.");

  // Set inputs for BLASLAPACK
  int length = length_;
  int increment = 1;

  // Return
  return ddot_(&length, values_, &increment, other_vector.values(), &increment);

} // end innerProduct

// Maximum element
double Vector::max()
{

  // Check if computed
  if (!max_computed_) {

    // Initialize maximum
    max_value_ = values_[0];

    // Determine maximum
    for (int i = 1; i < length_; i++) {
      max_value_ = fmax(max_value_, values_[i]);
    }

    // Set to computed
    max_computed_ = true;

  } // end if

  // Return maximum
  return max_value_;

} // end max

// Minimum element
double Vector::min()
{

  // Check if computed
  if (!min_computed_) {

    // Initialize minimum
    min_value_ = values_[0];

    // Determine minimum
    for (int i = 1; i < length_; i++) {
      min_value_ = fmin(min_value_, values_[i]);
    }

    // Set to computed
    min_computed_ = true;

  } // end if

  // Return minimum
  return min_value_;

} // end max

// 1-norm
double Vector::norm1()
{

  // Check if computed
  if (!norm1_computed_) {

    // Set inputs for BLASLAPACK
    int length = length_;
    int increment = 1;

    // Evaluate 1-norm
    norm1_value_ = dasum_(&length, values_, &increment);

    // Set to computed
    norm1_computed_ = true;

  } // end if

  // Return 1-norm
  return norm1_value_;

} // end norm1

// 2-norm
double Vector::norm2()
{

  // Check if computed
  if (!norm2_computed_) {

    // Set inputs for BLASLAPACK
    int length = length_;
    int increment = 1;

    // Evaluate 2-norm
    norm2_value_ = dnrm2_(&length, values_, &increment);

    // Set to computed
    norm2_computed_ = true;

  } // end if

  // Return 2-norm
  return norm2_value_;

} // end norm2

// inf-norm
double Vector::normInf()
{

  // Check if computed
  if (!normInf_computed_) {

    // Set inputs for BLASLAPACK
    int length = length_;
    int increment = 1;

    // Find index of element with maximum absolute value
    // (returns index from 1,...,length_)
    int i = idamax_(&length, values_, &increment);

    // Evaluate inf-norm
    normInf_value_ = fabs(values_[i - 1]);

    // Set to computed
    normInf_computed_ = true;

  } // end if

  // Return inf-norm
  return normInf_value_;

} // end normInf

} // namespace NonOpt
