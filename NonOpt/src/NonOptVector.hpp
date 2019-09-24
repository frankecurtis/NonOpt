// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTVECTOR_HPP__
#define __NONOPTVECTOR_HPP__

#include <string>

#include "NonOptReporter.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Reporter;

/**
 * Vector class
 */
class Vector
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor; length and values not initialized
   */
  Vector()
      : values_(nullptr),
        length_(-1){};
  /**
   * Constructor with given length; values initialized to zero
   * \param[in] length is length of Vector to construct
   */
  Vector(int length);
  /**
   * Constructor with given length; values initialized to given value
   * \param[in] length is length of Vector to construct
   * \param[in] value is value at which to set all elements
   */
  Vector(int length,
         double value);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor; values array deleted
   */
  ~Vector();
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print array
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of Vector to print
   */
  void print(const Reporter *reporter,
             std::string name) const;
  //@}

  /** @name Make methods */
  //@{
  /**
   * Make new Vector as a copy
   * \return is pointer to new Vector
   */
  std::shared_ptr<Vector> makeNewCopy() const;
  /**
   * Make new Vector by adding "scalar1" times this Vector to "scalar2" times other_vector
   * \param[in] scalar1 is scalar value for linear combination
   * \param[in] scalar2 is scalar value for linear combination
   * \param[in] other_vector is reference to other Vector
   * \return pointer to new Vector
   */
  std::shared_ptr<Vector> makeNewLinearCombination(double scalar1,
                                                   double scalar2,
                                                   const Vector &other_vector) const;
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get length
   * \return is length of Vector
   */
  inline int const length() const { return length_; };
  /**
   * Get values (const)
   * \return is pointer to array of Vector values
   */
  inline double *const values() const { return values_; };
  /**
   * Get values (modifiable)
   * \return is pointer to array of Vector values (to allow modification of array)
   */
  inline double *valuesModifiable() { return values_; };
  //@}

  /** @name Modify methods */
  //@{
  /**
   * Set length and initialize values to zero
   * \param[in] length is length of Vector to set
   */
  void setLength(int length);
  /**
   * Set element with given index to given value
   * \param[in] index is index of value to set
   * \param[in] value is value to set
   */
  void set(int index,
           double value);
  /**
   * Copy elements of other Vector
   * \param[in] other_vector is reference to other Vector to copy
   */
  void copy(const Vector &other_vector);
  /**
   * Copy elements of double array
   * \param[in] array is array of double values to copy
   */
  void copyArray(double *array);
  /**
   * Scale elements by given scalar
   * \param[in] scalar is scalar for scaling
   */
  void scale(double scalar);
  /**
   * Add to this Vector "scalar" times given vector
   * \param[in] scalar is scalar for value for linear combination
   * \param[in] other_vector is reference to other Vector
   */
  void addScaledVector(double scalar,
                       const Vector &other_vector);
  /**
   * Set values as linear combination (scalar1*vector1 + scalar2*vector2)
   * \param[in] scalar1 is scalar value for linear combination
   * \param[in] scalar2 is scalar value for linear combination
   * \param[in] vector1 is reference to other Vector
   * \param[in] vector2 is reference to other Vector
   */
  void linearCombination(double scalar1,
                         const Vector &vector1,
                         double scalar2,
                         const Vector &vector2);
  //@}

  /** @name Scalar functions */
  //@{
  /**
   * Inner product with given vector
   * \param[in] vector is reference to other Vector
   */
  double innerProduct(const Vector &vector) const;
  /**
   * maximum value
   */
  double min() const;
  /**
   * minimum value
   */
  double max() const;
  /**
   * 1-norm
   */
  double norm1() const;
  /**
   * 2-norm
   */
  double norm2() const;
  /**
   * inf-norm
   */
  double normInf() const;
  //@}

 private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Vector(const Vector &);
  /**
   * Overloaded equals operator
   */
  void operator=(const Vector &);
  //@}

  /** @name Private members */
  //@{
  double *values_; /**< Double array */
  int length_;     /**< Length of array */
  //@}

};  // end Vector

}  // namespace NonOpt

#endif /* __VECTOR_HPP__ */
