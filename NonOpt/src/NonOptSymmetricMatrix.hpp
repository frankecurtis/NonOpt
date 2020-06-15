// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSYMMETRICMATRIX_HPP__
#define __NONOPTSYMMETRICMATRIX_HPP__

#include <memory>
#include <string>
#include <vector>

#include "NonOptEnumerations.hpp"
#include "NonOptOptions.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptStrategy.hpp"
#include "NonOptVector.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;
class Strategy;
class Vector;

/**
 * SymmetricMatrix class
 */
class SymmetricMatrix : public Strategy
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  SymmetricMatrix(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~SymmetricMatrix(){};
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void addOptions(Options* options,
                          const Reporter* reporter) = 0;
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void setOptions(const Options* options,
                          const Reporter* reporter) = 0;
  //@}

  /** @name Initialization method */
  //@{
  /**
   * Initialize strategy
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void initialize(const Options* options,
                          Quantities* quantities,
                          const Reporter* reporter) = 0;
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get column of symmetric matrix
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  virtual void const column(int column_index,
                            Vector& column) = 0;
  /**
   * Get column of symmetric matrix inverse
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  virtual void const columnOfInverse(int column_index,
                                     Vector& column) = 0;
  /**
    * Get element of symmetric matrix
    * \param[in] row_index is row index number
    * \param[in] column_index is column index number
    * \return (row_index,column_index) element of matrix
    */
  virtual double const element(int row_index,
                               int column_index) = 0;
  /**
    * Get element of symmetric matrix inverse
    * \param[in] row_index is row index number
    * \param[in] column_index is column index number
    * \return (row_index,column_index) element of matrix
    */
  virtual double const elementOfInverse(int row_index,
                                        int column_index) = 0;
  /**
   * Get inner product of symmetric matrix with vector
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  virtual double innerProduct(const Vector& vector) = 0;
  /**
   * Get inner product of symmetric matrix inverse with vector
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  virtual double innerProductOfInverse(const Vector& vector) = 0;
  /**
   * Get product of symmetric matrix with vector
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  virtual void matrixVectorProduct(const Vector& vector,
                                   Vector& product) = 0;
  /**
   * Get product of symmetric matrix inverse with vector
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  virtual void matrixVectorProductOfInverse(const Vector& vector,
                                            Vector& product) = 0;
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  /**
    * Get number of rows
    * \return number of rows of the matrix
    */
  virtual int const size() const = 0;
  /**
   * Get status
   * \return current status of symmetric matrix
   */
  inline SM_Status status() { return status_; };
  //@}

  /** @name Modify methods */
  //@{
  /**
    * Set as diagonal matrix
    * \param[in] size is size of matrix to create
    * \param[in] value is value to set in diagonal elements (and set all else zero)
    */
  virtual void setAsDiagonal(int size,
                             double value) = 0;
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(SM_Status status) { status_ = status; };
  /**
    * Update approximation
    * \param[in] s is reference to Vector representing iteration displacement
    * \param[in] y is reference to Vector representing gradient displacement
    */
  virtual void update(const Vector& s,
                      const Vector& y) = 0;
  //@}

  /** @name Print methods */
  //@{
  /**
    * Print array
    * \param[in] reporter is pointer to Reporter object from NonOpt
    * \param[in] name is name of Symmetric Matrix to print
    */
  virtual void print(const Reporter* reporter,
                     std::string name) const = 0;
  //@}

protected:
  /** @name Protected members */
  //@{
  bool initial_scaling_; /**< Indicator of initial scaling */
  std::string type_;     /**< Type of update */
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  SymmetricMatrix(const SymmetricMatrix&);
  /**
   * Overloaded equals operator
   */
  void operator=(const SymmetricMatrix&);
  //@}

  /** @name Private members */
  //@{
  SM_Status status_; /**< Termination status */
  //@}

}; // end SymmetricMatrix

} // namespace NonOpt

#endif /* __NONOPTSYMMETRICMATRIX_HPP__ */
