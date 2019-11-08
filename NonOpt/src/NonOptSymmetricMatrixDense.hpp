// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSYMMETRICMATRIXDENSE_HPP__
#define __NONOPTSYMMETRICMATRIXDENSE_HPP__

#include "NonOptSymmetricMatrix.hpp"

namespace NonOpt
{

/**
  * SymmetricMatrixDense class
  */
class SymmetricMatrixDense : public SymmetricMatrix
{

 public:
  /** @name Constructors */
  //@{
  /**
    * Construct SymmetricMatrixDense
    */
  SymmetricMatrixDense()
      : size_(-1),
		    hessian_inverse_values_(nullptr),
		    hessian_values_(nullptr) {};

  //@}

  /** @name Destructor */
  //@{
  /**
    * Delete array
    */
  ~SymmetricMatrixDense();
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void addOptions(Options* options,
                  const Reporter* reporter);
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void setOptions(const Options* options,
                  const Reporter* reporter);
  //@}

  /** @name Initialize method */
  //@{
  /**
   * Initialize strategy
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void initialize(const Options* options,
                  Quantities* quantities,
                  const Reporter* reporter);
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get column of Hessian
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  void const columnHessian(int column_index,
                           Vector& column);
  /**
   * Get column of inverse Hessian
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  void const columnHessianInverse(int column_index,
                                  Vector& column);
  /**
    * Get element of Hessian
    * \param[in] row_index is row index number
    * \param[in] column_index is column index number
    * \return (row_index,column_index) element of matrix
    */
  double const elementHessian(int row_index,
                              int column_index);
  /**
    * Get element of inverse Hessian
    * \param[in] row_index is row index number
    * \param[in] column_index is column index number
    * \return (row_index,column_index) element of matrix
    */
  double const elementHessianInverse(int row_index,
                                     int column_index);
  /**
   * Get inner product of Hessian with vector
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  double innerProductHessian(const Vector& vector);
  /**
   * Get inner product of inverse Hessian with vector
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  double innerProductHessianInverse(const Vector& vector);
  /**
   * Get product of Hessian with vector
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  void matrixVectorProductHessian(const Vector& vector,
                                  Vector& product);
  /**
   * Get product of inverse Hessian with vector
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  void matrixVectorProductHessianInverse(const Vector& vector,
                                         Vector& product);
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "Dense"; };
  /**
    * Get number of rows
    * \return number of rows of the matrix
    */
  inline int const size() const { return size_; };
  /**
   * Get values (const) of Hessian
   * \return is pointer to array of SymmetricMatrixDense values
   */
  inline double* const valuesHessian() const { return hessian_values_; };
  /**
   * Get values (const) of inverse Hessian
   * \return is pointer to array of SymmetricMatrixDense values
   */
  inline double* const valuesHessianInverse() const { return hessian_inverse_values_; };
  /**
   * Get values (modifiable) of Hessian
   * \return is pointer to array of SymmetricMatrixDense values (to allow modification of array)
   */
  inline double* valuesModifiableHessian() { return hessian_values_; };
  /**
   * Get values (modifiable) of inverse Hessian
   * \return is pointer to array of SymmetricMatrixDense values (to allow modification of array)
   */
  inline double* valuesModifiableHessianInverse() { return hessian_inverse_values_; };
  //@}

  /** @name Modify methods */
  //@{
  /**
    * Set as diagonal matrix
    * \param[in] size is size of matrix to create
    * \param[in] value is value to set in diagonal elements (and set all else zero)
    */
  void setAsDiagonal(int size,
                     double value);
  /**
    * BFGS update
    * \param[in] s is reference to Vector representing iteration displacement
    * \param[in] y is reference to Vector representing gradient displacement
    */
  void updateBFGS(const Vector& s,
                  const Vector& y);
  //@}

  /** @name Print methods */
  //@{
  /**
    * Print array
    * \param[in] name is name of Symmetric Matrix to print
    * \param[in] reporter is pointer to Reporter object from NonOpt
    */
  void print(const Reporter* reporter,
             std::string name) const;
  //@}

 private:
  /** @name Default compiler generated methods
    * (Hidden to avoid implicit creation/calling.)
    */
  //@{
  /**
    * Copy constructor
    */
  SymmetricMatrixDense(const SymmetricMatrixDense&);
  /**
    * Overloaded equals operator
    */
  void operator=(const SymmetricMatrixDense&);
  //@}

  /** @name Private members */
  //@{
  int size_;                       /**< Number of rows and number of columns */
  int length_;                     /**< Number of rows *   number of columns */
  double* hessian_inverse_values_; /**< Double array */
  double* hessian_values_;         /**< Double array */
  //@}

  /** @name Indexing methods */
  //@{
  /**
    * Row index
    * \param[in] i is matrix array index
    * \return row corresponding to array index i
    */
  inline int const row_(int i) const { return i / size_; };
  /**
    * Column index
    * \param[in] i is matrix array index
    * \return column corresponding to array index i
    */
  inline int const col_(int i) const { return i % size_; };
  //@}

};  // end SymmetricMatrixDense

}  // namespace NonOpt

#endif /* __NONOPTSYMMETRICMATRIXDENSE_HPP__ */
