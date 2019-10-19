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
		hessian_values_(nullptr){};

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
   * Get column
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  void const column(int column_index,
                    Vector& column);
  /**
   * Get column
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  void const column_Hessian(int column_index,
                    Vector& column);
  /**
   * Get element (const)
   * \param[in] row_index is row index number
   * \param[in] column_index is column index number
   * \return (row_index,column_index) element of matrix
   */
  double const element(int row_index,
                       int column_index);
  /**
   * Get element (const)
   * \param[in] row_index is row index number
   * \param[in] column_index is column index number
   * \return (row_index,column_index) element of matrix
   */
  double const element_Hessian(int row_index,
                       int column_index);
  /**
   * Get inner product with vector for Inverse Hessian approximation
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  double innerProduct_HessianInverse(const Vector& vector);
  /**
   * Get inner product with vector for Hessian approximation
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  double innerProduct_Hessian(const Vector& vector);
  /**
   * Get product with vector for Inverse Hessian approximation
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  void matrixVectorProduct_HessianInverse(const Vector& vector,
                           Vector& product);
  /**
   * Get product with vector for Hessian approximation
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  void matrixVectorProduct_Hessian(const Vector& vector,
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
   * Get values (const)
   * \return is pointer to array of SymmetricMatrixDense values
   */
  inline double* const values() const { return hessian_inverse_values_; };
  /**
   * Get values (modifiable)
   * \return is pointer to array of SymmetricMatrixDense values (to allow modification of array)
   */
  inline double* valuesModifiable() { return hessian_inverse_values_; };
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
  int size_;       /**< Number of rows and number of columns */
  int length_;     /**< Number of rows *   number of columns */
  double* hessian_inverse_values_; /**< Double array   */
  double* hessian_values_; /**< Double array
                      */

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
