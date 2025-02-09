// Copyright (C) 2025 Frank E. Curtis
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
      values_(nullptr),
      values_of_inverse_(nullptr){};
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
   */
  void addOptions(Options* options);
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   */
  void setOptions(Options* options);
  //@}

  /** @name Initialization method */
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
   * Get column of symmetric matrix
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  void const column(int column_index,
                    Vector& column);
  /**
   * Get column of symmetric matrix inverse
   * \param[in] index is index of column to return
   * \param[out] column is Vector to store column values
   */
  void const columnOfInverse(int column_index,
                             Vector& column);
  /**
   * Get element of symmetric matrix
   * \param[in] row_index is row index number
   * \param[in] column_index is column index number
   * \return (row_index,column_index) element of matrix
   */
  double const element(int row_index,
                       int column_index);
  /**
   * Get element of symmetric matrix inverse
   * \param[in] row_index is row index number
   * \param[in] column_index is column index number
   * \return (row_index,column_index) element of matrix
   */
  double const elementOfInverse(int row_index,
                                int column_index);
  /**
   * Get inner product of symmetric matrix with vector
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  double innerProduct(const Vector& vector);
  /**
   * Get inner product of symmetric matrix inverse with vector
   * \param[in] vector is reference to a Vector
   * \return inner product of vector with this vector
   */
  double innerProductOfInverse(const Vector& vector);
  /**
   * Get product of symmetric matrix with vector
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  void matrixVectorProduct(const Vector& vector,
                           Vector& product);
  /**
   * Get product of symmetric matrix inverse with vector
   * \param[in] vector is reference to a Vector
   * \param[out] product is Vector to store product values
   */
  void matrixVectorProductOfInverse(const Vector& vector,
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
   * Get values (const) of symmetric matrix
   * \return is pointer to array of values
   */
  inline double* const values() const { return values_; };
  /**
   * Get values (const) of symmetric matrix inverse
   * \return is pointer to array of values
   */
  inline double* const valuesOfInverse() const { return values_of_inverse_; };
  /**
   * Get values (modifiable) of symmetric matrix
   * \return is pointer to array of values (to allow modification of array)
   */
  inline double* valuesModifiable() { return values_; };
  /**
   * Get values (modifiable) of symmetric matrix inverse
   * \return is pointer to array of values (to allow modification of array)
   */
  inline double* valuesOfInverseModifiable() { return values_of_inverse_; };
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
   * Set element
   * \param[in] row_index is row index number
   * \param[in] column_index is column index number
   * \param[in] value is value to set (row_index,column_index) element
   */
  void setElement(int row_index,
                  int column_index,
                  double value);                  
  /**
   * Update approximation
   * \param[in] s is reference to Vector representing iteration displacement
   * \param[in] y is reference to Vector representing gradient displacement
   */
  void update(const Vector& s,
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
  int size_;                  /**< Number of rows and number of columns */
  int length_;                /**< Number of rows *   number of columns */
  double* values_;            /**< Double array */
  double* values_of_inverse_; /**< Double array */
  //@}

  /** @name Private methods */
  //@{
  void updateBFGS(const Vector& s,
                  const Vector& y);
  void updateDFP(const Vector& s,
                 const Vector& y);
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

}; // end SymmetricMatrixDense

} // namespace NonOpt

#endif /* __NONOPTSYMMETRICMATRIXDENSE_HPP__ */
