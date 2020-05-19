// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSYMMETRICMATRIXLIMITEDMEMORY_HPP__
#define __NONOPTSYMMETRICMATRIXLIMITEDMEMORY_HPP__

#include "NonOptSymmetricMatrix.hpp"

namespace NonOpt
{

/**
  * SymmetricMatrixLimitedMemory class
  */
class SymmetricMatrixLimitedMemory : public SymmetricMatrix
{

public:
  /** @name Constructors */
  //@{
  /**
    * Construct SymmetricMatrixLimitedMemory
    */
  SymmetricMatrixLimitedMemory()
    : size_(-1),
      history_(-1),
      initial_diagonal_value_(1.0),
      values_(nullptr)
  {
    s_.clear();
    y_.clear();
    rho_.clear();
    computed_columns_.clear();
    computed_columns_of_inverse_.clear();
    computed_column_indices_.clear();
    computed_column_indices_of_inverse_.clear();
  };
  //@}

  /** @name Destructor */
  //@{
  /**
    * Delete array
    */
  ~SymmetricMatrixLimitedMemory();
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
   * \param[in options is pointer to Options object from NonOpt
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
  std::string name() { return "LimitedMemory"; };
  /**
    * Get number of rows
    * \return number of rows of the matrix
    */
  inline int const size() const { return size_; };
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
    * \param[in] reporter is pointer to Reporter object from NonOpt
    * \param[in] name is name of Symmetric Matrix to print
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
  SymmetricMatrixLimitedMemory(const SymmetricMatrixLimitedMemory&);
  /**
    * Overloaded equals operator
    */
  void operator=(const SymmetricMatrixLimitedMemory&);
  //@}

  /** @name Private members */
  //@{
  int size_;                                                         /**< Number of rows and number of columns */
  int history_;                                                      /**< Limited memory history length */
  double initial_diagonal_value_;                                    /**< Diagonal value of "initial" matrix */
  double* values_;                                                   /**< Double array, values of matrix (TEMPORARY) */
  std::vector<std::shared_ptr<Vector>> s_;                           /**< Vector vector, "s" values */
  std::vector<std::shared_ptr<Vector>> y_;                           /**< Vector vector, "y" values */
  std::vector<double> rho_;                                          /**< Double vector, "rho" values */
  std::vector<std::shared_ptr<Vector>> computed_columns_;            /**< Vector vector, computed columns */
  std::vector<std::shared_ptr<Vector>> computed_columns_of_inverse_; /**< Vector vector, computed columns of inverse */
  std::vector<int> computed_column_indices_;                         /**< Integer vector, computed column indices */
  std::vector<int> computed_column_indices_of_inverse_;              /**< Integer vector, computed column indices of inverse */
  //@}

}; // end SymmetricMatrixLimitedMemory

} // namespace NonOpt

#endif /* __NONOPTSYMMETRICMATRIXLIMITEDMEMORY_HPP__ */
