// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

/**
 *
 * QP solver for problems of the form
 *
 * min_(omega, gamma) 0.5*(G*omega + gamma)'*W*(G*omega + gamma) - b'*omega + r*||gamma||_1
 * s.t.               sum(omega) = 1 and omega >= 0
 *
 */

#ifndef __NONOPTQPSOLVER_HPP__
#define __NONOPTQPSOLVER_HPP__

#include <memory>
#include <string>
#include <vector>

#include "NonOptEnumerations.hpp"
#include "NonOptOptions.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptStrategy.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptVector.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Quantities;
class Options;
class Reporter;
class Strategy;
class SymmetricMatrix;
class Vector;

/**
 * QPSolver class
 */
class QPSolver : public Strategy
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  QPSolver(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~QPSolver(){};

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

  /** @name Initialize method */
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
   * Get combination of vectors' infinity norm
   * \return "||G*omega||_inf"
   */
  virtual double combinationNormInf() = 0;
  /**
   * Get translated combination of vectors' infinity norm
   * \return "||G*omega + gamma||_inf"
   */
  virtual double combinationTranslatedNormInf() = 0;
  /**
   * Get dual step
   * \param[out] vector is dual step given by "-W*(G*omega + gamma)"
   */
  virtual void dualStep(double vector[]) = 0;
  /**
   * Get dual step's infinity norm
   * \return "||W*(G*omega + gamma)||_inf"
   */
  virtual double dualStepNormInf() = 0;
  /**
   * Get objective quadratic value
   * \return (G*omega + gamma)'*W*(G*omega + gamma)
   */
  virtual double objectiveQuadraticValue() = 0;
  /**
   * Get gamma
   * \param[out] vector is "gamma" solution value
   */
  virtual void gamma(double vector[]) = 0;
  /**
   * Get KKT error
   * \return solver's KKT error corresponding to current solution (ignores certain conditions assumed to be satisfied during solve)
   */
  virtual double KKTError() = 0;
  /**
   * Get KKT error full
   * \return full KKT error corresponding to current solution
   */
  virtual double KKTErrorFull() = 0;
  /**
   * Get iteration count
   * \return number of iterations performed
   */
  virtual int numberOfIterations() = 0;
  /**
   * Get omega
   * \param[out] vector is "omega" solution value
   */
  virtual void omega(double vector[]) = 0;
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  /**
   * Get status
   * \return current status of direction computation
   */
  inline QP_Status status() { return status_; };
  //@}

  /** @name Set method */
  //@{
  /**
   * Set dual step to zero
   */
  virtual void setDualStepToZero() = 0;
  /**
   * Set matrix
   * \param[in] matrix is pointer to SymmetricMatrix to be set as QP "W" data
   */
  virtual void setMatrix(const std::shared_ptr<SymmetricMatrix> matrix) = 0;
  /**
   * Set vector list
   * \param[in] vector_list is vector of pointers to Vectors to be set as QP "G" data
   */
  virtual void setVectorList(const std::vector<std::shared_ptr<Vector> > vector_list) = 0;
  /**
   * Set vector
   * \param[in] vector is vector of double values to be set as QP "b" data
   */
  virtual void setVector(const std::vector<double> vector) = 0;
  /**
   * Set scalar
   * \param[in] scalar is double value to be set as QP "r" data
   */
  virtual void setScalar(double scalar) = 0;
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(QP_Status status) { status_ = status; };
  //@}

  /** @name Add data methods */
  //@{
  /**
   * Add data
   * \param[in] vector_list is vector of pointers to Vectors to be added to QP "G" data
   * \param[in] vector is vector of double values to be added to QP "b" data
   */
  virtual void addData(const std::vector<std::shared_ptr<Vector> > vector_list,
                       const std::vector<double> vector) = 0;
  //@}

  /** @name QP solve methods */
  //@{
  /**
   * Solve QP
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void solveQP(const Options* options,
                       const Reporter* reporter) = 0;
  /**
   * Solve QP hot, after new data added, re-using previous solution, factorization, etc.
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void solveQPHot(const Options* options,
                          const Reporter* reporter) = 0;
  //@}

 private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  QPSolver(const QPSolver&);
  /**
   * Overloaded equals operator
   */
  void operator=(const QPSolver&);
  //@}

  /** @name Private members */
  //@{
  QP_Status status_; /**< Termination status */
  //@}

};  // end QPSolver

}  // namespace NonOpt

#endif /* __NONOPTQPSOLVER_HPP__ */
