// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

/**
 *
 * Primal/dual QP solver for pair of the form
 *
 * (PRIMAL) max_d (max_i b_i + g_i'*d) + 0.5*d'*H*d
 *          s.t. ||d||_inf <= r
 *
 * and, with G = [g_1 ... g_m] and W = inv(H),
 *
 * (DUAL) min_(omega, gamma) 0.5*(G*omega + gamma)'*W*(G*omega + gamma) - b'*omega + r*||gamma||_1
 *        s.t.               sum(omega) = 1 and omega >= 0.
 *
 */

#ifndef NONOPTQPSOLVER_HPP_
#define NONOPTQPSOLVER_HPP_

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
   */
  virtual void addOptions(Options* options) = 0;
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   */
  virtual void setOptions(Options* options) = 0;
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
   * Get combination of vectors' infinity norm
   * \return "||G*omega||_inf"
   */
  virtual double combinationNormInf();
  /**
   * Get translated combination of vectors' infinity norm
   * \return "||G*omega + gamma||_inf"
   */
  virtual double combinationTranslatedNormInf();
  /**
   * Get translated combination of vectors' infinity norm
   * \return "||G*omega + gamma||_2^2"
   */
  virtual double combinationTranslatedNorm2Squared();
  /**
   * Get dual objective quadratic value
   * \return "(G*omega + gamma)'*W*(G*omega + gamma)"
   */
  virtual double dualObjectiveQuadraticValue();
  /**
   * Get dual solution
   * \param[out] omega is dual solution, omega part
   * \param[out] gamma is dual solution, gamma part
   */
  virtual void dualSolution(double omega[], double gamma[]);
  /**
   * Get dual solution, omega part
   * \param[out] omega is dual solution, omega part
   */
  virtual void dualSolutionOmega(double omega[]);
  /**
   * Get dual solution, omega part, length
   * \return dual solution, omega part, length
   */
  inline int dualSolutionOmegaLength() { return (int)vector_.size(); };
  /**
   * Get KKT error
   * \return solver's KKT error
   */
  virtual double KKTError() = 0;
  /**
   * Get KKT error full
   * \return full KKT error corresponding to dual solution
   */
  virtual double KKTErrorDual() = 0;
  /**
   * Get iteration count
   * \return number of iterations performed
   */
  virtual int numberOfIterations() = 0;
  /**
   * Get primal solution
   * \param[out] "d" (equal to "-W*(G*omega + gamma)" if solution is exact)
   */
  virtual void primalSolution(double d[]);
  /**
   * Get primal solution that is feasible
   * \param[out] feasible primal solution
   */
  virtual void primalSolutionFeasible(double d_feasible[]);
  /**
   * Get primal solution infinity norm
   * \return "||d||_inf"
   */
  virtual double primalSolutionNormInf();
  /**
   * Get primal solution 2-norm square
   * \return "||d||_2^2"
   */
  virtual double primalSolutionNorm2Squared();
  /**
   * Get feasible primal solution infinity norm
   * \return inf-norm of feasible primal solution
   */
  virtual double primalSolutionFeasibleNormInf();
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
  /**
   * Vector list length
   */
  inline int const vectorListLength() const { return (int)vector_list_.size(); };
  //@}

  /** @name Set method */
  //@{
  /**
   * Set inexact solution tolerance
   */
  virtual void setInexactSolutionTolerance(double tolerance) = 0;
  /**
   * Set "d" to zero
   */
  virtual void setPrimalSolutionToZero() = 0;
  /**
   * Set matrix
   * \param[in] matrix is pointer to SymmetricMatrix, for which "W" is the "Inverse"
   */
  virtual void setMatrix(const std::shared_ptr<SymmetricMatrix> matrix) = 0;
  /**
   * Set vector list
   * \param[in] vector_list is vector of pointers to Vectors to be set as QP "G" data
   */
  virtual void setVectorList(const std::vector<std::shared_ptr<Vector>> vector_list) = 0;
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
  virtual void addData(const std::vector<std::shared_ptr<Vector>> vector_list,
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
                       const Reporter* reporter,
                       Quantities* quantities) = 0;
  /**
   * Solve QP hot, after new data added, re-using previous solution, factorization, etc.
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void solveQPHot(const Options* options,
                          const Reporter* reporter,
                          Quantities* quantities) = 0;
  //@}

protected:
  /** @name Protected members
   * (Shared QP data/solution quantities common to all QPSolver implementations.)
   */
  //@{
  /**
   * QP data quantities
   */
  double scalar_;                                    /**< "r" */
  std::shared_ptr<SymmetricMatrix> matrix_;          /**< "W"/"H", matrix (or its inverse) defining the quadratic term */
  std::vector<std::shared_ptr<Vector>> vector_list_; /**< "G" */
  std::vector<double> vector_;                       /**< "b" */
  int gamma_length_;
  /**
   * Solution quantities
   */
  double primal_solution_projection_scalar_;
  Vector combination_;
  Vector combination_translated_;
  Vector gamma_;
  Vector omega_;
  Vector primal_solution_;
  Vector primal_solution_feasible_;
  //@}

  /** @name Protected methods
   * (Internal helpers shared by all QPSolver implementations.)
   */
  //@{
  /**
   * Sanity check
   */
  virtual bool checkQuantityCompatibility();
  /**
   * Dual objective quadratic value, scaled to correspond to feasible "d"
   */
  virtual double dualObjectiveQuadraticValueScaled();
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  QPSolver(const QPSolver&) = delete;
  /**
   * Overloaded equals operator
   */
  void operator=(const QPSolver&) = delete;
  //@}

  /** @name Private members */
  //@{
  QP_Status status_; /**< Termination status */
  //@}

}; // end QPSolver

} // namespace NonOpt

#endif /* NONOPTQPSOLVER_HPP_ */
