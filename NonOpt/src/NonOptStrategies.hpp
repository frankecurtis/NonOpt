// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSTRATEGIES_HPP__
#define __NONOPTSTRATEGIES_HPP__

#include <memory>
#include <string>

#include "NonOptApproximateHessianUpdate.hpp"
#include "NonOptDerivativeChecker.hpp"
#include "NonOptDirectionComputation.hpp"
#include "NonOptLineSearch.hpp"
#include "NonOptPointSetUpdate.hpp"
#include "NonOptQPSolver.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptTermination.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class ApproximateHessianUpdate;
class DerivativeChecker;
class DirectionComputation;
class LineSearch;
class PointSetUpdate;
class QPSolver;
class SymmetricMatrix;
class Termination;

/**
 * Strategies class
 */
class Strategies
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Strategies(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~Strategies(){};
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
   * Initialize strategies
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
   * Get pointer to ApproximateHessianUpdate
   * \return pointer to ApproximateHessianUpdate object
   */
  inline std::shared_ptr<ApproximateHessianUpdate> approximateHessianUpdate() { return approximate_hessian_update_; }
  /**
   * Get pointer to DerivativeChecker
   * \return pointer to DerivativeChecker object
   */
  inline std::shared_ptr<DerivativeChecker> derivativeChecker() { return derivative_checker_; }
  /**
   * Get pointer to DirectionComputation
   * \return pointer to DirectionComputation object
   */
  inline std::shared_ptr<DirectionComputation> directionComputation() { return direction_computation_; }
  /**
   * Get pointer to LineSearch
   * \return pointer to LineSearch object
   */
  inline std::shared_ptr<LineSearch> lineSearch() { return line_search_; }
  /**
   * Get pointer to PointSetUpdate
   * \return pointer to PointSetUpdate object
   */
  inline std::shared_ptr<PointSetUpdate> pointSetUpdate() { return point_set_update_; }
  /**
   * Get pointer to QPSolver
   * \return pointer QPSolver object
   */
  inline std::shared_ptr<QPSolver> qpSolver(bool small) { return (small) ? qp_solver_small_scale_ : qp_solver_large_scale_; }
  /**
   * Get pointer to QPSolverTermination
   * \return pointer QPSolver object
   */
  inline std::shared_ptr<QPSolver> qpSolverTermination(bool small) { return (small) ? qp_solver_small_scale_termination_ : qp_solver_large_scale_termination_; }
  /**
   * Get pointer to SymmetricMatrix
   * \return pointer to SymmetricMatrix object
   */
  inline std::shared_ptr<SymmetricMatrix> symmetricMatrix() { return symmetric_matrix_; }
  /**
   * Get pointer to SymmetricMatrixTermination
   * \return pointer to SymmetricMatrix object
   */
  inline std::shared_ptr<SymmetricMatrix> symmetricMatrixTermination() { return symmetric_matrix_termination_; }
  /**
   * Get pointer to Termination
   * \return pointer to Termination object
   */
  inline std::shared_ptr<Termination> termination() { return termination_; }
  /**
   * Get iteration header
   * \return iteration header as string
   */
  inline std::string iterationHeader() const { return iteration_header_; }
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set iteration header
   */
  void setIterationHeader();
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print header
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void printHeader(const Reporter* reporter);
  /**
   * Print footer
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void printFooter(const Reporter* reporter);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Strategies(const Strategies&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Strategies&);
  //@}

  /** @name * Private members */
  //@{
  std::shared_ptr<ApproximateHessianUpdate> approximate_hessian_update_;
  std::shared_ptr<DerivativeChecker> derivative_checker_;
  std::shared_ptr<DirectionComputation> direction_computation_;
  std::shared_ptr<LineSearch> line_search_;
  std::shared_ptr<PointSetUpdate> point_set_update_;
  std::shared_ptr<QPSolver> qp_solver_small_scale_;
  std::shared_ptr<QPSolver> qp_solver_small_scale_termination_;
  std::shared_ptr<QPSolver> qp_solver_large_scale_;
  std::shared_ptr<QPSolver> qp_solver_large_scale_termination_;
  std::shared_ptr<SymmetricMatrix> symmetric_matrix_;
  std::shared_ptr<SymmetricMatrix> symmetric_matrix_termination_;
  std::shared_ptr<Termination> termination_;
  std::string iteration_header_;
  //@}

}; // end Strategies

} // namespace NonOpt

#endif /* __NONOPTSTRATEGIES_HPP__ */
