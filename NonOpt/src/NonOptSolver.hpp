// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSOLVER_HPP__
#define __NONOPTSOLVER_HPP__

#include <ctime>
#include <memory>

#include "NonOptEnumerations.hpp"
#include "NonOptOptions.hpp"
#include "NonOptProblem.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptStrategies.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Problem;
class Quantities;
class Reporter;
class Strategies;

/**
 * NonOptSolver class
 */
class NonOptSolver
{

public:
  /** @name Constructors */
  //@{
  /**
   * Construct NonOptSolver
   */
  NonOptSolver();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Delete NonOptSolver
   */
  ~NonOptSolver();

  /** @name Get methods */
  //@{
  /**
   * Get function evaluation counter
   * \return function evaluations so far
   */
  inline int const functionEvaluations() const { return quantities_.functionCounter(); };
  /**
   * Get gradient evaluation counter
   * \return gradient evaluations so far
   */
  inline int const gradientEvaluations() const { return quantities_.gradientCounter(); };
  /**
   * Get iteration counter
   * \return iterations performed so far
   */
  inline int const iterations() const { return quantities_.iterationCounter(); };
  /**
   * Get number of variables
   * \return number of variables
   */
  inline int const numberOfVariables() const { return quantities_.numberOfVariables(); };
  /**
   * Get inner iteration counter
   * \return total inner iterations performed so far
   */
  inline int const totalInnerIterations() const { return quantities_.totalInnerIterationCounter(); };
  /**
   * Get inner iteration counter
   * \return total inner iterations performed so far
   */
  inline int const totalQPIterations() const { return quantities_.totalQPIterationCounter(); };
  /**
   * Get objective value
   * \return objective value of current iterate
   */
  inline double const objective() { return quantities_.currentIterate()->objectiveUnscaled(); };
  /**
   * Get stationarity radius
   * \return current stationarity radius
   */
  inline double const stationarityRadius() const { return quantities_.stationarityRadius(); };
  /**
   * Get time in evaluations
   * \return seconds between start and end time
   */
  inline double const time() const { return (quantities_.endTime() - quantities_.startTime()) / (double)CLOCKS_PER_SEC; };
  /**
   * Get time in evaluations
   * \return seconds to perform problem function evaluations
   */
  inline double const timeEvaluations() const { return quantities_.evaluationTime() / (double)CLOCKS_PER_SEC; };
  /**
   * Get time in NonOpt
   * \return seconds between start and end time not including problem function evaluation time
   */
  inline double const timeNonOpt() const { return (quantities_.endTime() - quantities_.startTime() - quantities_.evaluationTime()) / (double)CLOCKS_PER_SEC; };
  /**
   * Get status
   * \return current status of algorithm
   */
  inline NonOpt_Status const status() const { return status_; };
  /**
   * Get options
   * \return pointer to Options object
   */
  inline Options* options() { return &options_; };
  /**
   * Get reporter
   * \return pointer to Reporter object
   */
  inline Reporter* reporter() { return &reporter_; };
  /**
   * Get solution
   * \param[out] vector is the current iterate
   */
  void solution(double vector[]);
  //@}

  /** @name Set method */
  //@{
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(NonOpt_Status status) { status_ = status; };
  //@}

  /** @name Optimize method */
  //@{
  /**
   * Optimize
   * \param[in] problem is a pointer to a Problem object
   */
  void optimize(const std::shared_ptr<Problem> problem);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  NonOptSolver(const NonOptSolver&);
  /**
   * Overloaded equals operator
   */
  void operator=(const NonOptSolver&);
  //@}

  /** @name Private members */
  //@{
  NonOpt_Status status_;
  //@}

  /** @name Private members, objects */
  //@{
  Options options_;
  Quantities quantities_;
  Reporter reporter_;
  Strategies strategies_;
  //@}

  /** @name Private members, options */
  //@{
  int print_level_;
  int print_level_file_;
  int qp_print_level_;
  int qp_print_level_file_;
  std::string print_file_name_;
  std::string qp_print_file_name_;
  //@}

  /** @name Private methods */
  //@{
  void addOptions();
  void initialize(const std::shared_ptr<Problem> problem);
  void printFooter();
  void printHeader();
  void printIterationHeader();
  void setOptions();
  //@}

}; // end NonOptSolver

} // namespace NonOpt

#endif /* __NONOPTSOLVER_HPP__ */
