// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTITERATIONQUANTITIES_HPP__
#define __NONOPTITERATIONQUANTITIES_HPP__

#include <ctime>
#include <memory>
#include <string>
#include <vector>

#include "NonOptOptions.hpp"
#include "NonOptPoint.hpp"
#include "NonOptProblem.hpp"
#include "NonOptReporter.hpp"
#include "NonOptVector.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Point;
class Problem;
class Reporter;
class Vector;

/**
 * Quantities class
 */
class Quantities
{

public:
  /** @name Constructors */
  //@{
  /**
   * Declare Quantities
   */
  Quantities();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Delete data
   */
  ~Quantities();
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

  /** @name Initialization methods */
  //@{
  /**
   * Initialize quantities
   * \param[in] problem is pointer to Problem object
   */
  bool initialize(const std::shared_ptr<Problem> problem);
  /**
   * Initialize radii
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void initializeRadii(const Options* options,
                       const Reporter* reporter);
  /**
   * Initialize inexact termination factor
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  inline void initializeInexactTerminationFactor(const Options* options,
                                                 const Reporter* reporter)
  {
    inexact_termination_factor_ = inexact_termination_factor_initial_;
  };
  //@}

  /** @name Get methods */
  //@{
  /**
   * Start time
   * \return start time that was set
   */
  inline clock_t const startTime() const { return start_time_; };
  /**
   * \return end time that was set
   */
  inline clock_t const endTime() const { return end_time_; };
  /**
   * Evaluation time
   * \return problem function evaluation time that was set
   */
  inline clock_t const evaluationTime() const { return evaluation_time_; };
  /**
   * Get inexact termination factor
   * \return current inexact termination factor
   */
  inline double const inexactTerminationFactor() const { return inexact_termination_factor_; };
  /**
   * Get scaling threshold
   * \return scaling threshold
   */
  inline double const scalingThreshold() const { return scaling_threshold_; };
  /**
   * Get stationarity radius
   * \return current stationarity radius
   */
  inline double const stationarityRadius() const { return stationarity_radius_; };
  /**
   * Get stepsize
   * \return current stepsize
   */
  inline double const stepsize() const { return stepsize_; };
  /**
   * Get trust region radius
   * \return current trust region radius
   */
  inline double const trustRegionRadius() const { return trust_region_radius_; };
  /**
   * Function evaluation counter
   * \return function evaluations performed so far
   */
  inline int const functionCounter() const { return function_counter_; };
  /**
   * Function evaluation limit
   * \return function evaluation limit
   */
  inline int const functionEvaluationLimit() const { return function_evaluation_limit_; };
  /**
   * Gradient evaluation counter
   * \return gradient evaluations performed so far
   */
  inline int const gradientCounter() const { return gradient_counter_; };
  /**
   * Gradient evaluation limit
   * \return gradient evaluation limit
   */
  inline int const gradientEvaluationLimit() const { return gradient_evaluation_limit_; };
  /**
   * Iteration counter
   * \return iterations performed so far
   */
  inline int const iterationCounter() const { return iteration_counter_; };
  /**
   * Inner iteration counter
   * \return inner iterations performed so far (during current iteration)
   */
  inline int const innerIterationCounter() const { return inner_iteration_counter_; };
  /**
   * QP iteration counter
   * \return QP iterations performed so far (during current iteration)
   */
  inline int const QPIterationCounter() const { return qp_iteration_counter_; };
  /**
   * Inner iteration counter
   * \return total inner iterations performed so far (over all iterations so far)
   */
  inline int const totalInnerIterationCounter() const { return total_inner_iteration_counter_; };
  /**
   * QP iteration counter
   * \return total qp iterations performed so far (over all iterations so far)
   */
  inline int const totalQPIterationCounter() const { return total_qp_iteration_counter_; };
  /**
   * Get problem size
   * \return number of variables
   */
  inline int const numberOfVariables() const { return number_of_variables_; };
  /**
   * Get pointer to current iterate
   * \return pointer to Point representing current iterate
   */
  inline std::shared_ptr<Point> currentIterate() { return current_iterate_; };
  /**
   * Get pointer to trial iterate
   * \return pointer to Point representing trial iterate
   */
  inline std::shared_ptr<Point> trialIterate() { return trial_iterate_; };
  /**
   * Get direction
   * \return pointer to Vector representing search direction
   */
  inline std::shared_ptr<Vector> direction() { return direction_; };
  /**
   * Get point set
   * \return pointer to vector of pointers to Points representing current point set
   */
  inline std::shared_ptr<std::vector<std::shared_ptr<Point>>> pointSet() { return point_set_; };
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set current iterate pointer
   * \param[in] iterate is pointer to Point to represent current iterate
   */
  inline void setCurrentIterate(const std::shared_ptr<Point> iterate) { current_iterate_ = iterate; };
  /**
   * Set trial iterate pointer
   * \param[in] trial_iterate is pointer to Point to represent trial iterate
   */
  inline void setTrialIterate(const std::shared_ptr<Point> trial_iterate) { trial_iterate_ = trial_iterate; };
  /**
   * Set trial iterate pointer to current iterate pointer
   */
  inline void setTrialIterateToCurrentIterate() { trial_iterate_ = current_iterate_; };
  /**
   * Set stepsize
   * \param[in] stepsize is new value to represent stepsize
   */
  inline void setStepsize(double stepsize) { stepsize_ = stepsize; };
  /**
   * Update radii
   * \param[in] stationarity_tolerance is stationarity tolerance for reference
   */
  void updateRadii(double stationarity_tolerance);
  /**
   * Update inexact termination factor
   */
  inline void updateInexactTerminationFactor()
  {
    inexact_termination_factor_ = inexact_termination_factor_ * inexact_termination_update_factor_;
  };
  //@}

  /** @name Increment methods */
  //@{
  /**
   * Increment evaluation time
   * \param[in] evaluation_time is amount to add to total problem function evaluation time
   */
  inline void incrementEvaluationTime(clock_t evaluation_time) { evaluation_time_ += evaluation_time; };
  /**
   * Increment function evaluation counter
   */
  inline void incrementFunctionCounter() { function_counter_++; };
  /**
   * Increment gradient evaluation counter
   */
  inline void incrementGradientCounter() { gradient_counter_++; };
  /**
   * Increment iteration counter
   */
  inline void incrementIterationCounter() { iteration_counter_++; };
  /**
   * Increments inner iteration counter by given amount
   * \param[in] amount is amount to increment inner iteration counter
   */
  inline void incrementInnerIterationCounter(int amount) { inner_iteration_counter_ += amount; };
  /**
   * Increments QP iteration counter by given amount
   * \param[in] amount is amount to increment QP iteration counter
   */
  inline void incrementQPIterationCounter(int amount) { qp_iteration_counter_ += amount; };
  /**
   * Increment total inner iteration counter
   */
  inline void incrementTotalInnerIterationCounter() { total_inner_iteration_counter_ += inner_iteration_counter_; };
  /**
   * Increment total qp iteration counter
   */
  inline void incrementTotalQPIterationCounter() { total_qp_iteration_counter_ += qp_iteration_counter_; };
  //@}

  /** @name Reset methods */
  //@{
  /**
   * Reset inner iteration counter
   */
  inline void resetInnerIterationCounter() { inner_iteration_counter_ = 0; };
  /**
   * Reset QP iteration counter
   */
  inline void resetQPIterationCounter() { qp_iteration_counter_ = 0; };
  //@}

  /** @name Print methods */
  //@{
  /**
   * Get iteration header string
   * \return string of header values
   */
  std::string iterationHeader();
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues();
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print header
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void printHeader(const Reporter* reporter);
  /**
   * Print iteration values
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void printIterationValues(const Reporter* reporter);
  /**
   * Print footer
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void printFooter(const Reporter* reporter);
  //@}

  /** @name Finalization method */
  //@{
  /**
   * Finalize
   */
  void finalize();
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Quantities(const Quantities&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Quantities&);
  //@}

  /** @name Private members */
  //@{
  clock_t end_time_;
  clock_t evaluation_time_;
  clock_t start_time_;
  double inexact_termination_factor_;
  double stationarity_radius_;
  double stepsize_;
  double trust_region_radius_;
  int function_counter_;
  int gradient_counter_;
  int iteration_counter_;
  int inner_iteration_counter_;
  int number_of_variables_;
  int qp_iteration_counter_;
  int total_inner_iteration_counter_;
  int total_qp_iteration_counter_;
  std::shared_ptr<Point> current_iterate_;
  std::shared_ptr<Point> trial_iterate_;
  std::shared_ptr<Vector> direction_;
  std::shared_ptr<std::vector<std::shared_ptr<Point>>> point_set_;
  //@}

  /** @name Private members (options) */
  //@{
  double inexact_termination_factor_initial_;
  double inexact_termination_update_factor_;
  double scaling_threshold_;
  double stationarity_radius_initialization_factor_;
  double stationarity_radius_initialization_minimum_;
  double stationarity_radius_update_factor_;
  double trust_region_radius_initialization_factor_;
  double trust_region_radius_initialization_minimum_;
  double trust_region_radius_update_factor_;
  int function_evaluation_limit_;
  int gradient_evaluation_limit_;
  //@}

}; // end Quantities

} // namespace NonOpt

#endif /* __NONOPTITERATIONQUANTITIES_HPP__ */
