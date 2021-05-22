// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptDefinitions.hpp"
#include "NonOptQuantities.hpp"

namespace NonOpt
{

// Constructor
Quantities::Quantities()
  : direction_computation_time_(0),
    evaluation_time_(0),
    line_search_time_(0),
    inexact_termination_factor_(0.0),
    iterate_norm_initial_(0.0),
    stationarity_radius_(0.0),
    stepsize_(0.0),
    trust_region_radius_(0.0),
    function_counter_(0),
    gradient_counter_(0),
    iteration_counter_(0),
    inner_iteration_counter_(0),
    number_of_variables_(0),
    qp_iteration_counter_(0),
    total_inner_iteration_counter_(0),
    total_qp_iteration_counter_(0),
    approximate_hessian_initial_scaling_(false),
    evaluate_function_with_gradient_(false),
    cpu_time_limit_(NONOPT_DOUBLE_INFINITY),
    inexact_termination_factor_initial_(1.0),
    inexact_termination_update_factor_(1.0),
    inexact_termination_update_stepsize_threshold_(1.0),
    iterate_norm_tolerance_(1.0),
    scaling_threshold_(1.0),
    stationarity_radius_initialization_factor_(1.0),
    stationarity_radius_initialization_minimum_(1.0),
    stationarity_radius_update_factor_(1.0),
    stationarity_tolerance_(0.0),
    trust_region_radius_initialization_factor_(1.0),
    trust_region_radius_initialization_minimum_(1.0),
    trust_region_radius_update_factor_(1.0),
    function_evaluation_limit_(10),
    gradient_evaluation_limit_(10),
    iteration_limit_(1)
{
  start_time_ = clock();
  end_time_ = start_time_;
  current_iterate_.reset();
  trial_iterate_.reset();
  direction_.reset();
  direction_termination_.reset();
  point_set_.reset();
}

// Destructor
Quantities::~Quantities(){}

// Add options
void Quantities::addOptions(Options* options)
{

  // Add bool options
  options->addBoolOption("approximate_hessian_initial_scaling",
                         false,
                         "Indicator of whether to scale initial matrix for approximate Hessian.\n"
                         "Default     : false");
  options->addBoolOption("evaluate_function_with_gradient",
                         false,
                         "Determines whether to evaluate function and gradient\n"
                         "              at the same time (or separately).\n"
                         "Default     : false");

  // Add double options
  options->addDoubleOption("cpu_time_limit",
                           1e+04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Limit on the number of CPU seconds.  This limit is only checked\n"
                           "              at the beginning of an iteration, so the true CPU time limit\n"
                           "              also depends on the time required to a complete an iteration.\n"
                           "Default     : 1e+04");
  options->addDoubleOption("inexact_termination_factor_initial",
                           sqrt(2.0) - 1.0,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Initial inexact termination factor.  Factor by which norm of\n"
                           "              inexact subproblem solution needs to be within norm of true\n"
                           "              (unknown) projection of origin onto convex hull of gradients.\n"
                           "Default     : sqrt(2.0)-1.0");
  options->addDoubleOption("inexact_termination_update_factor",
                           0.9999,
                           0.0,
                           1.0,
                           "Factor for updating the inexact termination factor.  If the\n"
                           "              conditions for updating the inexact termination factor are met,\n"
                           "              then the inexact termination factor is multiplied by this fraction.\n"
                           "Default     : 0.9999");
  options->addDoubleOption("inexact_termination_update_stepsize_threshold",
                           1e-10,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Stepsize threshold for determining whether to reduce the\n"
                           "              inexact termination factor.\n"
                           "Default     : 1e-10");
  options->addDoubleOption("iterate_norm_tolerance",
                           1e+20,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining divergence of the algorithm iterates.\n"
                           "              If the norm of an iterate is larger than this tolerance times\n"
                           "              the maximum of 1.0 and the norm of the initial iterate, then\n"
                           "              the algorithm terminates with a message of divergence.\n"
                           "Default     : 1e+20");
  options->addDoubleOption("scaling_threshold",
                           1e+02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Threshold for determining objective scaling.  If norm of gradient\n"
                           "              at the initial point is greater than this value, then the objective\n"
                           "              is scaled so that the initial gradient norm is at this value.\n"
                           "Default     : 1e+02");
  options->addDoubleOption("stationarity_radius_initialization_factor",
                           1e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for initializing the stationarity radius.  The initial\n"
                           "              stationarity radius is the maximum of this value times the\n"
                           "              inf-norm of the gradient at the initial point and\n"
                           "              stationarity_radius_initialization_minimum.\n"
                           "Default     : 1e-01");
  options->addDoubleOption("stationarity_radius_initialization_minimum",
                           1e-02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Minimum initial value for stationarity radius.\n"
                           "Default     : 1e-02");
  options->addDoubleOption("stationarity_radius_update_factor",
                           1e-01,
                           0.0,
                           1.0,
                           "Factor for updating the stationarity radius.  If the conditions\n"
                           "              for updating the stationarity and trust region radii are met,\n"
                           "              then the stationarity radius is multiplied by this fraction.\n"
                           "Default     : 1e-01");
  options->addDoubleOption("stationarity_tolerance",
                           1e-04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining stationarity.  If the stationarity\n"
                           "              radius falls below this tolerance and other factors determined\n"
                           "              by the termination strategy are met, then the algorithm\n"
                           "              terminates with a message of stationarity.\n"
                           "Default     : 1e-04");
  options->addDoubleOption("trust_region_radius_initialization_factor",
                           1e+04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for initializing the trust region radius.  The initial\n"
                           "              trust region radius is the maximum of this value times the\n"
                           "              inf-norm of the gradient at the initial point and\n"
                           "              trust_region_radius_initialization_minimum.\n"
                           "Default     : 1e+04");
  options->addDoubleOption("trust_region_radius_initialization_minimum",
                           1e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Minimum initial value for trust region radius.\n"
                           "Default     : 1e-01");
  options->addDoubleOption("trust_region_radius_update_factor",
                           1e-01,
                           0.0,
                           1.0,
                           "Factor for updating the trust region radius.  If the conditions\n"
                           "              for updating the stationarity and trust region radii are met,\n"
                           "              then the trust region radius is multiplied by this fraction.\n"
                           "Default     : 1e-01");

  // Add integer options
  options->addIntegerOption("function_evaluation_limit",
                            1e+05,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of function evaluations performed.\n"
                            "Default     : 1e+05");
  options->addIntegerOption("gradient_evaluation_limit",
                            1e+05,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of gradient evaluations performed.\n"
                            "Default     : 1e+05");
  options->addIntegerOption("iteration_limit",
                            1e+04,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of iterations that will be performed.\n"
                            "              Note that each iteration might involve inner iterations.\n"
                            "Default     : 1e+04");

} // end addOptions

// Set options
void Quantities::setOptions(Options* options)
{

  // Read bool options
  options->valueAsBool("approximate_hessian_initial_scaling", approximate_hessian_initial_scaling_);
  options->valueAsBool("evaluate_function_with_gradient", evaluate_function_with_gradient_);

  // Read double options
  options->valueAsDouble("cpu_time_limit", cpu_time_limit_);
  options->valueAsDouble("inexact_termination_factor_initial", inexact_termination_factor_initial_);
  options->valueAsDouble("inexact_termination_update_factor", inexact_termination_update_factor_);
  options->valueAsDouble("inexact_termination_update_stepsize_threshold", inexact_termination_update_stepsize_threshold_);
  options->valueAsDouble("iterate_norm_tolerance", iterate_norm_tolerance_);
  options->valueAsDouble("scaling_threshold", scaling_threshold_);
  options->valueAsDouble("stationarity_radius_initialization_factor", stationarity_radius_initialization_factor_);
  options->valueAsDouble("stationarity_radius_initialization_minimum", stationarity_radius_initialization_minimum_);
  options->valueAsDouble("stationarity_radius_update_factor", stationarity_radius_update_factor_);
  options->valueAsDouble("stationarity_tolerance", stationarity_tolerance_);
  options->valueAsDouble("trust_region_radius_initialization_factor", trust_region_radius_initialization_factor_);
  options->valueAsDouble("trust_region_radius_initialization_minimum", trust_region_radius_initialization_minimum_);
  options->valueAsDouble("trust_region_radius_update_factor", trust_region_radius_update_factor_);

  // Read integer options
  options->valueAsInteger("function_evaluation_limit", function_evaluation_limit_);
  options->valueAsInteger("gradient_evaluation_limit", gradient_evaluation_limit_);
  options->valueAsInteger("iteration_limit", iteration_limit_);

} // end setOptions

// Initialization
void Quantities::initialize(const std::shared_ptr<Problem> problem)
{

  // Start times
  start_time_ = clock();
  end_time_ = start_time_;
  direction_computation_time_ = 0;
  evaluation_time_ = 0;
  line_search_time_ = 0;

  // Initialize counters
  function_counter_ = 0;
  gradient_counter_ = 0;
  iteration_counter_ = 0;
  inner_iteration_counter_ = 0;
  qp_iteration_counter_ = 0;
  total_inner_iteration_counter_ = 0;
  total_qp_iteration_counter_ = 0;

  // Declare integer
  int n;

  // Get number of variables
  bool evaluation_success = problem->numberOfVariables(n);

  // Check for evaluation success
  if (!evaluation_success) {
    THROW_EXCEPTION(NONOPT_PROBLEM_DATA_FAILURE_EXCEPTION, "Read of number of variables failed.");
  }

  // Set number of variables
  number_of_variables_ = n;

  // Declare vector
  std::shared_ptr<Vector> v(new Vector(number_of_variables_));

  // Get initial point
  evaluation_success = problem->initialPoint(number_of_variables_, v->valuesModifiable());

  // Check for success
  if (!evaluation_success) {
    THROW_EXCEPTION(NONOPT_PROBLEM_DATA_FAILURE_EXCEPTION, "Read of initial point failed.");
  }

  // Declare iterate
  std::shared_ptr<Point> initial_iterate(new Point(problem, v, 1.0));

  // Set initial point
  current_iterate_ = initial_iterate;

  // Initialize direction
  direction_ = std::make_shared<Vector>(number_of_variables_);

  // Initialize direction for termination check
  direction_termination_ = std::make_shared<Vector>(number_of_variables_);

  // Initialize point set
  point_set_ = std::make_shared<std::vector<std::shared_ptr<Point>>>();

  // Initialize stepsize
  stepsize_ = 0.0;

  // Evaluate all functions at current iterate
  evaluateFunctionsAtCurrentIterate();

  // Determine problem scaling
  current_iterate_->determineScale(*this);

  // Scale evaluated objective
  current_iterate_->scaleObjective();

  // Scale evaluated gradient
  current_iterate_->scaleGradient();

  // Initialize stationarity radius
  stationarity_radius_ = fmax(stationarity_radius_initialization_minimum_, stationarity_radius_initialization_factor_ * current_iterate_->gradient()->normInf());

  // Initialize trust region radius
  trust_region_radius_ = fmax(trust_region_radius_initialization_minimum_, trust_region_radius_initialization_factor_ * current_iterate_->gradient()->normInf());

  // Initialize inexact termination factor
  resetInexactTerminationFactor();

  // Store norm of initial point (for termination check)
  iterate_norm_initial_ = current_iterate_->vector()->norm2();

} // end initialize

// Iteration header string
std::string Quantities::iterationHeader()
{
  return "  Iter.  Objective   St. Rad.  Tr. Rad.  |Pts|";
}

// Iteration null values
std::string Quantities::iterationNullValues()
{
  return " ------ ----------- --------- --------- ------";
}

// Evaluate all functions at current iterate
void Quantities::evaluateFunctionsAtCurrentIterate()
{

  // Evaluate objective
  bool evaluation_success;

  // Check whether to evaluate function with gradient
  if (evaluate_function_with_gradient_) {

    // Evaluate function
    evaluation_success = current_iterate_->evaluateObjectiveAndGradient(*this);

    // Check for evaluation success
    if (!evaluation_success) {
      THROW_EXCEPTION(NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION, "Function + gradient evaluation failed.");
    }

  }
  else {

    // Evaluate function
    evaluation_success = current_iterate_->evaluateObjective(*this);

    // Check for evaluation success
    if (!evaluation_success) {
      THROW_EXCEPTION(NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION, "Function evaluation failed.");
    }

    // Evaluate gradient
    evaluation_success = current_iterate_->evaluateGradient(*this);

    // Check for evaluation success
    if (!evaluation_success) {
      THROW_EXCEPTION(NONOPT_GRADIENT_EVALUATION_FAILURE_EXCEPTION, "Gradient evaluation failed.");
    }

  } // end else

} // end evaluateFunctionsAtCurrentIterate

// Reset inexact termination factor
void Quantities::resetInexactTerminationFactor()
{

  // Set termination factor to initial value
  inexact_termination_factor_ = inexact_termination_factor_initial_;

} // end resetInexactTerminationFactor

// Update inexact termination factor
void Quantities::updateInexactTerminationFactor()
{

  // Update inexact termination factor if stepsize is small
  if (stepsize_ < inexact_termination_update_stepsize_threshold_) {
    inexact_termination_factor_ = inexact_termination_factor_ * inexact_termination_update_factor_;
  }

} // end updateInexactTerminationFactor

// Update radii
void Quantities::updateRadii()
{

  // Decrease radii
  stationarity_radius_ = fmax(stationarity_tolerance_, stationarity_radius_update_factor_ * stationarity_radius_);
  trust_region_radius_ = trust_region_radius_update_factor_ * trust_region_radius_;

  // Reset step size
  stepsize_ = 1.0;

} // end updateRadii

// Print header
void Quantities::printHeader(const Reporter* reporter)
{

  // Print header
  reporter->printf(R_NL, R_BASIC, "Number of variables.................. : %d\n", number_of_variables_);

} // end printHeader

// Print iteration values
void Quantities::printIterationValues(const Reporter* reporter)
{

  // Print iteration values
  reporter->printf(R_NL, R_PER_ITERATION, " %6d %+.4e %+.2e %+.2e %6d", iteration_counter_, current_iterate_->objective(), stationarity_radius_, trust_region_radius_, (int)point_set_->size());

} // end printIterationValues

// Print footer
void Quantities::printFooter(const Reporter* reporter)
{

  // Print quantities footer
  reporter->printf(R_NL, R_BASIC, "\n\n"
                                  "Objective............................ : %e\n"
                                  "Objective (unscaled)................. : %e\n"
                                  "\n"
                                  "Number of iterations................. : %d\n"
                                  "Number of inner iterations........... : %d\n"
                                  "Number of QP iterations.............. : %d\n"
                                  "Number of function evaluations....... : %d\n"
                                  "Number of gradient evaluations....... : %d\n"
                                  "\n"
                                  "CPU seconds.......................... : %f\n"
                                  "CPU seconds in evaluations........... : %f\n"
                                  "CPU seconds in direction computations : %f\n"
                                  "CPU seconds in line searches......... : %f\n",
                   current_iterate_->objective(),
                   current_iterate_->objectiveUnscaled(),
                   iteration_counter_,
                   total_inner_iteration_counter_,
                   total_qp_iteration_counter_,
                   function_counter_,
                   gradient_counter_,
                   (end_time_ - start_time_) / (double)CLOCKS_PER_SEC,
                   evaluation_time_ / (double)CLOCKS_PER_SEC,
                   direction_computation_time_ / (double)CLOCKS_PER_SEC,
                   line_search_time_ / (double)CLOCKS_PER_SEC);

} // end printFooter

// Finalization
void Quantities::finalize()
{

  // Set end time
  end_time_ = clock();

} // end finalize

} // namespace NonOpt
