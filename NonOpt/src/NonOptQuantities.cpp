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
    : evaluation_time_(0),
      function_counter_(0),
      gradient_counter_(0),
      iteration_counter_(0),
      inner_iteration_counter_(0),
      total_inner_iteration_counter_(0),
      total_qp_iteration_counter_(0),
      number_of_variables_(0),
      stationarity_radius_(0.0),
      trust_region_radius_(0.0),
      stepsize_(0.0)
{
  start_time_ = clock();
  end_time_ = start_time_;
  current_iterate_.reset();
  trial_iterate_.reset();
  direction_.reset();
  point_set_.reset();
}

// Destructor
Quantities::~Quantities() {}

// Add options
void Quantities::addOptions(Options* options,
                            const Reporter* reporter)
{

  // Add double options
  options->addDoubleOption(reporter,
                           "scaling_threshold",
                           1e+02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Threshold for determining objective scaling.  If norm of gradient\n"
                           "at the initial point is greater than this value, then the objective\n"
                           "is scaled so that the initial gradient norm is at this value.\n"
                           "Default value: 1e+02.");
  options->addDoubleOption(reporter,
                           "stationarity_radius_initialization_minimum",
                           1e-02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Minimum initial value for stationarity radius.\n"
                           "Default value: 1e-02.");
  options->addDoubleOption(reporter,
                           "stationarity_radius_initialization_factor",
                           1e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for initializing the stationarity radius.  The initial\n"
                           "stationarity radius is the maximum of this value times the\n"
                           "inf-norm of the gradient at the initial point and\n"
                           "stationarity_radius_initialization_minimum.\n"
                           "Default value: 1e-01.");
  options->addDoubleOption(reporter,
                           "stationarity_radius_update_factor",
                           1e-01,
                           0.0,
                           1.0,
                           "Factor for updating the stationarity radius.  If the conditions\n"
                           "for updating the stationarity and trust region radii are met,\n"
                           "then the stationarity radius is multiplied by this fraction.\n"
                           "Default value: 1e-01.");
  options->addDoubleOption(reporter,
                           "trust_region_radius_initialization_minimum",
                           1e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Minimum initial value for trust region radius.\n"
                           "Default value: 1e-01.");
  options->addDoubleOption(reporter,
                           "trust_region_radius_initialization_factor",
                           1e+04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for initializing the trust region radius.  The initial\n"
                           "trust region radius is the maximum of this value times the\n"
                           "inf-norm of the gradient at the initial point and\n"
                           "trust_region_radius_initialization_minimum.\n"
                           "Default value: 1e+04.");
  options->addDoubleOption(reporter,
                           "trust_region_radius_update_factor",
                           1e-01,
                           0.0,
                           1.0,
                           "Factor for updating the trust region radius.  If the conditions\n"
                           "for updating the stationarity and trust region radii are met,\n"
                           "then the trust region radius is multiplied by this fraction.\n"
                           "Default value: 1e-01.");
  options->addDoubleOption(reporter,
                           "inexact_termination_factor_initial",
                           1.5,
                           0.0,
                           4.0,
                           "Initial inexact termination, if allowed.  Factor by which\n"
                           "norm of inexact solution needs to be within true norm of true\n"
                           "(unknown) projection of origin onto convex hull of gradients.\n"
                           "Default value: 1.5.");
  options->addDoubleOption(reporter,
                           "inexact_termination_update_factor",
                           1.0,
                           0.0,
                           1.0,
                           "Factor for updating the inexact termination factor.  If the conditions\n"
                           "for updating the inexact termination factor is met,\n"
                           "then the inexact termination factor is multiplied by this fraction.\n"
                           "Default value: 0.9999.");

  // Add integer options
  options->addIntegerOption(reporter,
                            "function_evaluation_limit",
                            1e+05,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of function evaluations that will be.\n"
                            "performed.\n"
                            "Default value: 1e+04.");
  options->addIntegerOption(reporter,
                            "gradient_evaluation_limit",
                            5e+04,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of gradient evaluations that will be.\n"
                            "performed.\n"
                            "Default value: 1e+04.");

}  // end addOptions

// Set options
void Quantities::setOptions(const Options* options,
                            const Reporter* reporter)
{

  // Read double options
  options->valueAsDouble(reporter, "scaling_threshold", scaling_threshold_);
  options->valueAsDouble(reporter, "stationarity_radius_initialization_minimum", stationarity_radius_initialization_minimum_);
  options->valueAsDouble(reporter, "stationarity_radius_initialization_factor", stationarity_radius_initialization_factor_);
  options->valueAsDouble(reporter, "stationarity_radius_update_factor", stationarity_radius_update_factor_);
  options->valueAsDouble(reporter, "trust_region_radius_initialization_minimum", trust_region_radius_initialization_minimum_);
  options->valueAsDouble(reporter, "trust_region_radius_initialization_factor", trust_region_radius_initialization_factor_);
  options->valueAsDouble(reporter, "trust_region_radius_update_factor", trust_region_radius_update_factor_);
  options->valueAsDouble(reporter, "inexact_termination_factor_initial", inexact_termination_factor_initial_);
  options->valueAsDouble(reporter, "inexact_termination_update_factor", inexact_termination_update_factor_);
  // Read integer options
  options->valueAsInteger(reporter, "function_evaluation_limit", function_evaluation_limit_);
  options->valueAsInteger(reporter, "gradient_evaluation_limit", gradient_evaluation_limit_);

}  // end setOptions

// Initialization
bool Quantities::initialize(const std::shared_ptr<Problem> problem)
{

  // Start clock
  start_time_ = clock();
  end_time_ = start_time_;

  // Initialize counters
  evaluation_time_ = 0;
  function_counter_ = 0;
  gradient_counter_ = 0;
  iteration_counter_ = 0;
  inner_iteration_counter_ = 0;
  total_inner_iteration_counter_ = 0;
  total_qp_iteration_counter_ = 0;

  // Declare success boolean
  bool success = true;

  // Declare integer
  int n;

  // Get number of variables
  success = problem->numberOfVariables(n);

  // Check for success
  if (!success) {
    return false;
  }

  // Set number of variables
  number_of_variables_ = n;

  // Declare vector
  std::shared_ptr<Vector> v(new Vector(number_of_variables_));

  // Get initial point
  success = problem->initialPoint(number_of_variables_, v->valuesModifiable());

  // Check for success
  if (!success) {
    return false;
  }

  // Declare iterate
  std::shared_ptr<Point> initial_iterate(new Point(problem, v, 1.0));

  // Set initial point
  current_iterate_ = initial_iterate;

  // Initialize direction
  direction_ = std::make_shared<Vector>(number_of_variables_);

  // Initialize point set
  point_set_ = std::make_shared<std::vector<std::shared_ptr<Point> > >();

  // Initialize radii
  stationarity_radius_ = 0.0;
  trust_region_radius_ = 0.0;

  // Initialize inexact termination factor
  inexact_termination_factor_= 0.0;

  // Return
  return success;

}  // end initialize

// Iteration header string
std::string Quantities::iterationHeader()
{
  return "  Iter.   Objective    Stat. Rad.   Trust Rad.  |Points|";
}

// Iteration null values
std::string Quantities::iterationNullValues()
{
  return " ------  -----------  -----------  -----------  --------";
}

// Initialization of radii
void Quantities::initializeRadii(const Options* options,
                                 const Reporter* reporter)
{

  // Initialize stationarity radius
  stationarity_radius_ = fmax(stationarity_radius_initialization_minimum_, stationarity_radius_initialization_factor_ * current_iterate_->gradient()->normInf());

  // Initialize trust region radius
  trust_region_radius_ = fmax(trust_region_radius_initialization_minimum_, trust_region_radius_initialization_factor_ * current_iterate_->gradient()->normInf());

}  // end initializeRadii

// Update radii
void Quantities::updateRadii(double stationarity_tolerance)
{

  // Decrease radii
  stationarity_radius_ = fmax(stationarity_tolerance, stationarity_radius_update_factor_ * stationarity_radius_);
  trust_region_radius_ = trust_region_radius_update_factor_ * trust_region_radius_;

  // reset step size
  stepsize_ = 1.0;
}  // end updateRadii

// Print header
void Quantities::printHeader(const Reporter* reporter)
{

  // Print header
  reporter->printf(R_NL, R_BASIC,
                   "Number of variables.............. : %d\n",
                   number_of_variables_);

}  // end printHeader

// Print iteration values
void Quantities::printIterationValues(const Reporter* reporter)
{

  // Print iteration values
  reporter->printf(R_NL, R_PER_ITERATION,
                   " %6d  %+.4e  %+.4e  %+.4e  %8d",
                   iteration_counter_,
                   current_iterate_->objective(),
                   stationarity_radius_,
                   trust_region_radius_,
                   (int)point_set_->size());

}  // end printIterationValues

// Print footer
void Quantities::printFooter(const Reporter* reporter)
{

  // Print quantities footer
  reporter->printf(R_NL, R_BASIC,
                   "\n\n"
                   "Objective........................ : %e\n"
                   "Objective (unscaled)............. : %e\n"
                   "\n"
                   "Number of iterations............. : %d\n"
                   "Number of inner iterations....... : %d\n"
                   "Number of QP iterations.......... : %d\n"
                   "Number of function evaluations... : %d\n"
                   "Number of gradient evaluations... : %d\n"
                   "\n"
                   "CPU seconds...................... : %f\n"
                   "CPU seconds in NonOpt............ : %f\n"
                   "CPU seconds in evaluations....... : %f\n",
                   current_iterate_->objective(),
                   current_iterate_->objectiveUnscaled(),
                   iteration_counter_,
                   total_inner_iteration_counter_,
                   total_qp_iteration_counter_,
                   function_counter_,
                   gradient_counter_,
                   (end_time_ - start_time_) / (double)CLOCKS_PER_SEC,
                   (end_time_ - start_time_ - evaluation_time_) / (double)CLOCKS_PER_SEC,
                   evaluation_time_ / (double)CLOCKS_PER_SEC);

}  // end printFooter

// Finalization
void Quantities::finalize()
{

  // Set end time
  end_time_ = clock();

}  // end finalize

}  // namespace NonOpt
