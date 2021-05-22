// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>
#include <ctime>

#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptLineSearchWeakWolfe.hpp"

namespace NonOpt
{

// Add options
void LineSearchWeakWolfe::addOptions(Options* options)
{

  // Add bool options
  options->addBoolOption("LSWW_fail_on_small_interval",
                         false,
                         "Indicator for whether to indicate failure on small interval.\n"
                         "Default     : false.");

  // Add double options
  options->addDoubleOption("LSWW_stepsize_initial",
                           1.0,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Initial stepsize to be used in the first iteration.  Note that\n"
                           "              the initial stepsize used in the line search in subsequent\n"
                           "              iterations is set the minimum of this value and a factor times\n"
                           "              the stepsize accepted in the previous iteration.\n"
                           "Default     : 1.0.");
  options->addDoubleOption("LSWW_stepsize_minimum",
                           1e-20,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining an insufficient stepsize.  If the\n"
                           "              line search yields a stepsize below this tolerance, then the\n"
                           "              algorithm may terminate with a message of a small stepsize.\n"
                           "Default     : 1e-20.");
  options->addDoubleOption("LSWW_stepsize_maximum",
                           1e+02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Maximum stepsize allowed by the line search.\n"
                           "Default     : 1e+02.");
  options->addDoubleOption("LSWW_stepsize_sufficient_decrease_threshold",
                           1e-10,
                           0.0,
                           1.0,
                           "Sufficient decrease constant for the weak Wolfe line search.\n"
                           "Default     : 1e-10.");
  options->addDoubleOption("LSWW_stepsize_sufficient_decrease_fudge_factor",
                           1e-10,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Sufficient decrease fudge factor.\n"
                           "Default     : 1e-10.");
  options->addDoubleOption("LSWW_stepsize_curvature_threshold",
                           9e-01,
                           0.0,
                           1.0,
                           "Curvature condition constant for the weak Wolfe line search.\n"
                           "Default     : 9e-01.");
  options->addDoubleOption("LSWW_stepsize_curvature_fudge_factor",
                           1e-10,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Curvature condition fudge factor.\n"
                           "Default     : 1e-10.");
  options->addDoubleOption("LSWW_stepsize_decrease_factor",
                           5e-01,
                           0.0,
                           1.0,
                           "Factor for updating the stepsize during the line search.\n"
                           "Default     : 5e-01.");
  options->addDoubleOption("LSWW_stepsize_increase_factor",
                           1e+01,
                           1.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for updating the stepsize before the line search.\n"
                           "Default     : 1e+01.");
  options->addDoubleOption("LSWW_stepsize_bound_tolerance",
                           1e-20,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for terminating line search.  If the stepsize is\n"
                           "              greater than the maximum stepsize minus this tolerance,\n"
                           "              then the line search terminates.\n"
                           "Default     : 1e-20.");

} // end addOptions

// Set options
void LineSearchWeakWolfe::setOptions(Options* options)
{

  // Read bool options
  options->valueAsBool("LSWW_fail_on_small_interval", fail_on_small_interval_);

  // Read options
  options->valueAsDouble("LSWW_stepsize_initial", stepsize_initial_);
  options->valueAsDouble("LSWW_stepsize_minimum", stepsize_minimum_);
  options->valueAsDouble("LSWW_stepsize_maximum", stepsize_maximum_);
  options->valueAsDouble("LSWW_stepsize_sufficient_decrease_threshold", stepsize_sufficient_decrease_threshold_);
  options->valueAsDouble("LSWW_stepsize_sufficient_decrease_fudge_factor", stepsize_sufficient_decrease_fudge_factor_);
  options->valueAsDouble("LSWW_stepsize_curvature_threshold", stepsize_curvature_threshold_);
  options->valueAsDouble("LSWW_stepsize_curvature_fudge_factor", stepsize_curvature_fudge_factor_);
  options->valueAsDouble("LSWW_stepsize_decrease_factor", stepsize_decrease_factor_);
  options->valueAsDouble("LSWW_stepsize_increase_factor", stepsize_increase_factor_);
  options->valueAsDouble("LSWW_stepsize_bound_tolerance", stepsize_bound_tolerance_);

} // end setOptions

// Initialize
void LineSearchWeakWolfe::initialize(const Options* options,
                                     Quantities* quantities,
                                     const Reporter* reporter)
{

  // Initialize stepsize
  quantities->setStepsize(fmax(stepsize_minimum_, fmin(stepsize_initial_, stepsize_maximum_)));

} // end initialize

// Run line search
void LineSearchWeakWolfe::runLineSearch(const Options* options,
                                        Quantities* quantities,
                                        const Reporter* reporter,
                                        Strategies* strategies)
{

  // Initialize values
  setStatus(LS_UNSET);
  quantities->setTrialIterateToCurrentIterate();
  clock_t start_time = clock();

  // try line search, terminate on any exception
  try {

    // Declare bool for evaluations
    bool evaluation_success;

    // Check whether to evaluate function with gradient
    if (quantities->evaluateFunctionWithGradient()) {

      // Evaluate function
      evaluation_success = quantities->currentIterate()->evaluateObjectiveAndGradient(*quantities);

      // Check for evaluation success
      if (!evaluation_success) {
        quantities->setStepsize(0.0);
        THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, "Line search unsuccessful. Evaluation failed.");
      }

    }
    else {

      // Evaluate function
      evaluation_success = quantities->currentIterate()->evaluateObjective(*quantities);

      // Check for evaluation success
      if (!evaluation_success) {
        quantities->setStepsize(0.0);
        THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, "Line search unsuccessful. Evaluation failed.");
      }

      // Evaluate gradient
      evaluation_success = quantities->currentIterate()->evaluateGradient(*quantities);

      // Check for evaluation success
      if (!evaluation_success) {
        quantities->setStepsize(0.0);
        THROW_EXCEPTION(LS_EVALUATION_FAILURE_EXCEPTION, "Line search unsuccessful. Evaluation failed.");
      }

    } // end else

    // Compute directional derivative
    double directional_derivative = quantities->currentIterate()->gradient()->innerProduct(*quantities->direction());

    // Initialize stepsize bounds for search
    double stepsize_minimum = stepsize_minimum_;
    double stepsize_maximum = stepsize_maximum_;

    // Initialize stepsize
    quantities->setStepsize(fmax(stepsize_minimum_, fmin(stepsize_increase_factor_ * quantities->stepsize(), fmin(stepsize_initial_, stepsize_maximum_))));

    // Loop
    while (true) {

      // Initialize booleans
      bool sufficient_decrease = false;
      bool curvature_condition = false;

      // Declare new point
      quantities->setTrialIterate(quantities->currentIterate()->makeNewLinearCombination(1.0, quantities->stepsize(), *quantities->direction()));

      // Evaluate trial objective
      if (quantities->evaluateFunctionWithGradient()) {
        evaluation_success = quantities->trialIterate()->evaluateObjectiveAndGradient(*quantities);
      }
      else {
        evaluation_success = quantities->trialIterate()->evaluateObjective(*quantities);
      }

      // Check for evaluation success
      if (evaluation_success) {

        // Check for sufficient decrease
        sufficient_decrease = (quantities->trialIterate()->objective() - quantities->currentIterate()->objective() <= -stepsize_sufficient_decrease_threshold_ * quantities->stepsize() * fmin(strategies->qpSolver()->dualObjectiveQuadraticValue(), fmax(strategies->qpSolver()->combinationTranslatedNorm2Squared(), strategies->qpSolver()->primalSolutionNorm2Squared())) + stepsize_sufficient_decrease_fudge_factor_);

        // Check Armijo condition
        if (sufficient_decrease) {

          // Evaluate new gradient
          evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

          // Check for evaluation success
          if (evaluation_success) {

            // Check for curvature condition
            curvature_condition = (quantities->trialIterate()->gradient()->innerProduct(*quantities->direction()) >= stepsize_curvature_threshold_ * directional_derivative - stepsize_curvature_fudge_factor_);

            // Check curvature condition
            if (curvature_condition) {
              THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");
            }

          } // end if

        } // end if

      } // end if

      // Check if stepsize near bound
      if (quantities->stepsize() <= stepsize_minimum + stepsize_bound_tolerance_ ||
          quantities->stepsize() >= stepsize_maximum - stepsize_bound_tolerance_) {

        // Check for failure on interval
        if (fail_on_small_interval_) {
          THROW_EXCEPTION(LS_INTERVAL_TOO_SMALL_EXCEPTION, "Line search unsuccessful.  Interval too small.");
          printf("small interval");
        }

        // Evaluate objective at trial iterate
        if (quantities->evaluateFunctionWithGradient()) {
          evaluation_success = quantities->trialIterate()->evaluateObjectiveAndGradient(*quantities);
        }
        else {
          evaluation_success = quantities->trialIterate()->evaluateObjective(*quantities);
        }

        // Check for evaluation success
        if (evaluation_success) {

          // Check for decrease
          if (quantities->trialIterate()->objective() < quantities->currentIterate()->objective()) {

            // Evaluate gradient at trial iterate
            evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

            // Check for successful evaluation
            if (evaluation_success) {
              THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");
            }

          } // end if

        } // end if

        // Set null step
        quantities->setStepsize(0.0);

        // Set new point
        quantities->setTrialIterateToCurrentIterate();

        // Terminate
        THROW_EXCEPTION(LS_SUCCESS_EXCEPTION, "Line search successful.");

      } // end if

      // Update lower or upper bound
      if (sufficient_decrease && evaluation_success) {
        stepsize_minimum = quantities->stepsize();
      }
      else {
        stepsize_maximum = quantities->stepsize();
      }

      // Update stepsize
      quantities->setStepsize((1 - stepsize_decrease_factor_) * stepsize_minimum + stepsize_decrease_factor_ * stepsize_maximum);

    } // end while

  } // end try

  // catch exceptions
  catch (LS_SUCCESS_EXCEPTION& exec) {
    setStatus(LS_SUCCESS);
  } catch (LS_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(LS_EVALUATION_FAILURE);
  } catch (LS_STEPSIZE_TOO_SMALL_EXCEPTION& exec) {
    setStatus(LS_STEPSIZE_TOO_SMALL);
  } // end catch

  // Print iterate information
  reporter->printf(R_NL, R_PER_ITERATION, " %+.2e", quantities->stepsize());

  // Increment line search time
  quantities->incrementLineSearchTime(clock() - start_time);

} // end runLineSearch

} // namespace NonOpt
