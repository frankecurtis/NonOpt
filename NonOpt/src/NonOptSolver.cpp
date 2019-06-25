// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <iostream>

#include "NonOptBLAS.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptException.hpp"
#include "NonOptSolver.hpp"
#include "NonOptVersion.hpp"

namespace NonOpt
{

// Constructor
NonOptSolver::NonOptSolver()
{

  // Declare stream report
  std::shared_ptr<StreamReport> s(new StreamReport("s", R_NL, R_PER_ITERATION));

  // Set stream report to standard output
  s->setStream(&std::cout);

  // Add stream report to reporter
  reporter_.addReport(s);

  // Add options
  addOptions();

}  // end constructor

// Destructor
NonOptSolver::~NonOptSolver()
{

  // Delete reporter
  reporter_.deleteReports();

}  // end destructor

// Add options
void NonOptSolver::addOptions()
{

  // Add double options (name, default value, lower bound, upper bound, description)
  options_.addDoubleOption(&reporter_,
                           "cpu_time_limit",
                           1e+04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Limit on the number of CPU seconds.  This limit is only checked\n"
                           "at the beginning of an iteration, so the true CPU time limit\n"
                           "also depends on the time required to a complete an iteration.\n"
                           "Default value: 1e+04.");
  options_.addDoubleOption(&reporter_,
                           "iterate_norm_tolerance",
                           1e+20,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining divergence of the algorithm iterates.\n"
                           "If the norm of an iterate is larger than this tolerance times\n"
                           "the maximum of 1.0 and the norm of the initial iterate, then\n"
                           "the algorithm terminates with a message of divergence.\n"
                           "Default value: 1e+20.");
  options_.addDoubleOption(&reporter_,
                           "stationarity_tolerance",
                           1e-04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining stationarity.  If the stationarity\n"
                           "radius falls below this tolerance and a computed convex\n"
                           "combination of gradients has norm below this tolerance times\n"
                           "the parameter stationarity_tolerance_factor, then the algorithm\n"
                           "terminates with a message of stationarity.\n"
                           "Default value: 1e-04.");
  options_.addDoubleOption(&reporter_,
                           "stationarity_tolerance_factor",
                           1e+00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for checking termination with respect to stationarity.\n"
                           "For further explanation, see description for the parameter\n"
                           "stationarity_tolerance.\n"
                           "Default value: 1e+00.");

  // Add integer options (name, default value, lower bound, upper bound, description)
  options_.addIntegerOption(&reporter_,
                            "iteration_limit",
                            1e+04,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of iterations that will be performed.\n"
                            "Note that each iteration might involve inner iterations.\n"
                            "Default value: 1e+04.");

  // Add options for quantities
  quantities_.addOptions(&options_, &reporter_);

  // Add options for strategies
  strategies_.addOptions(&options_, &reporter_);

}  // end addOptions

// Set options
void NonOptSolver::setOptions()
{

  // Set double options
  options_.valueAsDouble(&reporter_, "cpu_time_limit", cpu_time_limit_);
  options_.valueAsDouble(&reporter_, "iterate_norm_tolerance", iterate_norm_tolerance_);
  options_.valueAsDouble(&reporter_, "stationarity_tolerance", stationarity_tolerance_);
  options_.valueAsDouble(&reporter_, "stationarity_tolerance_factor", stationarity_tolerance_factor_);
  options_.valueAsDouble(&reporter_, "stationarity_radius_update_factor", stationarity_radius_update_factor_);
  options_.valueAsDouble(&reporter_, "trust_region_radius_update_factor", trust_region_radius_update_factor_);

  // Set integer options
  options_.valueAsInteger(&reporter_, "iteration_limit", iteration_limit_);

  // Set quantities options
  quantities_.setOptions(&options_, &reporter_);

  // Set strategies options
  strategies_.setOptions(&options_, &reporter_);

}  // end setOptions

// Solution
void NonOptSolver::solution(double vector[])
{

  // Set inputs for blas
  int length = (int)quantities_.numberOfVariables();
  int increment = 1;

  // Copy values
  dcopy_(&length, quantities_.currentIterate()->vector()->values(), &increment, vector, &increment);

}  // end solution

// Optimize
void NonOptSolver::optimize(const std::shared_ptr<Problem> problem)
{

  // Initialize solver status
  setStatus(NONOPT_UNSET);

  // (Re)set options
  setOptions();

  // try to run algorithm, terminate on any exception
  try {

    // (Re)initialize quantities
    bool initialization_success = quantities_.initialize(problem);

    // Check for initialization success
    if (!initialization_success) {
      THROW_EXCEPTION(NONOPT_INITIALIZATION_FAILURE_EXCEPTION, "Initialization failed.");
    }

    // Evaluate all functions at current iterate
    evaluateFunctionsAtCurrentIterate();

    // Determine problem scaling
    quantities_.currentIterate()->determineScale(quantities_);

    // Scale evaluated gradient
    quantities_.currentIterate()->scaleGradient();

    // Initialize radii
    quantities_.initializeRadii(&options_, &reporter_);

    // Store norm of initial point (for termination check)
    double initial_iterate_norm = quantities_.currentIterate()->vector()->norm2();

    // Initialize strategies
    strategies_.initialize(&options_, &quantities_, &reporter_);

    // Print header
    printHeader();

    // Set iteration header
    strategies_.setIterationHeader();

    // (Outer) Loop
    while (true) {

      // Print iteration header
      printIterationHeader();

      // Print quantities iteration values
      quantities_.printIterationValues(&reporter_);

      // Flush buffer
      reporter_.flushBuffer();

      // Check termination conditions
      if (quantities_.iterationCounter() >= iteration_limit_) {
        THROW_EXCEPTION(NONOPT_ITERATION_LIMIT_EXCEPTION, "Iteration limit has been reached.");
      }
      if ((clock() - quantities_.startTime()) / (double)CLOCKS_PER_SEC >= cpu_time_limit_) {
        THROW_EXCEPTION(NONOPT_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
      }
      if (quantities_.currentIterate()->vector()->norm2() >= iterate_norm_tolerance_ * fmax(1.0, initial_iterate_norm)) {
        THROW_EXCEPTION(NONOPT_ITERATE_NORM_LIMIT_EXCEPTION, "Iterates appear to be diverging.");
      }

      // Compute direction
      strategies_.directionComputation()->computeDirection(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.directionComputation()->status() != DC_SUCCESS) {
        THROW_EXCEPTION(NONOPT_DIRECTION_COMPUTATION_FAILURE_EXCEPTION, "Direction computation failed.");
      }

      // Check termination conditions
      if (quantities_.stationarityRadius() <= stationarity_tolerance_ &&
          strategies_.qpSolver()->dualStepNormInf() <= stationarity_tolerance_ * stationarity_tolerance_factor_ &&
          strategies_.qpSolver()->combinationNormInf() <= stationarity_tolerance_ * stationarity_tolerance_factor_ &&
          strategies_.qpSolver()->combinationTranslatedNormInf() <= stationarity_tolerance_ * stationarity_tolerance_factor_) {
        THROW_EXCEPTION(NONOPT_SUCCESS_EXCEPTION, "Stationary point found.");
      }

      // Check radius update conditions
      if (stationarity_tolerance_ < quantities_.stationarityRadius() &&
          strategies_.qpSolver()->dualStepNormInf() <= quantities_.stationarityRadius() * stationarity_tolerance_factor_ &&
          strategies_.qpSolver()->combinationNormInf() <= quantities_.stationarityRadius() * stationarity_tolerance_factor_ &&
          strategies_.qpSolver()->combinationTranslatedNormInf() <= quantities_.stationarityRadius() * stationarity_tolerance_factor_) {
        quantities_.updateRadii(stationarity_tolerance_);
      }

      // Run line search
      strategies_.lineSearch()->runLineSearch(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.lineSearch()->status() != LS_SUCCESS) {
        THROW_EXCEPTION(NONOPT_LINE_SEARCH_FAILURE_EXCEPTION, "Line search failed.");
      }

      // Update inverse Hessian
      strategies_.inverseHessianUpdate()->updateInverseHessian(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.inverseHessianUpdate()->status() != IH_SUCCESS) {
        THROW_EXCEPTION(NONOPT_INVERSE_HESSIAN_UPDATE_FAILURE_EXCEPTION, "Inverse Hessian update failed.");
      }

      // Add current iterate to point set
      if (quantities_.stepsize() > 0.0) {
        quantities_.pointSet()->push_back(quantities_.currentIterate());
      }

      // Update iterate
      quantities_.setCurrentIterate(quantities_.trialIterate());

      // Increment iteration counter
      quantities_.incrementIterationCounter();

      // Evaluate all functions at current iterate
      evaluateFunctionsAtCurrentIterate();

      // Update point set
      strategies_.pointSetUpdate()->updatePointSet(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.pointSetUpdate()->status() != PS_SUCCESS) {
        THROW_EXCEPTION(NONOPT_POINT_SET_UPDATE_FAILURE_EXCEPTION, "Point set update failed.");
      }

      // Print end of line
      reporter_.printf(R_NL, R_PER_ITERATION, "\n");

    }  // end while

  }  // end try

  // catch exceptions
  catch (NONOPT_SUCCESS_EXCEPTION& exec) {
    setStatus(NONOPT_SUCCESS);
  }
  catch (NONOPT_CPU_TIME_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_CPU_TIME_LIMIT);
  }
  catch (NONOPT_ITERATE_NORM_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_ITERATE_NORM_LIMIT);
  }
  catch (NONOPT_ITERATION_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_ITERATION_LIMIT);
  }
  catch (NONOPT_FUNCTION_EVALUATION_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_LIMIT);
  }
  catch (NONOPT_GRADIENT_EVALUATION_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_LIMIT);
  }
  catch (NONOPT_INITIALIZATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_INITIALIZATION_FAILURE);
  }
  catch (NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_FAILURE);
  }
  catch (NONOPT_GRADIENT_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_FAILURE);
  }
  catch (NONOPT_FUNCTION_EVALUATION_ASSERT_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_ASSERT);
  }
  catch (NONOPT_GRADIENT_EVALUATION_ASSERT_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_ASSERT);
  }
  catch (NONOPT_DIRECTION_COMPUTATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_DIRECTION_COMPUTATION_FAILURE);
  }
  catch (NONOPT_LINE_SEARCH_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_LINE_SEARCH_FAILURE);
  }
  catch (NONOPT_INVERSE_HESSIAN_UPDATE_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_INVERSE_HESSIAN_UPDATE_FAILURE);
  }
  catch (NONOPT_POINT_SET_UPDATE_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_POINT_SET_UPDATE_FAILURE);
  }

  // Print end of line
  reporter_.printf(R_NL, R_PER_ITERATION, "\n");

  // Check whether to finalize problem solution
  if (status() != NONOPT_INITIALIZATION_FAILURE &&
      status() != NONOPT_FUNCTION_EVALUATION_FAILURE &&
      status() != NONOPT_GRADIENT_EVALUATION_FAILURE) {

    // Finalize problem solution
    problem->finalizeSolution(quantities_.numberOfVariables(),
                              quantities_.currentIterate()->vector()->values(),
                              quantities_.currentIterate()->objective(),
                              quantities_.currentIterate()->gradient()->values());
  }

  // Finalize
  quantities_.finalize();

  // Print footer
  printFooter();

}  // end optimize

// Evaluate all functions at current iterate
void NonOptSolver::evaluateFunctionsAtCurrentIterate()
{

  // Evaluate objective
  bool evaluation_success = quantities_.currentIterate()->evaluateObjective(quantities_);

  // Check for evaluation success
  if (!evaluation_success) {
    THROW_EXCEPTION(NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION, "Function evaluation failed.");
  }

  // Evaluate gradient
  evaluation_success = quantities_.currentIterate()->evaluateGradient(quantities_);

  // Check for evaluation success
  if (!evaluation_success) {
    THROW_EXCEPTION(NONOPT_GRADIENT_EVALUATION_FAILURE_EXCEPTION, "Initialization failed.");
  }

}  // end evaluateFunctionsAtCurrentIterate

// Print footer
void NonOptSolver::printFooter()
{

  // Print footer
  reporter_.printf(R_NL, R_BASIC, "\nEXIT: ");

  // Print exit status
  switch (status()) {
    case NONOPT_UNSET:
      reporter_.printf(R_NL, R_BASIC, "Exit status wasn't set! This wasn't supposed to happen!");
      break;
    case NONOPT_SUCCESS:
      reporter_.printf(R_NL, R_BASIC, "Stationary point found.");
      break;
    case NONOPT_CPU_TIME_LIMIT:
      reporter_.printf(R_NL, R_BASIC, "CPU time limit reached.");
      break;
    case NONOPT_ITERATE_NORM_LIMIT:
      reporter_.printf(R_NL, R_BASIC, "Iterates seem to be diverging.");
      break;
    case NONOPT_ITERATION_LIMIT:
      reporter_.printf(R_NL, R_BASIC, "Iteration limit reached.");
      break;
    case NONOPT_FUNCTION_EVALUATION_LIMIT:
      reporter_.printf(R_NL, R_BASIC, "Function evaluation limit reached.");
      break;
    case NONOPT_GRADIENT_EVALUATION_LIMIT:
      reporter_.printf(R_NL, R_BASIC, "Gradient evaluation limit reached.");
      break;
    case NONOPT_INITIALIZATION_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Initialization failure! Check definition of problem.");
      break;
    case NONOPT_FUNCTION_EVALUATION_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Function evaluation failure! Check definition of problem.");
      break;
    case NONOPT_GRADIENT_EVALUATION_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Gradient evaluation failure! Check definition of problem.");
      break;
    case NONOPT_FUNCTION_EVALUATION_ASSERT:
      reporter_.printf(R_NL, R_BASIC, "Function evaluation assert failure! This wasn't supposed to happen!");
      break;
    case NONOPT_GRADIENT_EVALUATION_ASSERT:
      reporter_.printf(R_NL, R_BASIC, "Gradient evaluation assert failure! This wasn't supposed to happen!");
      break;
    case NONOPT_DIRECTION_COMPUTATION_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Direction computation failure.");
      break;
    case NONOPT_LINE_SEARCH_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Line search failure.");
      break;
    case NONOPT_INVERSE_HESSIAN_UPDATE_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Inverse Hessian update failure.");
      break;
    case NONOPT_POINT_SET_UPDATE_FAILURE:
      reporter_.printf(R_NL, R_BASIC, "Point set update failure.");
      break;
    default:
      reporter_.printf(R_NL, R_BASIC, "Unknown exit status! This wasn't supposed to happen!");
      break;
  }  // end switch

  // Check whether to print final data
  if (status() != NONOPT_INITIALIZATION_FAILURE &&
      status() != NONOPT_FUNCTION_EVALUATION_FAILURE) {

    // Print quantities footer
    quantities_.printFooter(&reporter_);

    // Print strategies footer
    strategies_.printFooter(&reporter_);

  }  // end if

}  // end printFooter

// Print header
void NonOptSolver::printHeader()
{

  // Print header
  reporter_.printf(R_NL, R_BASIC,
                   "+--------------------------------------------------------------+\n"
                   "|            NonOpt = Nonsmooth Optimization Solver            |\n"
                   "| NonOpt is released as open source code under the MIT License |\n"
                   "| Please visit http://coral.ise.lehigh.edu/frankecurtis/nonopt |\n"
                   "+--------------------------------------------------------------+\n"
                   "\n"
                   "This is NonOpt version %s\n"
                   "\n",
                   NONOPT_VERSION);

  // Print quantities header
  quantities_.printHeader(&reporter_);

  // Print strategies header
  strategies_.printHeader(&reporter_);

}  // end printHeader

// Print iteration header
void NonOptSolver::printIterationHeader()
{

  if (quantities_.iterationCounter() == 0) {
    reporter_.printf(R_NL, R_PER_ITERATION, "\n");
  }
  if (quantities_.iterationCounter() % 20 == 0) {
    std::string b(quantities_.iterationHeader().length() + strategies_.iterationHeader().length(), '-');
    reporter_.printf(R_NL, R_PER_ITERATION, "%s\n", (b + "\n" + quantities_.iterationHeader() + strategies_.iterationHeader() + "\n" + b).c_str());
  }  // end if

}  // end printIterationHeader

}  // namespace NonOpt