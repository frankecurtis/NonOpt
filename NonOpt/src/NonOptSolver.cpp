// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>
#include <iostream>

#include "NonOptBLASLAPACK.hpp"
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
  std::shared_ptr<StreamReport> s(new StreamReport("default", R_NL, R_BASIC));

  // Set stream report to standard output
  s->setStream(&std::cout);

  // Add stream report to reporter
  reporter_.addReport(s);

  // Add options
  addOptions();

} // end constructor

// Destructor
NonOptSolver::~NonOptSolver()
{

  // Delete reporter
  reporter_.deleteReports();

} // end destructor

// Add options
void NonOptSolver::addOptions()
{

  // Add options for quantities
  quantities_.addOptions(&options_, &reporter_);

  // Add options for strategies
  strategies_.addOptions(&options_, &reporter_);

} // end addOptions

// Set options
void NonOptSolver::setOptions()
{

  // Set quantities options
  quantities_.setOptions(&options_, &reporter_);

  // Set strategies options
  strategies_.setOptions(&options_, &reporter_);

} // end setOptions

// Solution
void NonOptSolver::solution(double vector[])
{

  // Set inputs for BLASLAPACK
  int length = (int)quantities_.numberOfVariables();
  int increment = 1;

  // Copy values
  dcopy_(&length, quantities_.currentIterate()->vector()->values(), &increment, vector, &increment);

} // end solution

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

    // Scale evaluated objective
    quantities_.currentIterate()->scaleObjective();

    // Scale evaluated gradient
    quantities_.currentIterate()->scaleGradient();

    // Initialize radii
    quantities_.initializeRadii(&options_, &reporter_);

    // Initialize inexact termination factor
    quantities_.initializeInexactTerminationFactor(&options_, &reporter_);

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
      if (quantities_.iterationCounter() >= quantities_.iterationLimit()) {
        THROW_EXCEPTION(NONOPT_ITERATION_LIMIT_EXCEPTION, "Iteration limit has been reached.");
      }
      if ((clock() - quantities_.startTime()) / (double)CLOCKS_PER_SEC >= quantities_.cpuTimeLimit()) {
        THROW_EXCEPTION(NONOPT_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
      }
      if (quantities_.currentIterate()->vector()->norm2() >= quantities_.iterateNormTolerance() * fmax(1.0, initial_iterate_norm)) {
        THROW_EXCEPTION(NONOPT_ITERATE_NORM_LIMIT_EXCEPTION, "Iterates appear to be diverging.");
      }

      // Check derivatives
      strategies_.derivativeChecker()->checkDerivatives(&options_, &quantities_, &reporter_);

      // Check status
      if (strategies_.derivativeChecker()->status() != DE_SUCCESS) {
        THROW_EXCEPTION(NONOPT_DERIVATIVE_CHECKER_FAILURE_EXCEPTION, "Derivative checker failed.");
      }

      // Compute direction
      strategies_.directionComputation()->computeDirection(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.directionComputation()->status() != DC_SUCCESS) {
        THROW_EXCEPTION(NONOPT_DIRECTION_COMPUTATION_FAILURE_EXCEPTION, "Direction computation failed.");
      }

      // Check radius update and termination conditions
      strategies_.termination()->checkConditions(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.termination()->status() != TE_SUCCESS) {
        THROW_EXCEPTION(NONOPT_TERMINATION_FAILURE_EXCEPTION, "Termination check failed.");
      }

      // Check final termination conditions
      if (strategies_.termination()->terminateStationary()) {
        THROW_EXCEPTION(NONOPT_SUCCESS_EXCEPTION, "Stationary point found.");
      }
      else if (strategies_.termination()->terminateObjective()) {
        THROW_EXCEPTION(NONOPT_OBJECTIVE_SIMILARITY_EXCEPTION, "Insufficient objective improvement.");
      }
      else if (strategies_.termination()->updateRadii()) {
        quantities_.updateRadii();
        quantities_.initializeInexactTerminationFactor(&options_, &reporter_);
      }

      // Run line search
      strategies_.lineSearch()->runLineSearch(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.lineSearch()->status() != LS_SUCCESS) {
        THROW_EXCEPTION(NONOPT_LINE_SEARCH_FAILURE_EXCEPTION, "Line search failed.");
      }

      // Update inexact termination factor (depends on stepsize from line search)
      quantities_.updateInexactTerminationFactor();

      // Update approximate Hessian
      strategies_.approximateHessianUpdate()->updateApproximateHessian(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.approximateHessianUpdate()->status() != AH_SUCCESS) {
        THROW_EXCEPTION(NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE_EXCEPTION, "Approximate Hessian update failed.");
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

    } // end while

  } // end try

  // catch exceptions
  catch (NONOPT_SUCCESS_EXCEPTION& exec) {
    setStatus(NONOPT_SUCCESS);
  } catch (NONOPT_OBJECTIVE_SIMILARITY_EXCEPTION& exec) {
    setStatus(NONOPT_OBJECTIVE_SIMILARITY);
  } catch (NONOPT_CPU_TIME_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_CPU_TIME_LIMIT);
  } catch (NONOPT_ITERATE_NORM_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_ITERATE_NORM_LIMIT);
  } catch (NONOPT_ITERATION_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_ITERATION_LIMIT);
  } catch (NONOPT_FUNCTION_EVALUATION_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_LIMIT);
  } catch (NONOPT_GRADIENT_EVALUATION_LIMIT_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_LIMIT);
  } catch (NONOPT_INITIALIZATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_INITIALIZATION_FAILURE);
  } catch (NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_FAILURE);
  } catch (NONOPT_GRADIENT_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_FAILURE);
  } catch (NONOPT_FUNCTION_EVALUATION_ASSERT_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_ASSERT);
  } catch (NONOPT_GRADIENT_EVALUATION_ASSERT_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_ASSERT);
  } catch (NONOPT_DERIVATIVE_CHECKER_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_DERIVATIVE_CHECKER_FAILURE);
  } catch (NONOPT_DIRECTION_COMPUTATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_DIRECTION_COMPUTATION_FAILURE);
  } catch (NONOPT_LINE_SEARCH_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_LINE_SEARCH_FAILURE);
  } catch (NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE);
  } catch (NONOPT_POINT_SET_UPDATE_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_POINT_SET_UPDATE_FAILURE);
  } catch (NONOPT_TERMINATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_TERMINATION_FAILURE);
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

} // end optimize

// Evaluate all functions at current iterate
void NonOptSolver::evaluateFunctionsAtCurrentIterate()
{

  // Evaluate objective
  bool evaluation_success;

  // Check whether to evaluate function with gradient
  if (quantities_.evaluateFunctionWithGradient()) {

    // Evaluate function
    evaluation_success = quantities_.currentIterate()->evaluateObjectiveAndGradient(quantities_);

    // Check for evaluation success
    if (!evaluation_success) {
      THROW_EXCEPTION(NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION, "Function + gradient evaluation failed.");
    }

  }
  else {

    // Evaluate function
    evaluation_success = quantities_.currentIterate()->evaluateObjective(quantities_);

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

  } // end else

} // end evaluateFunctionsAtCurrentIterate

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
  case NONOPT_OBJECTIVE_SIMILARITY:
    reporter_.printf(R_NL, R_BASIC, "Objective not improving sufficiently.");
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
  case NONOPT_DERIVATIVE_CHECKER_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Derivative checker failure.");
    break;
  case NONOPT_DIRECTION_COMPUTATION_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Direction computation failure.");
    break;
  case NONOPT_LINE_SEARCH_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Line search failure.");
    break;
  case NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Approximate Hessian update failure.");
    break;
  case NONOPT_POINT_SET_UPDATE_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Point set update failure.");
    break;
  case NONOPT_TERMINATION_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Termination check failure.");
    break;
  default:
    reporter_.printf(R_NL, R_BASIC, "Unknown exit status! This wasn't supposed to happen!");
    break;
  } // end switch

  // Check whether to print final data
  if (status() != NONOPT_INITIALIZATION_FAILURE &&
      status() != NONOPT_FUNCTION_EVALUATION_FAILURE) {

    // Print quantities footer
    quantities_.printFooter(&reporter_);

    // Print strategies footer
    strategies_.printFooter(&reporter_);

  } // end if

} // end printFooter

// Print header
void NonOptSolver::printHeader()
{

  // Print header
  reporter_.printf(R_NL, R_BASIC, "+--------------------------------------------------------------+\n"
                                  "|       NonOpt = Nonlinear/Nonconvex/Nonsmooth Optimizer       |\n"
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

} // end printHeader

// Print iteration header
void NonOptSolver::printIterationHeader()
{

  if (quantities_.iterationCounter() == 0) {
    reporter_.printf(R_NL, R_PER_ITERATION, "\n");
  }
  if (quantities_.iterationCounter() % 20 == 0) {
    std::string b(quantities_.iterationHeader().length() + strategies_.iterationHeader().length(), '-');
    reporter_.printf(R_NL, R_PER_ITERATION, "%s\n", (b + "\n" + quantities_.iterationHeader() + strategies_.iterationHeader() + "\n" + b).c_str());
  } // end if

} // end printIterationHeader

} // namespace NonOpt
