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

  // Add options
  addOptions();

} // end constructor

// Destructor
NonOptSolver::~NonOptSolver()
{

  // Delete reports
  reporter_.deleteReports();

} // end destructor

// Add options
void NonOptSolver::addOptions()
{

  // Add integer options
  options_.addIntegerOption("print_level",
                            R_BASIC,
                            R_NONE,
                            R_PER_INNER_ITERATION,
                            "Print level to standard output.\n"
                            "Default     : 1");
  options_.addIntegerOption("print_level_file",
                            R_NONE,
                            R_NONE,
                            R_PER_INNER_ITERATION,
                            "Print level to file specified by print_file_name.\n"
                            "Default     : 0");
  options_.addIntegerOption("qp_print_level",
                            R_NONE,
                            R_NONE,
                            R_PER_INNER_ITERATION,
                            "Print level for QP solver to standard output.\n"
                            "Default     : 0");
  options_.addIntegerOption("qp_print_level_file",
                            R_NONE,
                            R_NONE,
                            R_PER_INNER_ITERATION,
                            "Print level for QP solver to file specified qp_print_file_name.\n"
                            "Default     : 0");

  // Add string options
  options_.addStringOption("print_file_name",
                           "nonopt.out",
                           "File name for printing if print_level_file > 0.\n"
                           "Default     : nonopt.out");
  options_.addStringOption("qp_print_file_name",
                           "nonopt_qp.out",
                           "File name for printing for QP solver if qp_print_level_file > 0.\n"
                           "Default     : nonopt_qp.out");

  // Add options for quantities
  quantities_.addOptions(&options_);

  // Add options for strategies
  strategies_.addOptions(&options_);

} // end addOptions

// Set options
void NonOptSolver::setOptions()
{

  // Read integer options
  options_.valueAsInteger("print_level", print_level_);
  options_.valueAsInteger("print_level_file", print_level_file_);
  options_.valueAsInteger("qp_print_level", qp_print_level_);
  options_.valueAsInteger("qp_print_level_file", qp_print_level_file_);

  // Read string options
  options_.valueAsString("print_file_name", print_file_name_);
  options_.valueAsString("qp_print_file_name", qp_print_file_name_);

  // Set standard output stream
  std::shared_ptr<StreamReport> s_out(new StreamReport("default", R_NL,  static_cast<ReportLevel>(print_level_)));
  s_out->setStream(&std::cout);
  reporter_.addReport(s_out);

  // Set file output stream
  if (print_level_file_ > R_NONE) {
    std::shared_ptr<FileReport> f_out(new FileReport("default_file", R_NL, static_cast<ReportLevel>(print_level_file_)));
    f_out->open(print_file_name_.c_str());
    reporter_.addReport(f_out);
  } // end if

  // Set standard output stream for QP solver
  std::shared_ptr<StreamReport> s_qp_out(new StreamReport("default_qp", R_QP, static_cast<ReportLevel>(qp_print_level_)));
  s_qp_out->setStream(&std::cout);
  reporter_.addReport(s_qp_out);

  // Set file output stream for QP solver
  if (qp_print_level_file_ > R_NONE) {
    std::shared_ptr<FileReport> f_qp_out(new FileReport("default_qp_file", R_QP, static_cast<ReportLevel>(qp_print_level_file_)));
    f_qp_out->open(qp_print_file_name_.c_str());
    reporter_.addReport(f_qp_out);
  } // end if

  // Set quantities options
  quantities_.setOptions(&options_);

  // Set strategies options
  strategies_.setOptions(&options_);

  // Print message
  reporter_.printf(R_NL, R_BASIC, options_.message().c_str());

  // Clear message
  options_.resetMessage();

} // end setOptions

// Initialize
void NonOptSolver::initialize(const std::shared_ptr<Problem> problem)
{

  // Initialize quantities
  quantities_.initialize(problem);

  // Initialize strategies
  strategies_.initialize(&options_, &quantities_, &reporter_);

} // end initialize

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

    // Initialize
    initialize(problem);

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
      if (quantities_.currentIterate()->vector()->norm2() >= quantities_.iterateNormTolerance() * fmax(1.0, quantities_.iterateNormInitial())) {
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
        quantities_.resetInexactTerminationFactor();
      }

      // Run line search
      strategies_.lineSearch()->runLineSearch(&options_, &quantities_, &reporter_, &strategies_);

      // Check status
      if (strategies_.lineSearch()->status() != LS_SUCCESS) {
        THROW_EXCEPTION(NONOPT_LINE_SEARCH_FAILURE_EXCEPTION, "Line search failed.");
      }

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

      // Update inexact termination factor (depends on stepsize from line search)
      quantities_.updateInexactTerminationFactor();

      // Update iterate
      quantities_.setCurrentIterate(quantities_.trialIterate());

      // Evaluate all functions at current iterate
      quantities_.evaluateFunctionsAtCurrentIterate();

      // Increment iteration counter
      quantities_.incrementIterationCounter();

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
  } catch (NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE);
  } catch (NONOPT_DERIVATIVE_CHECKER_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_DERIVATIVE_CHECKER_FAILURE);
  } catch (NONOPT_DIRECTION_COMPUTATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_DIRECTION_COMPUTATION_FAILURE);
  } catch (NONOPT_FUNCTION_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_FAILURE);
  } catch (NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE);
  } catch (NONOPT_GRADIENT_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_FAILURE);
  } catch (NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE);
  } catch (NONOPT_LINE_SEARCH_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_LINE_SEARCH_FAILURE);
  } catch (NONOPT_POINT_SET_UPDATE_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_POINT_SET_UPDATE_FAILURE);
  } catch (NONOPT_PROBLEM_DATA_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_PROBLEM_DATA_FAILURE);
  } catch (NONOPT_SYMMETRIC_MATRIX_ASSERT_EXCEPTION& exec) {
    setStatus(NONOPT_SYMMETRIC_MATRIX_ASSERT_FAILURE);
  } catch (NONOPT_TERMINATION_FAILURE_EXCEPTION& exec) {
    setStatus(NONOPT_TERMINATION_FAILURE);
  } catch (NONOPT_VECTOR_ASSERT_EXCEPTION& exec) {
    setStatus(NONOPT_VECTOR_ASSERT_FAILURE);
  }

  // Print end of line
  reporter_.printf(R_NL, R_PER_ITERATION, "\n");

  // Check whether to finalize problem solution
  if (status() != NONOPT_FUNCTION_EVALUATION_FAILURE &&
      status() != NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE &&
      status() != NONOPT_GRADIENT_EVALUATION_FAILURE &&
      status() != NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE &&
      status() != NONOPT_PROBLEM_DATA_FAILURE &&
      status() != NONOPT_VECTOR_ASSERT_FAILURE) {

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
  case NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Approximate Hessian update failure.");
    break;
  case NONOPT_DERIVATIVE_CHECKER_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Derivative checker failure.");
    break;
  case NONOPT_DIRECTION_COMPUTATION_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Direction computation failure.");
    break;
  case NONOPT_FUNCTION_EVALUATION_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Function evaluation failure! Check definition of problem.");
    break;
  case NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Function evaluation assert failure! This wasn't supposed to happen!");
    break;
  case NONOPT_GRADIENT_EVALUATION_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Gradient evaluation failure! Check definition of problem.");
    break;
  case NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Gradient evaluation assert failure! This wasn't supposed to happen!");
    break;
  case NONOPT_LINE_SEARCH_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Line search failure.");
    break;
  case NONOPT_POINT_SET_UPDATE_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Point set update failure.");
    break;
  case NONOPT_PROBLEM_DATA_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Problem data read failure!  Check definition of problem.");
    break;
  case NONOPT_SYMMETRIC_MATRIX_ASSERT_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Symmetric matrix assert failure!  This wasn't supposed to happen!");
    break;
  case NONOPT_TERMINATION_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Termination check failure.");
    break;
  case NONOPT_VECTOR_ASSERT_FAILURE:
    reporter_.printf(R_NL, R_BASIC, "Vector assert failure!  This wasn't supposed to happen!");
    break;
  default:
    reporter_.printf(R_NL, R_BASIC, "Unknown exit status! This wasn't supposed to happen!");
    break;
  } // end switch

  // Check whether to print final data
  if (status() != NONOPT_FUNCTION_EVALUATION_FAILURE &&
      status() != NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE &&
      status() != NONOPT_GRADIENT_EVALUATION_FAILURE &&
      status() != NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE &&
      status() != NONOPT_PROBLEM_DATA_FAILURE &&
      status() != NONOPT_VECTOR_ASSERT_FAILURE) {

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
