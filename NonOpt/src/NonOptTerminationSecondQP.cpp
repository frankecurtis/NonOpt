// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptTerminationSecondQP.hpp"
#include "NonOptDefinitions.hpp"

namespace NonOpt
{

// Add options
void TerminationSecondQP::addOptions(Options* options)
{

  // Add bool options

  // Add double options
  options->addDoubleOption("TS_objective_similarity_tolerance",
                           1e-05,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining objective similarity.  If the\n"
                           "              difference between two consecutive objective values is less\n"
                           "              than this tolerance times max{1,objective value}, then a\n"
                           "              counter is increased.  If the counter exceeds\n"
                           "              objective_similarity_limit, then the stationarity radius\n"
                           "              is decreased or the algorithm terminates.\n"
                           "Default     : 1e-05");
  options->addDoubleOption("TS_objective_tolerance",
                           -NONOPT_DOUBLE_INFINITY,
                           -NONOPT_DOUBLE_INFINITY,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for objective function value.  Algorithm terminates\n"
                           "              if unscaled objective falls below this tolerance.\n"
                           "Default     : -Infinity");
  options->addDoubleOption("TS_stationarity_tolerance_factor",
                           1e+00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for checking termination with respect to stationarity.\n"
                           "              For further explanation, see description for the parameter\n"
                           "              stationarity_tolerance.\n"
                           "Default     : 1e+00");

  // Add integer options
  options->addIntegerOption("TS_objective_similarity_limit",
                            10,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on an objective similarity counter.  If the difference\n"
                            "              between two consecutive objective values is less than\n"
                            "              objective_similarity_tolerance times max{1,objective value},\n"
                            "              then a counter is increased.  If the counter exceeds\n"
                            "              objective_similarity_limit, then the stationarity\n"
                            "              radius is decreased or the algorithm terminates.\n"
                            "Default     : 10");
  options->addIntegerOption("TS_solve_QP_every",
                            10,
                            1,
                            NONOPT_INT_INFINITY,
                            "Dictates that QP should be solved every how many iterations.\n"
                            "Default     : 10");

} // end addOptions

// Set options
void TerminationSecondQP::setOptions(Options* options)
{

  // Read bool options

  // Read double options
  options->valueAsDouble("TS_objective_similarity_tolerance", objective_similarity_tolerance_);
  options->valueAsDouble("TS_objective_tolerance", objective_tolerance_);
  options->valueAsDouble("TS_stationarity_tolerance_factor", stationarity_tolerance_factor_);

  // Read integer options
  options->valueAsInteger("TS_objective_similarity_limit", objective_similarity_limit_);
  options->valueAsInteger("TS_solve_QP_every", solve_QP_every_);

} // end setOptions

// Initialize
void TerminationSecondQP::initialize(const Options* options,
                                     Quantities* quantities,
                                     const Reporter* reporter)
{

  // Initialize counters
  objective_similarity_counter_ = 0;
  solve_QP_counter_ = 0;

  // Set reference values
  objective_reference_ = quantities->currentIterate()->objective();
  stationarity_reference_ = quantities->currentIterate()->gradient()->normInf();

} // end initialize

// Check conditions
void TerminationSecondQP::checkConditions(const Options* options,
                                          Quantities* quantities,
                                          const Reporter* reporter,
                                          Strategies* strategies)
{

  // Initialize status
  setStatus(TE_SUCCESS);

  // Initialize local indicator
  bool direction_conditions = false;

  // Initialize indicators
  terminate_objective_ = false;
  terminate_objective_similarity_ = false;
  terminate_stationary_ = false;
  update_radii_ = false;

  //////////////
  // SOLVE QP //
  //////////////

  // Increment counter
  solve_QP_counter_++;

  // Initialize solve boolean
  bool qp_solved = false;

  // Solve QP
  if (solve_QP_counter_ % solve_QP_every_ == 0) {
    qp_solved = true;
    solveQP(options, quantities, reporter, strategies);
  }

  /////////////////////////////////
  // TERMINATION BY STATIONARITY //
  /////////////////////////////////

  // Evaluate stationarity reference value
  double reference = fmax(1.0, fmax(stationarity_reference_, quantities->currentIterate()->gradient()->normInf()));

  // Check for termination based on stationarity
  if (quantities->stationarityRadius() <= quantities->stationarityTolerance() &&
      strategies->qpSolver(quantities->qpIsSmall())->combinationTranslatedNormInf() <= quantities->stationarityTolerance() * stationarity_tolerance_factor_ * reference) {
    terminate_stationary_ = true;
  }
  else if (quantities->stationarityRadius() <= quantities->stationarityTolerance() &&
           qp_solved &&
           strategies->qpSolverTermination(quantities->qpIsSmall())->combinationTranslatedNormInf() <= quantities->stationarityTolerance() * stationarity_tolerance_factor_ * reference) {
    terminate_stationary_ = true;
  }

  /////////////////////////////////////////
  // TERMINATION BY OBJECTIVE SIMILARITY //
  /////////////////////////////////////////

  // Check objective similarity
  if (objective_reference_ - quantities->currentIterate()->objective() <= objective_similarity_tolerance_ * fmax(1.0, fabs(objective_reference_))) {
    objective_similarity_counter_++;
  }
  else {
    objective_similarity_counter_--;
    objective_similarity_counter_ = fmax(0, objective_similarity_counter_);
  } // end else

  // Update objective reference value
  objective_reference_ = quantities->currentIterate()->objective();

  // Check for termination based on objective changes
  if (quantities->stationarityRadius() <= quantities->stationarityTolerance() &&
      objective_similarity_counter_ > objective_similarity_limit_) {
    terminate_objective_similarity_ = true;
  } // end if

  ////////////////////////////////////////
  // TERMINATION BY OBJECTIVE TOLERANCE //
  ////////////////////////////////////////

  // Check for termination based on objective tolerance
  if (quantities->currentIterate()->objectiveUnscaled() <= objective_tolerance_) {
    terminate_objective_ = true;
  }

  //////////////////
  // RADII UPDATE //
  //////////////////

  // Check radii update conditions
  if (strategies->qpSolver(quantities->qpIsSmall())->primalSolutionNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference &&
      strategies->qpSolver(quantities->qpIsSmall())->combinationNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference &&
      strategies->qpSolver(quantities->qpIsSmall())->combinationTranslatedNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference) {
    direction_conditions = true;
  } // end if
  else if (qp_solved &&
           strategies->qpSolverTermination(quantities->qpIsSmall())->combinationTranslatedNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference) {
    direction_conditions = true;
  }

  // Check for radius update conditions
  if (quantities->stationarityRadius() > quantities->stationarityTolerance() &&
      (direction_conditions || objective_similarity_counter_ > objective_similarity_limit_)) {
    objective_similarity_counter_ = 0;
    update_radii_ = true;
  }

  // Print check values
  if (qp_solved) {
    reporter->printf(R_NL, R_PER_ITERATION, " %+.2e %8d %8d %+.2e %+.2e", quantities->currentIterate()->gradient()->normInf(), strategies->qpSolverTermination(quantities->qpIsSmall())->vectorListLength(), strategies->qpSolverTermination(quantities->qpIsSmall())->numberOfIterations(), strategies->qpSolver(quantities->qpIsSmall())->combinationTranslatedNormInf(), strategies->qpSolverTermination(quantities->qpIsSmall())->combinationTranslatedNormInf());
  }
  else {
    reporter->printf(R_NL, R_PER_ITERATION, " %+.2e %8s %8s %+.2e %9s", quantities->currentIterate()->gradient()->normInf(), "--------", "--------", strategies->qpSolver(quantities->qpIsSmall())->combinationTranslatedNormInf(), "---------");
  }

} // end checkConditions

// Check conditions for use in direction computation
void TerminationSecondQP::checkConditionsDirectionComputation(const Options* options,
                                                              Quantities* quantities,
                                                              const Reporter* reporter,
                                                              Strategies* strategies)
{

  // Initialize indicators
  update_radii_direction_computation_ = false;

  //////////////////
  // RADII UPDATE //
  //////////////////

  // Evaluate stationarity reference value
  double reference = fmax(1.0, fmax(stationarity_reference_, quantities->currentIterate()->gradient()->normInf()));

  // Check radii update conditions
  if (strategies->qpSolver(quantities->qpIsSmall())->primalSolutionNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference &&
      strategies->qpSolver(quantities->qpIsSmall())->combinationNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference &&
      strategies->qpSolver(quantities->qpIsSmall())->combinationTranslatedNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference) {
    update_radii_direction_computation_ = true;
  } // end if

} // end checkConditionsDirectionComputation

// Solve QP
void TerminationSecondQP::solveQP(const Options* options,
                                  Quantities* quantities,
                                  const Reporter* reporter,
                                  Strategies* strategies)
{

  // Initialize values
  strategies->qpSolverTermination(quantities->qpIsSmall())->setPrimalSolutionToZero();
  quantities->resetQPIterationCounter();

  // try QP solve, terminate on any exception
  try {

    // Declare bool for evaluations
    bool evaluation_success;

    // Check whether to evaluate function with gradient
    if (quantities->evaluateFunctionWithGradient()) {

      // Evaluate current objective
      evaluation_success = quantities->currentIterate()->evaluateObjectiveAndGradient(*quantities);

      // Check for successful evaluation
      if (!evaluation_success) {
        THROW_EXCEPTION(TE_EVALUATION_FAILURE_EXCEPTION, "Termination check unsuccessful. Evaluation failed.");
      }
    }
    else {

      // Evaluate current objective
      evaluation_success = quantities->currentIterate()->evaluateObjective(*quantities);

      // Check for successful evaluation
      if (!evaluation_success) {
        THROW_EXCEPTION(TE_EVALUATION_FAILURE_EXCEPTION, "Termination check unsuccessful. Evaluation failed.");
      }

      // Evaluate current gradient
      evaluation_success = quantities->currentIterate()->evaluateGradient(*quantities);

      // Check for successful evaluation
      if (!evaluation_success) {
        THROW_EXCEPTION(TE_EVALUATION_FAILURE_EXCEPTION, "Termination check unsuccessful. Evaluation failed.");
      }

    } // end else

    // Set QP scalar values
    strategies->qpSolverTermination(quantities->qpIsSmall())->setScalar(NONOPT_DOUBLE_INFINITY);
    strategies->qpSolverTermination(quantities->qpIsSmall())->setInexactSolutionTolerance(quantities->stationarityRadius());

    // Declare QP quantities
    std::vector<std::shared_ptr<Vector>> QP_gradient_list;
    std::vector<double> QP_vector;

    // Add pointer to current gradient to list
    QP_gradient_list.push_back(quantities->currentIterate()->gradient());

    // Add linear term value
    QP_vector.push_back(quantities->currentIterate()->objective());

    // Loop through point set
    for (int point_count = 0; point_count < (int)quantities->pointSet()->size(); point_count++) {

      // Create difference vector
      std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *(*quantities->pointSet())[point_count]->vector());

      // Check distance between point and current iterate
      if (difference->normInf() <= quantities->stationarityRadius()) {

        // Evaluate gradient
        evaluation_success = (*quantities->pointSet())[point_count]->evaluateGradient(*quantities);

        // Check for successful evaluation
        if (evaluation_success) {

          // Add pointer to gradient in point set to list
          QP_gradient_list.push_back((*quantities->pointSet())[point_count]->gradient());

          // Add linear term value
          QP_vector.push_back(quantities->currentIterate()->objective());

        } // end if

      } // end if

    } // end for

    // Set QP vector list and linear term
    strategies->qpSolverTermination(quantities->qpIsSmall())->setVectorList(QP_gradient_list);
    strategies->qpSolverTermination(quantities->qpIsSmall())->setVector(QP_vector);

    // Solve QP
    strategies->qpSolverTermination(quantities->qpIsSmall())->solveQP(options, reporter, quantities);

    // Increment QP iteration counter
    quantities->incrementQPIterationCounter(strategies->qpSolverTermination(quantities->qpIsSmall())->numberOfIterations());

    // Get primal solution
    strategies->qpSolverTermination(quantities->qpIsSmall())->primalSolution(quantities->directionTermination()->valuesModifiable());

    // Increment total QP iteration counter
    quantities->incrementTotalQPIterationCounter();

  } // end try

  // catch exceptions
  catch (TE_SUCCESS_EXCEPTION& exec) {
    setStatus(TE_SUCCESS);
  } catch (TE_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(TE_EVALUATION_FAILURE);
  }

} // end solveQP

} // namespace NonOpt
