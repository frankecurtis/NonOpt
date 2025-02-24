// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptTerminationBasic.hpp"
#include "NonOptDefinitions.hpp"

namespace NonOpt
{

// Add options
void TerminationBasic::addOptions(Options* options)
{

  // Add bool options

  // Add double options
  options->addDoubleOption("TB_objective_similarity_tolerance",
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
  options->addDoubleOption("TB_objective_tolerance",
                           -NONOPT_DOUBLE_INFINITY,
                           -NONOPT_DOUBLE_INFINITY,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for objective function value.  Algorithm terminates\n"
                           "              if unscaled objective falls below this tolerance.\n"
                           "Default     : -Infinity");
  options->addDoubleOption("TB_stationarity_tolerance_factor",
                           1e+00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for checking termination with respect to stationarity.\n"
                           "              For further explanation, see description for the parameter\n"
                           "              stationarity_tolerance.\n"
                           "Default     : 1e+00");

  // Add integer options
  options->addIntegerOption("TB_objective_similarity_limit",
                            10,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Limit on an objective similarity counter.  If the difference\n"
                            "              between two consecutive objective values is less than\n"
                            "              objective_similarity_tolerance times max{1,objective value},\n"
                            "              then a counter is increased.  If the counter exceeds\n"
                            "              objective_similarity_limit, then the stationarity\n"
                            "              radius is decreased or the algorithm terminates.\n"
                            "Default     : 10");

} // end addOptions

// Set options
void TerminationBasic::setOptions(Options* options)
{

  // Read bool options

  // Read double options
  options->valueAsDouble("TB_objective_similarity_tolerance", objective_similarity_tolerance_);
  options->valueAsDouble("TB_objective_tolerance", objective_tolerance_);
  options->valueAsDouble("TB_stationarity_tolerance_factor", stationarity_tolerance_factor_);

  // Read integer options
  options->valueAsInteger("TB_objective_similarity_limit", objective_similarity_limit_);

} // end setOptions

// Initialize
void TerminationBasic::initialize(const Options* options,
                                  Quantities* quantities,
                                  const Reporter* reporter)
{

  // Initialize objective similarity counter
  objective_similarity_counter_ = 0;

  // Set reference values
  objective_reference_ = quantities->currentIterate()->objective();
  stationarity_reference_ = quantities->currentIterate()->gradient()->normInf();

} // end initialize

// Check conditions
void TerminationBasic::checkConditions(const Options* options,
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

  // Check for radius update conditions
  if (quantities->stationarityRadius() > quantities->stationarityTolerance() &&
      (direction_conditions || objective_similarity_counter_ > objective_similarity_limit_)) {
    objective_similarity_counter_ = 0;
    update_radii_ = true;
  }

  // Print check values
  reporter->printf(R_NL, R_PER_ITERATION, " %+.2e %+.2e", quantities->currentIterate()->gradient()->normInf(), strategies->qpSolver(quantities->qpIsSmall())->combinationTranslatedNormInf());

} // end checkConditions

// Check conditions for use in direction computation
void TerminationBasic::checkConditionsDirectionComputation(const Options* options,
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

} // namespace NonOpt
