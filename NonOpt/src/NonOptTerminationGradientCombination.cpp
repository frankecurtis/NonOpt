// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptTerminationGradientCombination.hpp"
#include "NonOptDefinitions.hpp"

namespace NonOpt
{

// Add options
void TerminationGradientCombination::addOptions(Options* options,
                                                const Reporter* reporter)
{

  // Add bool options

  // Add double options
  options->addDoubleOption(reporter,
                           "stationarity_tolerance",
                           1e-04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining stationarity.  If the stationarity\n"
                           "              radius falls below this tolerance and a computed convex\n"
                           "              combination of gradients has norm below this tolerance times\n"
                           "              the parameter stationarity_tolerance_factor, then the algorithm\n"
                           "              terminates with a message of stationarity.\n"
                           "Default     : 1e-04.");
  options->addDoubleOption(reporter,
                           "stationarity_tolerance_factor",
                           1e+00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for checking termination with respect to stationarity.\n"
                           "              For further explanation, see description for the parameter\n"
                           "              stationarity_tolerance.\n"
                           "Default     : 1e+00.");

  // Add integer options

} // end addOptions

// Set options
void TerminationGradientCombination::setOptions(const Options* options,
                                                const Reporter* reporter)
{

  // Read bool options

  // Read double options
  options->valueAsDouble(reporter, "stationarity_tolerance", stationarity_tolerance_);
  options->valueAsDouble(reporter, "stationarity_tolerance_factor", stationarity_tolerance_factor_);

  // Read integer options

} // end setOptions

// Initialize
void TerminationGradientCombination::initialize(const Options* options,
                                                Quantities* quantities,
                                                const Reporter* reporter)
{

  // Set stationarity reference
  stationarity_reference_ = quantities->currentIterate()->gradient()->normInf();

} // end initialize

// Check final conditions
bool TerminationGradientCombination::checkFinal(const Options* options,
                                                Quantities* quantities,
                                                const Reporter* reporter,
                                                Strategies* strategies) const
{

  // Initialize check
  bool check = false;

  // Evaluate reference value
  double reference = fmax(1.0, fmax(stationarity_reference_, quantities->currentIterate()->gradient()->normInf()));

  // Check final conditions
  if (quantities->stationarityRadius() <= stationarity_tolerance_ &&
      strategies->qpSolver()->combinationTranslatedNormInf() <= stationarity_tolerance_ * stationarity_tolerance_factor_ * reference) {
    check = true;
  }

  // Print check values
  reporter->printf(R_NL, R_PER_ITERATION, " %+.2e %+.2e", quantities->currentIterate()->gradient()->normInf(), strategies->qpSolver()->combinationTranslatedNormInf());

  // Return
  return check;

} // end checkFinal

// Check radii not final condition
bool TerminationGradientCombination::checkRadiiNotFinal(const Options* options,
                                                        Quantities* quantities,
                                                        const Reporter* reporter,
                                                        Strategies* strategies) const
{

  // Initialize check
  bool check = false;

  // Check radii not final condition
  if (stationarity_tolerance_ < quantities->stationarityRadius()) {
    check = true;
  }

  // Return
  return check;
}

// Check radii update conditions
bool TerminationGradientCombination::checkRadiiUpdate(const Options* options,
                                                      Quantities* quantities,
                                                      const Reporter* reporter,
                                                      Strategies* strategies) const
{

  // Initialize check
  bool check = false;

  // Evaluate reference value
  double reference = fmax(1.0, fmax(stationarity_reference_, quantities->currentIterate()->gradient()->normInf()));

  // Check radii update conditions
  if (strategies->qpSolver()->primalSolutionNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference &&
      strategies->qpSolver()->combinationNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference &&
      strategies->qpSolver()->combinationTranslatedNormInf() <= quantities->stationarityRadius() * stationarity_tolerance_factor_ * reference) {
    check = true;
  }

  // Return
  return check;

} // end checkFinal

} // namespace NonOpt
