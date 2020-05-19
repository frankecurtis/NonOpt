// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptPointSetUpdateProximity.hpp"
#include "NonOptDefinitions.hpp"
#include <iostream>

namespace NonOpt
{

// Add options
void PointSetUpdateProximity::addOptions(Options* options,
                                         const Reporter* reporter)
{

  // Add double options
  options->addDoubleOption(reporter,
                           "PSP_envelope_factor",
                           1e+02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Magnitude of envelope to include around stationarity radius\n"
                           "when updating set.  If the difference between a point in the\n"
                           "set and the current iterate is larger in norm than this\n"
                           "radius times this envelope, then the point is removed;\n"
                           "otherwise, it is kept in the point set.\n"
                           "Default value: 1e+02.");
  options->addDoubleOption(reporter,
                           "PSP_size_factor",
                           1e+02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Size factor for removing points from point set.  If size of\n"
                           "point set exceeds this factor times the number of variables,\n"
                           "then the oldest members are removed.\n"
                           "Default value: 1e+02.");

  options->addDoubleOption(reporter,
                           "PSP_lower_size",
                           5e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Lower size factor for determine unadulterated BFGS, alpha_lower will be "
                           "PSP_lower_size*radius ");

  options->addDoubleOption(reporter,
                           "PSP_lower_epsilon",
                           1e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Lower size factor for determine unadulterated BFGS, alpha_lower will be "
                           "PSP_lower_size*radius ");

  // Add integer options

} // end addOptions

// Set options
void PointSetUpdateProximity::setOptions(const Options* options,
                                         const Reporter* reporter)
{

  // Read double options
  options->valueAsDouble(reporter, "PSP_envelope_factor", envelope_factor_);
  options->valueAsDouble(reporter, "PSP_size_factor", size_factor_);
  options->valueAsDouble(reporter, "PSP_lower_size", lower_size_);
  options->valueAsDouble(reporter, "PSP_lower_epsilon", lower_epsilon_);

  // Read integer options

} // end setOptions

// Initialize
void PointSetUpdateProximity::initialize(const Options* options,
                                         Quantities* quantities,
                                         const Reporter* reporter) {}

// Update point set
void PointSetUpdateProximity::updatePointSet(const Options* options,
                                             Quantities* quantities,
                                             const Reporter* reporter,
                                             Strategies* strategies)
{

  // Initialize status
  setStatus(PS_UNSET);

  // Check if Unadulterated BFGS, if is, then keep only current iterate in point set
  //  if(quantities->stepsize()> lower_size_*quantities->stationarityRadius() && lower_epsilon_*strategies->qpSolver()->primalSolutionNormInf()*strategies->qpSolver()->primalSolutionNormInf()<= strategies->qpSolver()->dualObjectiveQuadraticValue()){
  //	  quantities->pointSet()->clear();
  //	  quantities->pointSet()->push_back(quantities->currentIterate());
  //
  //  }

  // Remove old points
  while ((double)quantities->pointSet()->size() > std::min(5000.0, size_factor_ * (double)quantities->numberOfVariables())) {
    quantities->pointSet()->erase(quantities->pointSet()->begin());
  }

  // Loop through points
  for (int i = 0; i < (int)quantities->pointSet()->size(); i++) {

    // Create difference vector
    std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *(*quantities->pointSet())[i]->vector());

    // Check distance between point and current iterate
    if (difference->normInf() > std::max(2.0, envelope_factor_ / (double)quantities->numberOfVariables()) * quantities->stationarityRadius()) {

      // Erase element
      quantities->pointSet()->erase(quantities->pointSet()->begin() + i);

      // Decrement iterator
      i--;

    } // end if

  } // end for

  // Update status
  setStatus(PS_SUCCESS);

} // end updatePointSet

} // namespace NonOpt
