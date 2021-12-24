// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptPointSetUpdateProximity.hpp"
#include "NonOptDefinitions.hpp"

namespace NonOpt
{

// Add options
void PointSetUpdateProximity::addOptions(Options* options)
{

  // Add bool options

  // Add double options
  options->addDoubleOption("PSP_envelope_factor",
                           1e+02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Magnitude of envelope to include around stationarity radius\n"
                           "              when updating set.  If the difference between a point in the\n"
                           "              set and the current iterate is larger in norm than this\n"
                           "              radius times this envelope, then the point is removed;\n"
                           "              otherwise, it is kept in the point set.\n"
                           "Default     : 1e+02.");
  options->addDoubleOption("PSP_size_factor",
                           2e+00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Size factor for removing points from point set.  If size of\n"
                           "              point set exceeds this factor times the number of variables,\n"
                           "              then the oldest members are removed.\n"
                           "Default     : 2e+00.");

  // Add integer options
  options->addIntegerOption("PSP_size_maximum",
                           NONOPT_INT_INFINITY,
                           0.0,
                           NONOPT_INT_INFINITY,
                           "Size maximum for removing points from point set.  If size of\n"
                           "              point set exceeds this number, then the oldest members are\n"
                           "              removed.\n"
                           "Default     : Infinity.");

} // end addOptions

// Set options
void PointSetUpdateProximity::setOptions(Options* options)
{

  // Read bool options

  // Read double options
  options->valueAsDouble("PSP_envelope_factor", envelope_factor_);
  options->valueAsDouble("PSP_size_factor", size_factor_);

  // Read integer options
  options->valueAsInteger("PSP_size_maximum", size_maximum_);

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

  // Remove old points
  while ((double)quantities->pointSet()->size() > size_factor_ * (double)quantities->numberOfVariables() || quantities->pointSet()->size() > size_maximum_) {
    quantities->pointSet()->erase(quantities->pointSet()->begin());
  }

  // Loop through points
  for (int i = 0; i < (int)quantities->pointSet()->size(); i++) {

    // Create difference vector
    std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *(*quantities->pointSet())[i]->vector());

    // Check distance between point and current iterate
    if (difference->norm2() > envelope_factor_ * quantities->stationarityRadius()) {

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
