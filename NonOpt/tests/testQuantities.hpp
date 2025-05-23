// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTITERATIONQUANTITIES_HPP__
#define __TESTITERATIONQUANTITIES_HPP__

#include <iostream>

#include "MaxQ.hpp"
#include "NonOptProblem.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptVector.hpp"

using namespace NonOpt;

// Main function
int testQuantitiesImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter reporter;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> s(new StreamReport("s", R_NL, R_BASIC));

    // Set stream report to standard output
    s->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(s);

  } // end if

  // Declare problem
  std::shared_ptr<MaxQ> problem(new MaxQ(50));

  // Declare quantities
  Quantities quantities;

  // Declare Options
  Options options;

  // Add and set quantities options
  quantities.addOptions(&options);
  quantities.setOptions(&options);

  // Initialize quantities
  quantities.initialize(problem);

  // Evaluate gradient (for radii initialization)
  quantities.currentIterate()->evaluateGradient(quantities);

  // Check quantities
  if (quantities.numberOfVariables() != 50) {
    result = 1;
  }
  for (int i = 0; i < 50; i++) {
    if (i < 25) {
      if (quantities.currentIterate()->vector()->values()[i] < i + 1 - 1e-12 || quantities.currentIterate()->vector()->values()[i] > i + 1 + 1e-12) {
        result = 1;
      }
    }
    else {
      if (quantities.currentIterate()->vector()->values()[i] < -i - 1 - 1e-12 || quantities.currentIterate()->vector()->values()[i] > -i - 1 + 1e-12) {
        result = 1;
      }
    }
  } // end for
  if (quantities.stationarityRadius() < 1e+01 - 1e-12 || quantities.stationarityRadius() > 1e+01 + 1e-12) {
    result = 1;
  }
  if (quantities.trustRegionRadius() < 1e+12 - 1e-12 || quantities.trustRegionRadius() > 1e+12 + 1e-12) {
    result = 1;
  }

  // Print quantities
  reporter.printf(R_NL, R_BASIC, "Problem size: %d\n", quantities.numberOfVariables());
  quantities.currentIterate()->print(&reporter, "Current iterate");
  reporter.printf(R_NL, R_BASIC, "Stationarity radius: %+23.16e\n", quantities.stationarityRadius());
  reporter.printf(R_NL, R_BASIC, "Trust region radius: %+23.16e\n", quantities.trustRegionRadius());

  ///////////////////////
  // Update quantities //
  ///////////////////////

  // Declare problem size
  int n;

  // Get number of variables
  problem->numberOfVariables(n);

  // Declare ones vector
  std::shared_ptr<Vector> v(new Vector(n, 99.0));

  // Declare point
  std::shared_ptr<Point> p(new Point(problem, v, 1.0));

  // Set new iterate
  quantities.setCurrentIterate(p);

  // Modify radii
  quantities.updateRadii();

  // Check quantities
  if (quantities.numberOfVariables() != 50) {
    result = 1;
  }
  for (int i = 0; i < 50; i++) {
    if (quantities.currentIterate()->vector()->values()[i] < 99.0 - 1e-12 || quantities.currentIterate()->vector()->values()[i] > 99.0 + 1e-12) {
      result = 1;
    }
  } // end for
  if (quantities.stationarityRadius() < 1e+00 - 1e-12 || quantities.stationarityRadius() > 1e+00 + 1e-12) {
    result = 1;
  }
  if (quantities.trustRegionRadius() < 1e+11 - 1e-12 || quantities.trustRegionRadius() > 1e+11 + 1e-12) {
    result = 1;
  }

  // Print quantities
  quantities.currentIterate()->print(&reporter, "Current iterate (updated)");
  reporter.printf(R_NL, R_BASIC, "Stationarity radius (updated): %+23.16e\n", quantities.stationarityRadius());
  reporter.printf(R_NL, R_BASIC, "Trust region radius (updated): %+23.16e\n", quantities.trustRegionRadius());

  //////////////////////
  // Reset quantities //
  //////////////////////

  // Reset quantities
  quantities.initialize(problem);

  // Check quantities
  if (quantities.numberOfVariables() != 50) {
    result = 1;
  }
  for (int i = 0; i < 50; i++) {
    if (i < 25) {
      if (quantities.currentIterate()->vector()->values()[i] < i + 1 - 1e-12 || quantities.currentIterate()->vector()->values()[i] > i + 1 + 1e-12) {
        result = 1;
      }
    }
    else {
      if (quantities.currentIterate()->vector()->values()[i] < -i - 1 - 1e-12 || quantities.currentIterate()->vector()->values()[i] > -i - 1 + 1e-12) {
        result = 1;
      }
    }
  } // end for
  if (quantities.stationarityRadius() < 1e+01 - 1e-12 || quantities.stationarityRadius() > 1e+01 + 1e-12) {
    result = 1;
  }
  if (quantities.trustRegionRadius() < 1e+12 - 1e-12 || quantities.trustRegionRadius() > 1e+12 + 1e-12) {
    result = 1;
  }

  // Print quantities
  quantities.currentIterate()->print(&reporter, "Current iterate (reset)");
  reporter.printf(R_NL, R_BASIC, "Stationarity radius (reset): %+23.16e\n", quantities.stationarityRadius());
  reporter.printf(R_NL, R_BASIC, "Trust region radius (reset): %+23.16e\n", quantities.trustRegionRadius());

  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      reporter.printf(R_NL, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    }
    else {
      reporter.printf(R_NL, R_BASIC, "TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testQuantitiesImplementation

#endif /* __TESTITERATIONQUANTITIES_HPP__ */
