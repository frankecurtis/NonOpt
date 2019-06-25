// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTPOINT_HPP__
#define __TESTPOINT_HPP__

#include <iostream>

#include "MaxQ.hpp"
#include "NonOptOptions.hpp"
#include "NonOptPoint.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptRandomNumberGenerator.hpp"
#include "NonOptReporter.hpp"
#include "NonOptVector.hpp"

using namespace NonOpt;

// Implementation of test
int testPointImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter reporter;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> sr(new StreamReport("s", R_NL, R_BASIC));

    // Set stream report to standard output
    sr->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(sr);

  }  // end if

  // Declare problem
  std::shared_ptr<MaxQ> problem(new MaxQ);

  // Declare problem size
  int n;

  // Get number of variables
  problem->numberOfVariables(n);

  // Declare ones vector
  std::shared_ptr<Vector> v(new Vector(n, 1.0));

  // Declare point
  Point p(problem, v, 1.0);

  // Check values
  std::shared_ptr<Vector> vp = p.vector();
  for (int i = 0; i < n; i++) {
    if (vp->values()[i] < 1.0 - 1e-12 || vp->values()[i] > 1.0 + 1e-12) {
      result = 1;
    }
  }  // end for

  // Print point
  reporter.printf(R_NL, R_BASIC, "Testing constructor... should be ones vector:\n");
  p.print(&reporter, "OnesVector");

  // Declare quantities (for counter purposes)
  Quantities quantities;

  // Declare Options
  Options options;

  // Add and set quantities options
  quantities.addOptions(&options, &reporter);
  quantities.setOptions(&options, &reporter);

  // Get function value
  p.evaluateObjective(quantities);

  // Check objective value
  if (p.objective() < 1.0 - 1e-12 || p.objective() > 1.0 + 1e-12) {
    result = 1;
  }

  // Print function value
  reporter.printf(R_NL, R_BASIC, "Testing function evaluation... should be 1: %+23.16e\n", p.objective());

  // Get gradient value
  p.evaluateGradient(quantities);
  std::shared_ptr<Vector> g(p.gradient());

  // Check values
  for (int i = 0; i < n; i++) {
    if (i == 0) {
      if (g->values()[i] < 2.0 - 1e-12 || g->values()[i] > 2.0 + 1e-12) {
        result = 1;
      }
    }
    else {
      if (g->values()[i] < 0.0 - 1e-12 || g->values()[i] > 0.0 + 1e-12) {
        result = 1;
      }
    }
  }  // end for

  // Print gradient value
  reporter.printf(R_NL, R_BASIC, "Testing gradient evaluation... should be [2,0,...,0]:\n");
  g->print(&reporter, "GradientAtOnesVector");

  // Declare linear combination of Point
  std::shared_ptr<Point> r = p.makeNewLinearCombination(1.0, 2.0, *v);

  // Check values
  std::shared_ptr<Vector> rp = r->vector();
  for (int i = 0; i < n; i++) {
    if (rp->values()[i] < 3.0 - 1e-12 || rp->values()[i] > 3.0 + 1e-12) {
      result = 1;
    }
  }  // end for

  // Print linear combination of Point
  reporter.printf(R_NL, R_BASIC, "Testing new linear combination method... should be vector of threes:\n");
  r->print(&reporter, "LinearCombination");

  // Get function value
  r->evaluateObjective(quantities);

  // Check objective value
  if (r->objective() < 9.0 - 1e-12 || r->objective() > 9.0 + 1e-12) {
    result = 1;
  }

  // Print function value
  reporter.printf(R_NL, R_BASIC, "Testing function evaluation of linear combination... should be 9: %+23.16e\n", r->objective());

  // Get gradient value
  r->evaluateGradient(quantities);
  std::shared_ptr<Vector> gr(r->gradient());

  // Check values
  for (int i = 0; i < n; i++) {
    if (i == 0) {
      if (gr->values()[i] < 6.0 - 1e-12 || gr->values()[i] > 6.0 + 1e-12) {
        result = 1;
      }
    }
    else {
      if (gr->values()[i] < 0.0 - 1e-12 || gr->values()[i] > 0.0 + 1e-12) {
        result = 1;
      }
    }
  }  // end for

  // Print gradient value
  reporter.printf(R_NL, R_BASIC, "Testing gradient evaluation of linear combination... should be [6,0,...,0]:\n");
  gr->print(&reporter, "GradientAtLinearCombination");

  // Declare sampling radius
  double epsilon = 0.01;

  // Declare random new Point
  RandomNumberGenerator rng;
  std::shared_ptr<Point> s = p.makeNewRandom(epsilon, &rng);

  // Check values
  std::shared_ptr<Vector> sp = s->vector();
  for (int i = 0; i < n; i++) {
    if (sp->values()[i] < 1.0 - epsilon - 1e-12 || sp->values()[i] > 1.0 + epsilon + 1e-12) {
      result = 1;
    }
  }  // end for

  // Print random point
  reporter.printf(R_NL, R_BASIC, "Testing new random point method... should be close to ones vector:\n");
  s->print(&reporter, "RandomVector");

  // Get function value
  s->evaluateObjective(quantities);

  // Check objective value
  if (s->objective() < (1.0 - epsilon) * (1.0 - epsilon) - 1e-12 || s->objective() > (1.0 + epsilon) * (1.0 + epsilon) + 1e-12) {
    result = 1;
  }

  // Print function value
  reporter.printf(R_NL, R_BASIC, "Testing function evaluation of random point... should be near 1: %+23.16e\n", s->objective());

  // Get gradient value
  s->evaluateGradient(quantities);
  std::shared_ptr<Vector> gs(s->gradient());

  // Check values
  int index = 0;
  int value = 0.0;
  for (int i = 0; i < n; i++) {
    if (gs->values()[i] > value) {
      index = i;
      value = gs->values()[i];
    }
  }  // end for
  for (int i = 0; i < n; i++) {
    if (i == index) {
      if (gs->values()[i] < 2.0 * (1.0 - epsilon) - 1e-12 || gs->values()[i] > 2.0 * (1.0 + epsilon) + 1e-12) {
        result = 1;
      }
    }
    else {
      if (gs->values()[i] < 0.0 - 1e-12 || gs->values()[i] > 0.0 + 1e-12) {
        result = 1;
      }
    }
  }  // end for

  // Print gradient value
  reporter.printf(R_NL, R_BASIC, "Testing gradient evaluation of random point... should be near [0,...,0,2,0,...,0]:\n");
  gs->print(&reporter, "GradientAtRandomVector");

  // Print values again
  reporter.printf(R_NL, R_BASIC, "Print values again, in same order, to make sure none have changed:\n");
  reporter.printf(R_NL, R_BASIC, "... should be ones vector:\n");
  p.print(&reporter, "OnesVector");
  reporter.printf(R_NL, R_BASIC, "... should be 1: %+23.16e\n", p.objective());
  reporter.printf(R_NL, R_BASIC, "... should be [2,0,...,0]:\n");
  g->print(&reporter, "GradientAtOnesVector");
  reporter.printf(R_NL, R_BASIC, "... should be vector of threes:\n");
  r->print(&reporter, "LinearCombination");
  reporter.printf(R_NL, R_BASIC, "... should be 9: %+23.16e\n", r->objective());
  reporter.printf(R_NL, R_BASIC, "... should be [6,0,...,0]:\n");
  gr->print(&reporter, "GradientAtLinearCombination");
  reporter.printf(R_NL, R_BASIC, "... should be vector of near-ones:\n");
  s->print(&reporter, "RandomVector");
  reporter.printf(R_NL, R_BASIC, "... should be near 1: %+23.16e\n", s->objective());
  reporter.printf(R_NL, R_BASIC, "... should be near [0,...,0,2,0,...,0]:\n");
  gs->print(&reporter, "GradientAtRandomVector");

  // Return
  return result;

}  // end testPointImplementation

#endif /* __TESTPOINT_HPP__ */
