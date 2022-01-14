// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTVECTOR_HPP__
#define __TESTVECTOR_HPP__

#include <iostream>

#include "NonOptReporter.hpp"
#include "NonOptVector.hpp"

using namespace NonOpt;

// Implementation of test
int testVectorImplementation(int option)
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

  // Declare zero vector
  Vector u(5);

  // Set values to 0.0
  for (int i = 0; i < 5; i++) {
    u.values()[i] = 0.0;
  }

  // Check values
  for (int i = 0; i < 5; i++) {
    if (u.values()[i] < -1e-12 || u.values()[i] > 1e-12) {
      result = 1;
    }
  } // end for

  // Print ones vector
  u.print(&reporter, "Testing constructor... should be zero vector:");

  // Check length
  if (u.length() != 5) {
    result = 1;
  }

  // Print length
  reporter.printf(R_NL, R_BASIC, "Testing length access... print length (should be 5): %d\n", u.length());

  // Check first element
  if (u.values()[0] < -1e-12 || u.values()[0] > 1e-12) {
    result = 1;
  }

  // Print first element
  reporter.printf(R_NL, R_BASIC, "Testing value access... print first element (should be 0): %+23.16e\n", u.values()[0]);

  // Set first element to 1.0
  u.set(0, 1.0);

  // Check first element
  if (u.values()[0] < 1.0 - 1e-12 || u.values()[0] > 1.0 + 1e-12) {
    result = 1;
  }

  // Print first element
  u.print(&reporter, "Testing set method... now first element should be 1:");

  // Declare ones vector
  Vector v(5, 1.0);

  // Check values
  for (int i = 0; i < 5; i++) {
    if (v.values()[i] < 1.0 - 1e-12 || v.values()[i] > 1.0 + 1e-12) {
      result = 1;
    }
  } // end for

  // Print ones vector
  v.print(&reporter, "Testing constructor with default value... should be ones vector:");

  // Declare copy of ones vector
  std::shared_ptr<Vector> w = v.makeNewCopy();

  // Check values
  for (int i = 0; i < 5; i++) {
    if (w->values()[i] < 1.0 - 1e-12 || w->values()[i] > 1.0 + 1e-12) {
      result = 1;
    }
  } // end for

  // Print copy of ones vector
  w->print(&reporter, "Testing copy method... should be copy of ones vector:");

  // Declare linear combination
  std::shared_ptr<Vector> x = v.makeNewLinearCombination(2.0, 3.0, *w);

  // Check values
  for (int i = 0; i < 5; i++) {
    if (x->values()[i] < 5.0 - 1e-12 || x->values()[i] > 5.0 + 1e-12) {
      result = 1;
    }
  } // end for

  // Print linear combination
  x->print(&reporter, "Testing new linear combination method (multiple not 1)... should be vector of fives:");

  // Declare linear combination
  std::shared_ptr<Vector> y = v.makeNewLinearCombination(1.0, 7.0, u);

  // Check values
  if (y->values()[0] < 8.0 - 1e-12 || y->values()[0] > 8.0 + 1e-12) {
    result = 1;
  }
  for (int i = 1; i < 5; i++) {
    if (y->values()[i] < 1.0 - 1e-12 || y->values()[i] > 1.0 + 1e-12) {
      result = 1;
    }
  } // end for

  // Print linear combination
  y->print(&reporter, "Testing new linear combination method (multiple is 1)... should be [8,1,...,1]:");

  // Declare array
  double a[5] = {0.1, 0.2, 0.3, 0.4, 0.5};

  // Copy array
  y->copyArray(a);

  // Check values
  for (int i = 0; i < 5; i++) {
    if (i == 0) {
      if (y->values()[i] < 0.1 - 1e-12 || y->values()[i] > 0.1 + 1e-12) {
        result = 1;
      }
    }
    if (i == 1) {
      if (y->values()[i] < 0.2 - 1e-12 || y->values()[i] > 0.2 + 1e-12) {
        result = 1;
      }
    }
    if (i == 2) {
      if (y->values()[i] < 0.3 - 1e-12 || y->values()[i] > 0.3 + 1e-12) {
        result = 1;
      }
    }
    if (i == 3) {
      if (y->values()[i] < 0.4 - 1e-12 || y->values()[i] > 0.4 + 1e-12) {
        result = 1;
      }
    }
    if (i == 4) {
      if (y->values()[i] < 0.5 - 1e-12 || y->values()[i] > 0.5 + 1e-12) {
        result = 1;
      }
    }
  } // end for

  // Print array
  y->print(&reporter, "Testing copy array method... should be [0.1,0.2,0.3,0.4,0.5]:");

  // Declare new vector
  Vector z(5, 0.1);

  // Add scaled vector
  y->addScaledVector(2.0, z);

  // Check values
  for (int i = 0; i < 5; i++) {
    if (i == 0) {
      if (y->values()[i] < 0.3 - 1e-12 || y->values()[i] > 0.3 + 1e-12) {
        result = 1;
      }
    }
    if (i == 1) {
      if (y->values()[i] < 0.4 - 1e-12 || y->values()[i] > 0.4 + 1e-12) {
        result = 1;
      }
    }
    if (i == 2) {
      if (y->values()[i] < 0.5 - 1e-12 || y->values()[i] > 0.5 + 1e-12) {
        result = 1;
      }
    }
    if (i == 3) {
      if (y->values()[i] < 0.6 - 1e-12 || y->values()[i] > 0.6 + 1e-12) {
        result = 1;
      }
    }
    if (i == 4) {
      if (y->values()[i] < 0.7 - 1e-12 || y->values()[i] > 0.7 + 1e-12) {
        result = 1;
      }
    }
  } // end for

  // Print vector plus scaled vector
  y->print(&reporter, "Testing add scaled vector... should be [0.3,0.4,0.5,0.6,0.7]:");

  // Declare scaled vector
  w->scale(4.0);

  // Check values
  for (int i = 0; i < 5; i++) {
    if (w->values()[i] < 4.0 - 1e-12 || w->values()[i] > 4.0 + 1e-12) {
      result = 1;
    }
  } // end for

  // Print scaled vector
  w->print(&reporter, "Testing scale method... should be vector of fours:");

  // Compute inner product
  double wx = w->innerProduct(*x);

  // Check inner product
  if (wx < 100.0 - 1e-12 || wx > 100.0 + 1e-12) {
    result = 1;
  }

  // Print inner product
  reporter.printf(R_NL, R_BASIC, "Testing inner product method... should be 100: %+23.16e\n", wx);

  // Compute norms
  double w1 = w->norm1();
  double w2 = w->norm2();
  double wInf = w->normInf();

  // Check norms
  if (w1 < 20.0 - 1e-12 || w1 > 20.0 + 1e-12) {
    result = 1;
  }
  if (w2 < 8.944271909999 - 1e-12 || w2 > 8.944271909999 + 1e-12) {
    result = 1;
  }
  if (wInf < 4.0 - 1e-12 || wInf > 4.0 + 1e-12) {
    result = 1;
  }

  // Evaluate norms
  reporter.printf(R_NL, R_BASIC, "Testing norm methods... for preceding vector:\n"
                                 "1-norm   (should be 20                 ) : %+23.16e\n"
                                 "2-norm   (should be  8.9442719099991592) : %+23.16e\n"
                                 "inf-norm (should be  4                 ) : %+23.16e\n",
                  w1,
                  w2,
                  wInf);

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

} // end testVectorImplementation

#endif /* __TESTVECTOR_HPP__ */
