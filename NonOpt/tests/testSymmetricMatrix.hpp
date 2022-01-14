// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Authors : Frank E. Curtis

#ifndef __TESTSYMMETRICMATRIX_HPP__
#define __TESTSYMMETRICMATRIX_HPP__

#include <iostream>

#include "MaxQ.hpp"
#include "NonOptOptions.hpp"
#include "NonOptProblem.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptSymmetricMatrixDense.hpp"
#include "NonOptSymmetricMatrixLimitedMemory.hpp"

using namespace NonOpt;

// Implementation of test
int testSymmetricMatrixImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare options
  Options options;

  // Declare problem
  std::shared_ptr<Problem> problem = std::make_shared<MaxQ>(10);

  // Declare quantities
  Quantities quantities;

  // Declare reporter
  Reporter reporter;

  // Initialize quantities
  quantities.initialize(problem);

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> sr(new StreamReport("s", R_NL, R_BASIC));

    // Set stream report to standard output
    sr->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(sr);

  } // end if

  // Declare pointer to symmetric matrix object
  std::shared_ptr<SymmetricMatrix> symmetric_matrix;
  symmetric_matrix = std::make_shared<SymmetricMatrixDense>();
  symmetric_matrix->addOptions(&options);
  symmetric_matrix = std::make_shared<SymmetricMatrixLimitedMemory>();
  symmetric_matrix->addOptions(&options);

  // Add approximate hessian update option
  options.addStringOption("approximate_hessian_update",
                          "BFGS",
                          "Approximate Hessian update strategy to use.\n"
                          "Default     : BFGS.");

  // Loop over symmetric matrix strategies
  for (int symmetric_matrix_number = 0; symmetric_matrix_number < 2; symmetric_matrix_number++) {

    // Loop over update strategies
    for (int approximate_hessian_update = 0; approximate_hessian_update < 2; approximate_hessian_update++) {

      // Set symmetric matrix
      if (symmetric_matrix_number == 0) {
        symmetric_matrix = std::make_shared<SymmetricMatrixDense>();
        reporter.printf(R_NL, R_BASIC, "TESTING SYMMETRIC MATRIX DENSE\n");
      } // end if
      else {
        symmetric_matrix = std::make_shared<SymmetricMatrixLimitedMemory>();
        reporter.printf(R_NL, R_BASIC, "TESTING SYMMETRIC MATRIX LIMITED-MEMORY\n");
      } // end else

      // Set update stratgy
      if (approximate_hessian_update == 0) {
        options.modifyStringValue("approximate_hessian_update", "BFGS");
        reporter.printf(R_NL, R_BASIC, "TESTING BFGS\n");
      } // end if
      else {
        options.modifyStringValue("approximate_hessian_update", "DFP");
        reporter.printf(R_NL, R_BASIC, "TESTING DFP\n");
      } // end else

      // Set options
      symmetric_matrix->setOptions(&options);

      // Initialize
      symmetric_matrix->initialize(&options, &quantities, &reporter);

      // Set size
      symmetric_matrix->setAsDiagonal(5, 1.0);

      // Use H as name
      std::shared_ptr<SymmetricMatrix> H = symmetric_matrix;

      // Check elements
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          if (i == j) {
            if (H->element(i, j) < 1.0 - 1e-08 || H->element(i, j) > 1.0 + 1e-08 ||
                H->elementOfInverse(i, j) < 1.0 - 1e-08 || H->elementOfInverse(i, j) > 1.0 + 1e-08) {
              result = 1;
            }
          }
          else {
            if (H->element(i, j) < 0.0 - 1e-08 || H->element(i, j) > 0.0 + 1e-08 ||
                H->elementOfInverse(i, j) < 0.0 - 1e-08 || H->elementOfInverse(i, j) > 0.0 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end for

      // Print identity matrix
      H->print(&reporter, "Testing constructor with default value... should be identity matrix:");

      // Declare ones vector
      Vector s(5, 1.0);

      // Print creation message
      s.print(&reporter, "Creating vector... should be ones vector:");

      // Declare [1,2,3,4,5] vector and set values
      Vector y(5, 1.0);
      y.set(1, 2.0);
      y.set(2, 3.0);
      y.set(3, 4.0);
      y.set(4, 5.0);

      // Print creation message
      y.print(&reporter, "Creating vector... should be [1,2,3,4,5]:");

      // Perform update
      H->update(s, y);

      // Print matrix
      H->print(&reporter, "Testing update... result with I, s=ones, y=1..5:");

      // Print element-by-element while checking values
      reporter.printf(R_NL, R_BASIC, "Testing element access... printing same matrices again:\n");
      reporter.printf(R_NL, R_BASIC, "Matrix:\n");
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->element(i, j));
            if (i == 0 && j == 0) {
              if (H->element(i, j) < 0.866666666666 - 1e-08 || H->element(i, j) > 0.866666666666 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->element(i, j) < -0.066666666666 - 1e-08 || H->element(i, j) > -0.066666666666 + 1e-08 ||
                  H->element(j, i) < -0.066666666666 - 1e-08 || H->element(j, i) > -0.066666666666 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->element(i, j) < 0.000000000000 - 1e-08 || H->element(i, j) > 0.000000000000 + 1e-08 ||
                  H->element(j, i) < 0.000000000000 - 1e-08 || H->element(j, i) > 0.000000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->element(i, j) < 0.066666666666 - 1e-08 || H->element(i, j) > 0.066666666666 + 1e-08 ||
                  H->element(j, i) < 0.066666666666 - 1e-08 || H->element(j, i) > 0.066666666666 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->element(i, j) < 0.133333333333 - 1e-08 || H->element(i, j) > 0.133333333333 + 1e-08 ||
                  H->element(j, i) < 0.133333333333 - 1e-08 || H->element(j, i) > 0.133333333333 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->element(i, j) < 1.066666666666 - 1e-08 || H->element(i, j) > 1.066666666666 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->element(i, j) < 0.200000000000 - 1e-08 || H->element(i, j) > 0.200000000000 + 1e-08 ||
                  H->element(j, i) < 0.200000000000 - 1e-08 || H->element(j, i) > 0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->element(i, j) < 0.333333333333 - 1e-08 || H->element(i, j) > 0.333333333333 + 1e-08 ||
                  H->element(j, i) < 0.333333333333 - 1e-08 || H->element(j, i) > 0.333333333333 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->element(i, j) < 0.466666666666 - 1e-08 || H->element(i, j) > 0.466666666666 + 1e-08 ||
                  H->element(j, i) < 0.466666666666 - 1e-08 || H->element(j, i) > 0.466666666666 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->element(i, j) < 1.400000000000 - 1e-08 || H->element(i, j) > 1.400000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->element(i, j) < 0.600000000000 - 1e-08 || H->element(i, j) > 0.600000000000 + 1e-08 ||
                  H->element(j, i) < 0.600000000000 - 1e-08 || H->element(j, i) > 0.600000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->element(i, j) < 0.800000000000 - 1e-08 || H->element(i, j) > 0.800000000000 + 1e-08 ||
                  H->element(j, i) < 0.800000000000 - 1e-08 || H->element(j, i) > 0.800000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->element(i, j) < 1.866666666666 - 1e-08 || H->element(i, j) > 1.866666666666 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->element(i, j) < 1.133333333333 - 1e-08 || H->element(i, j) > 1.133333333333 + 1e-08 ||
                  H->element(j, i) < 1.133333333333 - 1e-08 || H->element(j, i) > 1.133333333333 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->element(i, j) < 2.466666666666 - 1e-08 || H->element(i, j) > 2.466666666666 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->element(i, j));
            if (i == 0 && j == 0) {
              if (H->element(i, j) < 0.955555555555 - 1e-08 || H->element(i, j) > 0.955555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->element(i, j) < -0.022222222222 - 1e-08 || H->element(i, j) > -0.022222222222 + 1e-08 ||
                  H->element(j, i) < -0.022222222222 - 1e-08 || H->element(j, i) > -0.022222222222 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->element(i, j) < 0.000000000000 - 1e-08 || H->element(i, j) > 0.000000000000 + 1e-08 ||
                  H->element(j, i) < 0.000000000000 - 1e-08 || H->element(j, i) > 0.000000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->element(i, j) < 0.022222222222 - 1e-08 || H->element(i, j) > 0.022222222222 + 1e-08 ||
                  H->element(j, i) < 0.022222222222 - 1e-08 || H->element(j, i) > 0.022222222222 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->element(i, j) < 0.044444444444 - 1e-08 || H->element(i, j) > 0.044444444444 + 1e-08 ||
                  H->element(j, i) < 0.044444444444 - 1e-08 || H->element(j, i) > 0.044444444444 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->element(i, j) < 1.088888888888 - 1e-08 || H->element(i, j) > 1.088888888888 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->element(i, j) < 0.199999999999 - 1e-08 || H->element(i, j) > 0.199999999999 + 1e-08 ||
                  H->element(j, i) < 0.199999999999 - 1e-08 || H->element(j, i) > 0.199999999999 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->element(i, j) < 0.311111111111 - 1e-08 || H->element(i, j) > 0.311111111111 + 1e-08 ||
                  H->element(j, i) < 0.311111111111 - 1e-08 || H->element(j, i) > 0.311111111111 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->element(i, j) < 0.422222222222 - 1e-08 || H->element(i, j) > 0.422222222222 + 1e-08 ||
                  H->element(j, i) < 0.422222222222 - 1e-08 || H->element(j, i) > 0.422222222222 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->element(i, j) < 1.400000000000 - 1e-08 || H->element(i, j) > 1.400000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->element(i, j) < 0.599999999999 - 1e-08 || H->element(i, j) > 0.599999999999 + 1e-08 ||
                  H->element(j, i) < 0.599999999999 - 1e-08 || H->element(j, i) > 0.599999999999 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->element(i, j) < 0.799999999999 - 1e-08 || H->element(i, j) > 0.799999999999 + 1e-08 ||
                  H->element(j, i) < 0.799999999999 - 1e-08 || H->element(j, i) > 0.799999999999 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->element(i, j) < 1.888888888888 - 1e-08 || H->element(i, j) > 1.888888888888 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->element(i, j) < 1.177777777777 - 1e-08 || H->element(i, j) > 1.177777777777 + 1e-08 ||
                  H->element(j, i) < 1.177777777777 - 1e-08 || H->element(j, i) > 1.177777777777 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->element(i, j) < 2.555555555555 - 1e-08 || H->element(i, j) > 2.555555555555 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end else
      reporter.printf(R_NL, R_BASIC, "Matrix inverse:\n");
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->elementOfInverse(i, j));
            if (i == 0 && j == 0) {
              if (H->elementOfInverse(i, j) < 1.177777777777 - 1e-08 || H->elementOfInverse(i, j) > 1.177777777777 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.111111111111 - 1e-08 || H->elementOfInverse(i, j) > 0.111111111111 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.111111111111 - 1e-08 || H->elementOfInverse(j, i) > 0.111111111111 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.044444444444 - 1e-08 || H->elementOfInverse(i, j) > 0.044444444444 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.044444444444 - 1e-08 || H->elementOfInverse(j, i) > 0.044444444444 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.022222222222 - 1e-08 || H->elementOfInverse(i, j) > -0.022222222222 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.022222222222 - 1e-08 || H->elementOfInverse(j, i) > -0.022222222222 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.088888888888 - 1e-08 || H->elementOfInverse(i, j) > -0.088888888888 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.088888888888 - 1e-08 || H->elementOfInverse(j, i) > -0.088888888888 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->elementOfInverse(i, j) < 1.044444444444 - 1e-08 || H->elementOfInverse(i, j) > 1.044444444444 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.022222222222 - 1e-08 || H->elementOfInverse(i, j) > -0.022222222222 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.022222222222 - 1e-08 || H->elementOfInverse(j, i) > -0.022222222222 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.088888888888 - 1e-08 || H->elementOfInverse(i, j) > -0.088888888888 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.088888888888 - 1e-08 || H->elementOfInverse(j, i) > -0.088888888888 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.155555555555 - 1e-08 || H->elementOfInverse(i, j) > -0.155555555555 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.155555555555 - 1e-08 || H->elementOfInverse(j, i) > -0.155555555555 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->elementOfInverse(i, j) < 0.911111111111 - 1e-08 || H->elementOfInverse(i, j) > 0.911111111111 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.155555555555 - 1e-08 || H->elementOfInverse(i, j) > -0.155555555555 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.155555555555 - 1e-08 || H->elementOfInverse(j, i) > -0.155555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.222222222222 - 1e-08 || H->elementOfInverse(i, j) > -0.222222222222 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.222222222222 - 1e-08 || H->elementOfInverse(j, i) > -0.222222222222 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->elementOfInverse(i, j) < 0.777777777777 - 1e-08 || H->elementOfInverse(i, j) > 0.777777777777 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->elementOfInverse(i, j) < -0.288888888888 - 1e-08 || H->elementOfInverse(i, j) > -0.288888888888 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.288888888888 - 1e-08 || H->elementOfInverse(j, i) > -0.288888888888 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->elementOfInverse(i, j) < 0.644444444444 - 1e-08 || H->elementOfInverse(i, j) > 0.644444444444 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->elementOfInverse(i, j));
            if (i == 0 && j == 0) {
              if (H->elementOfInverse(i, j) < 1.048484848484 - 1e-08 || H->elementOfInverse(i, j) > 1.048484848484 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.030303030303 - 1e-08 || H->elementOfInverse(i, j) > 0.030303030303 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.030303030303 - 1e-08 || H->elementOfInverse(j, i) > 0.030303030303 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.012121212121 - 1e-08 || H->elementOfInverse(i, j) > 0.012121212121 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.012121212121 - 1e-08 || H->elementOfInverse(j, i) > 0.012121212121 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.006060606060 - 1e-08 || H->elementOfInverse(i, j) > -0.006060606060 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.006060606060 - 1e-08 || H->elementOfInverse(j, i) > -0.006060606060 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.024242424242 - 1e-08 || H->elementOfInverse(i, j) > -0.024242424242 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.024242424242 - 1e-08 || H->elementOfInverse(j, i) > -0.024242424242 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->elementOfInverse(i, j) < 0.993939393939 - 1e-08 || H->elementOfInverse(i, j) > 0.993939393939 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.042424242424 - 1e-08 || H->elementOfInverse(i, j) > -0.042424242424 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.042424242424 - 1e-08 || H->elementOfInverse(j, i) > -0.042424242424 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.078787878787 - 1e-08 || H->elementOfInverse(i, j) > -0.078787878787 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.078787878787 - 1e-08 || H->elementOfInverse(j, i) > -0.078787878787 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.115151515151 - 1e-08 || H->elementOfInverse(i, j) > -0.115151515151 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.115151515151 - 1e-08 || H->elementOfInverse(j, i) > -0.115151515151 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->elementOfInverse(i, j) < 0.903030303030 - 1e-08 || H->elementOfInverse(i, j) > 0.903030303030 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.151515151515 - 1e-08 || H->elementOfInverse(i, j) > -0.151515151515 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.151515151515 - 1e-08 || H->elementOfInverse(j, i) > -0.151515151515 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.206060606060 - 1e-08 || H->elementOfInverse(i, j) > -0.206060606060 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.206060606060 - 1e-08 || H->elementOfInverse(j, i) > -0.206060606060 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->elementOfInverse(i, j) < 0.775757575757 - 1e-08 || H->elementOfInverse(i, j) > 0.775757575757 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->elementOfInverse(i, j) < -0.296969696969 - 1e-08 || H->elementOfInverse(i, j) > -0.296969696969 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.296969696969 - 1e-08 || H->elementOfInverse(j, i) > -0.296969696969 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->elementOfInverse(i, j) < 0.612121212121 - 1e-08 || H->elementOfInverse(i, j) > 0.612121212121 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end else

      // Declare vector for column access
      Vector c(5, 0.0);

      // Get column
      H->column(3, c);

      // Check elements of column
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          if (i == 0) {
            if (c.values()[i] < 0.066666666666 - 1e-08 || c.values()[i] > 0.066666666666 + 1e-08) {
              result = 1;
            }
          }
          if (i == 1) {
            if (c.values()[i] < 0.333333333333 - 1e-08 || c.values()[i] > 0.333333333333 + 1e-08) {
              result = 1;
            }
          }
          if (i == 2) {
            if (c.values()[i] < 0.600000000000 - 1e-08 || c.values()[i] > 0.600000000000 + 1e-08) {
              result = 1;
            }
          }
          if (i == 3) {
            if (c.values()[i] < 1.866666666666 - 1e-08 || c.values()[i] > 1.866666666666 + 1e-08) {
              result = 1;
            }
          }
          if (i == 4) {
            if (c.values()[i] < 1.133333333333 - 1e-08 || c.values()[i] > 1.133333333333 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          if (i == 0) {
            if (c.values()[i] < 0.022222222222 - 1e-08 || c.values()[i] > 0.022222222222 + 1e-08) {
              result = 1;
            }
          }
          if (i == 1) {
            if (c.values()[i] < 0.311111111111 - 1e-08 || c.values()[i] > 0.311111111111 + 1e-08) {
              result = 1;
            }
          }
          if (i == 2) {
            if (c.values()[i] < 0.599999999999 - 1e-08 || c.values()[i] > 0.599999999999 + 1e-08) {
              result = 1;
            }
          }
          if (i == 3) {
            if (c.values()[i] < 1.888888888888 - 1e-08 || c.values()[i] > 1.888888888888 + 1e-08) {
              result = 1;
            }
          }
          if (i == 4) {
            if (c.values()[i] < 1.177777777777 - 1e-08 || c.values()[i] > 1.177777777777 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end else

      // Print column
      c.print(&reporter, "Column... should be vector of elements of column 3 of matrix:");

      // Get column
      H->columnOfInverse(3, c);

      // Check elements of column
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          if (i == 0) {
            if (c.values()[i] < -0.022222222222 - 1e-08 || c.values()[i] > -0.022222222222 + 1e-08) {
              result = 1;
            }
          }
          if (i == 1) {
            if (c.values()[i] < -0.088888888888 - 1e-08 || c.values()[i] > -0.088888888888 + 1e-08) {
              result = 1;
            }
          }
          if (i == 2) {
            if (c.values()[i] < -0.155555555555 - 1e-08 || c.values()[i] > -0.155555555555 + 1e-08) {
              result = 1;
            }
          }
          if (i == 3) {
            if (c.values()[i] < 0.777777777777 - 1e-08 || c.values()[i] > 0.777777777777 + 1e-08) {
              result = 1;
            }
          }
          if (i == 4) {
            if (c.values()[i] < -0.288888888888 - 1e-08 || c.values()[i] > -0.288888888888 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          if (i == 0) {
            if (c.values()[i] < -0.006060606060 - 1e-08 || c.values()[i] > -0.006060606060 + 1e-08) {
              result = 1;
            }
          }
          if (i == 1) {
            if (c.values()[i] < -0.078787878787 - 1e-08 || c.values()[i] > -0.078787878787 + 1e-08) {
              result = 1;
            }
          }
          if (i == 2) {
            if (c.values()[i] < -0.151515151515 - 1e-08 || c.values()[i] > -0.151515151515 + 1e-08) {
              result = 1;
            }
          }
          if (i == 3) {
            if (c.values()[i] < 0.775757575757 - 1e-08 || c.values()[i] > 0.775757575757 + 1e-08) {
              result = 1;
            }
          }
          if (i == 4) {
            if (c.values()[i] < -0.296969696969 - 1e-08 || c.values()[i] > -0.296969696969 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end else

      // Print column
      c.print(&reporter, "Column... should be vector of elements of column 3 of inverse:");

      // Declare vector for matrix-vector product
      Vector p(5, 0.0);

      // Compute matrix-vector product
      H->matrixVectorProduct(s, p);

      // Check elements of matrix-vector product
      for (int i = 0; i < 5; i++) {
        if (i == 0) {
          if (p.values()[i] < 1.000000000000 - 1e-08 || p.values()[i] > 1.000000000000 + 1e-08) {
            result = 1;
          }
        }
        if (i == 1) {
          if (p.values()[i] < 2.000000000000 - 1e-08 || p.values()[i] > 2.000000000000 + 1e-08) {
            result = 1;
          }
        }
        if (i == 2) {
          if (p.values()[i] < 3.000000000000 - 1e-08 || p.values()[i] > 3.000000000000 + 1e-08) {
            result = 1;
          }
        }
        if (i == 3) {
          if (p.values()[i] < 4.000000000000 - 1e-08 || p.values()[i] > 4.000000000000 + 1e-08) {
            result = 1;
          }
        }
        if (i == 4) {
          if (p.values()[i] < 5.000000000000 - 1e-08 || p.values()[i] > 5.000000000000 + 1e-08) {
            result = 1;
          }
        }
      } // end for

      // Print matrix-vector product
      p.print(&reporter, "Matrix-vector product... should be vector of sums of rows of matrix:");

      // Compute matrix-vector product
      H->matrixVectorProductOfInverse(s, p);

      // Check elements of matrix-vector product
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          if (i == 0) {
            if (p.values()[i] < 1.222222222222 - 1e-08 || p.values()[i] > 1.222222222222 + 1e-08) {
              result = 1;
            }
          }
          if (i == 1) {
            if (p.values()[i] < 0.888888888888 - 1e-08 || p.values()[i] > 0.888888888888 + 1e-08) {
              result = 1;
            }
          }
          if (i == 2) {
            if (p.values()[i] < 0.555555555555 - 1e-08 || p.values()[i] > 0.555555555555 + 1e-08) {
              result = 1;
            }
          }
          if (i == 3) {
            if (p.values()[i] < 0.222222222222 - 1e-08 || p.values()[i] > 0.222222222222 + 1e-08) {
              result = 1;
            }
          }
          if (i == 4) {
            if (p.values()[i] < -0.111111111111 - 1e-08 || p.values()[i] > -0.111111111111 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          if (i == 0) {
            if (p.values()[i] < 1.060606060606 - 1e-08 || p.values()[i] > 1.060606060606 + 1e-08) {
              result = 1;
            }
          }
          if (i == 1) {
            if (p.values()[i] < 0.787878787878 - 1e-08 || p.values()[i] > 0.787878787878 + 1e-08) {
              result = 1;
            }
          }
          if (i == 2) {
            if (p.values()[i] < 0.515151515151 - 1e-08 || p.values()[i] > 0.515151515151 + 1e-08) {
              result = 1;
            }
          }
          if (i == 3) {
            if (p.values()[i] < 0.242424242424 - 1e-08 || p.values()[i] > 0.242424242424 + 1e-08) {
              result = 1;
            }
          }
          if (i == 4) {
            if (p.values()[i] < -0.030303030303 - 1e-08 || p.values()[i] > -0.030303030303 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end else

      // Print matrix-vector product
      p.print(&reporter, "Matrix-vector product... should be vector of sums of rows of matrix inverse:");

      // Declare inner product
      double inner_product = H->innerProduct(s);

      // Check inner product
      if (inner_product < 15.000000000000 - 1e-08 || inner_product > 15.000000000000 + 1e-08) {
        result = 1;
      }

      // Compute inner product
      reporter.printf(R_NL, R_BASIC, "Inner product with matrix... should be 15.000000000000: %+23.16e\n", inner_product);

      // Declare inner product
      inner_product = H->innerProductOfInverse(s);

      // Check inner product
      if (approximate_hessian_update == 0) {
        if (inner_product < 2.777777777777 - 1e-08 || inner_product > 2.777777777777 + 1e-08) {
          result = 1;
        }
      }
      else {
        if (inner_product < 2.575757575757 - 1e-08 || inner_product > 2.575757575757 + 1e-08) {
          result = 1;
        }
      }

      // Compute inner product
      if (approximate_hessian_update == 0) {
        reporter.printf(R_NL, R_BASIC, "Inner product with matrix inverse... should be 2.777777777777: %+23.16e\n", inner_product);
      }
      else {
        reporter.printf(R_NL, R_BASIC, "Inner product with matrix inverse... should be 2.575757575757: %+23.16e\n", inner_product);
      }

      // Set vectors for second update
      s.set(0, -2.0);
      s.set(1, -1.0);
      s.set(2, 0.0);
      s.set(3, 1.0);
      s.set(4, 2.0);
      y.set(0, 0.0);
      y.set(1, 1.0);
      y.set(2, 2.0);
      y.set(3, 3.0);
      y.set(4, 4.0);

      // Print vector update message
      s.print(&reporter, "Updating vector... should be [-2,-1,0,1,2]:");

      // Print vector update message
      y.print(&reporter, "Updating vector... should be [0,1,2,3,4]:");

      // Perform update
      H->update(s, y);

      // Print matrix
      H->print(&reporter, "Testing update... result with H, s=-2..2, y=0..4:");

      // Print element-by-element while checking values
      reporter.printf(R_NL, R_BASIC, "Testing element access:\n");
      reporter.printf(R_NL, R_BASIC, "Matrix:\n");
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->element(i, j));
            if (i == 0 && j == 0) {
              if (H->element(i, j) < 0.760000000000 - 1e-08 || H->element(i, j) > 0.760000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->element(i, j) < -0.040000000000 - 1e-08 || H->element(i, j) > -0.040000000000 + 1e-08 ||
                  H->element(j, i) < -0.040000000000 - 1e-08 || H->element(j, i) > -0.040000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->element(i, j) < 0.160000000000 - 1e-08 || H->element(i, j) > 0.160000000000 + 1e-08 ||
                  H->element(j, i) < 0.160000000000 - 1e-08 || H->element(j, i) > 0.160000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->element(i, j) < 0.360000000000 - 1e-08 || H->element(i, j) > 0.360000000000 + 1e-08 ||
                  H->element(j, i) < 0.360000000000 - 1e-08 || H->element(j, i) > 0.360000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->element(i, j) < 0.560000000000 - 1e-08 || H->element(i, j) > 0.560000000000 + 1e-08 ||
                  H->element(j, i) < 0.560000000000 - 1e-08 || H->element(j, i) > 0.560000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->element(i, j) < 1.160000000000 - 1e-08 || H->element(i, j) > 1.160000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->element(i, j) < 0.360000000000 - 1e-08 || H->element(i, j) > 0.360000000000 + 1e-08 ||
                  H->element(j, i) < 0.360000000000 - 1e-08 || H->element(j, i) > 0.360000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->element(i, j) < 0.560000000000 - 1e-08 || H->element(i, j) > 0.560000000000 + 1e-08 ||
                  H->element(j, i) < 0.560000000000 - 1e-08 || H->element(j, i) > 0.560000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->element(i, j) < 0.760000000000 - 1e-08 || H->element(i, j) > 0.760000000000 + 1e-08 ||
                  H->element(j, i) < 0.760000000000 - 1e-08 || H->element(j, i) > 0.760000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->element(i, j) < 1.560000000000 - 1e-08 || H->element(i, j) > 1.560000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->element(i, j) < 0.760000000000 - 1e-08 || H->element(i, j) > 0.760000000000 + 1e-08 ||
                  H->element(j, i) < 0.760000000000 - 1e-08 || H->element(j, i) > 0.760000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->element(i, j) < 0.960000000000 - 1e-08 || H->element(i, j) > 0.960000000000 + 1e-08 ||
                  H->element(j, i) < 0.960000000000 - 1e-08 || H->element(j, i) > 0.960000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->element(i, j) < 1.960000000000 - 1e-08 || H->element(i, j) > 1.960000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->element(i, j) < 1.159999999999 - 1e-08 || H->element(i, j) > 1.159999999999 + 1e-08 ||
                  H->element(j, i) < 1.159999999999 - 1e-08 || H->element(j, i) > 1.159999999999 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->element(i, j) < 2.360000000000 - 1e-08 || H->element(i, j) > 2.360000000000 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->element(i, j));
            if (i == 0 && j == 0) {
              if (H->element(i, j) < 0.955555555555 - 1e-08 || H->element(i, j) > 0.955555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->element(i, j) < 0.155555555555 - 1e-08 || H->element(i, j) > 0.155555555555 + 1e-08 ||
                  H->element(j, i) < 0.155555555555 - 1e-08 || H->element(j, i) > 0.155555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->element(i, j) < 0.355555555555 - 1e-08 || H->element(i, j) > 0.355555555555 + 1e-08 ||
                  H->element(j, i) < 0.355555555555 - 1e-08 || H->element(j, i) > 0.355555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->element(i, j) < 0.555555555555 - 1e-08 || H->element(i, j) > 0.555555555555 + 1e-08 ||
                  H->element(j, i) < 0.555555555555 - 1e-08 || H->element(j, i) > 0.555555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->element(i, j) < 0.755555555555 - 1e-08 || H->element(i, j) > 0.755555555555 + 1e-08 ||
                  H->element(j, i) < 0.755555555555 - 1e-08 || H->element(j, i) > 0.755555555555 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->element(i, j) < 1.355555555555 - 1e-08 || H->element(i, j) > 1.355555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->element(i, j) < 0.555555555555 - 1e-08 || H->element(i, j) > 0.555555555555 + 1e-08 ||
                  H->element(j, i) < 0.555555555555 - 1e-08 || H->element(j, i) > 0.555555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->element(i, j) < 0.755555555555 - 1e-08 || H->element(i, j) > 0.755555555555 + 1e-08 ||
                  H->element(j, i) < 0.755555555555 - 1e-08 || H->element(j, i) > 0.755555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->element(i, j) < 0.955555555555 - 1e-08 || H->element(i, j) > 0.955555555555 + 1e-08 ||
                  H->element(j, i) < 0.955555555555 - 1e-08 || H->element(j, i) > 0.955555555555 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->element(i, j) < 1.755555555555 - 1e-08 || H->element(i, j) > 1.755555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->element(i, j) < 0.955555555555 - 1e-08 || H->element(i, j) > 0.955555555555 + 1e-08 ||
                  H->element(j, i) < 0.955555555555 - 1e-08 || H->element(j, i) > 0.955555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->element(i, j) < 1.155555555555 - 1e-08 || H->element(i, j) > 1.155555555555 + 1e-08 ||
                  H->element(j, i) < 1.155555555555 - 1e-08 || H->element(j, i) > 1.155555555555 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->element(i, j) < 2.155555555555 - 1e-08 || H->element(i, j) > 2.155555555555 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->element(i, j) < 1.355555555555 - 1e-08 || H->element(i, j) > 1.355555555555 + 1e-08 ||
                  H->element(j, i) < 1.355555555555 - 1e-08 || H->element(j, i) > 1.355555555555 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->element(i, j) < 2.555555555555 - 1e-08 || H->element(i, j) > 2.555555555555 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end else
      reporter.printf(R_NL, R_BASIC, "Matrix inverse:\n");
      if (approximate_hessian_update == 0) {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->elementOfInverse(i, j));
            if (i == 0 && j == 0) {
              if (H->elementOfInverse(i, j) < 1.800000000000 - 1e-08 || H->elementOfInverse(i, j) > 1.800000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.466666666666 - 1e-08 || H->elementOfInverse(i, j) > 0.466666666666 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.466666666666 - 1e-08 || H->elementOfInverse(j, i) > 0.466666666666 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.133333333333 - 1e-08 || H->elementOfInverse(i, j) > 0.133333333333 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.133333333333 - 1e-08 || H->elementOfInverse(j, i) > 0.133333333333 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.533333333333 - 1e-08 || H->elementOfInverse(i, j) > -0.533333333333 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.533333333333 - 1e-08 || H->elementOfInverse(j, i) > -0.533333333333 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->elementOfInverse(i, j) < 1.244444444444 - 1e-08 || H->elementOfInverse(i, j) > 1.244444444444 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->elementOfInverse(i, j) < 0.022222222222 - 1e-08 || H->elementOfInverse(i, j) > 0.022222222222 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.022222222222 - 1e-08 || H->elementOfInverse(j, i) > 0.022222222222 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.422222222222 - 1e-08 || H->elementOfInverse(i, j) > -0.422222222222 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.422222222222 - 1e-08 || H->elementOfInverse(j, i) > -0.422222222222 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->elementOfInverse(i, j) < 0.911111111111 - 1e-08 || H->elementOfInverse(i, j) > 0.911111111111 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.311111111111 - 1e-08 || H->elementOfInverse(i, j) > -0.311111111111 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.311111111111 - 1e-08 || H->elementOfInverse(j, i) > -0.311111111111 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->elementOfInverse(i, j) < 0.800000000000 - 1e-08 || H->elementOfInverse(i, j) > 0.800000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->elementOfInverse(i, j) < 0.911111111111 - 1e-08 || H->elementOfInverse(i, j) > 0.911111111111 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end if
      else {
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            reporter.printf(R_NL, R_BASIC, " %23.16e", H->elementOfInverse(i, j));
            if (i == 0 && j == 0) {
              if (H->elementOfInverse(i, j) < 1.448000000000 - 1e-08 || H->elementOfInverse(i, j) > 1.448000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.232000000000 - 1e-08 || H->elementOfInverse(i, j) > 0.232000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.232000000000 - 1e-08 || H->elementOfInverse(j, i) > 0.232000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 2) || (i == 2 && j == 0)) {
              if (H->elementOfInverse(i, j) < 0.016000000000 - 1e-08 || H->elementOfInverse(i, j) > 0.016000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < 0.016000000000 - 1e-08 || H->elementOfInverse(j, i) > 0.016000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 3) || (i == 3 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 0 && j == 4) || (i == 4 && j == 0)) {
              if (H->elementOfInverse(i, j) < -0.416000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.416000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.416000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.416000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 1 && j == 1) {
              if (H->elementOfInverse(i, j) < 1.088000000000 - 1e-08 || H->elementOfInverse(i, j) > 1.088000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.056000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.056000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.056000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.056000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 1 && j == 4) || (i == 4 && j == 1)) {
              if (H->elementOfInverse(i, j) < -0.344000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.344000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.344000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.344000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 2 && j == 2) {
              if (H->elementOfInverse(i, j) < 0.872000000000 - 1e-08 || H->elementOfInverse(i, j) > 0.872000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 2 && j == 4) || (i == 4 && j == 2)) {
              if (H->elementOfInverse(i, j) < -0.272000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.272000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.272000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.272000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 3 && j == 3) {
              if (H->elementOfInverse(i, j) < 0.799999999999 - 1e-08 || H->elementOfInverse(i, j) > 0.799999999999 + 1e-08) {
                result = 1;
              }
            }
            if ((i == 3 && j == 4) || (i == 4 && j == 3)) {
              if (H->elementOfInverse(i, j) < -0.200000000000 - 1e-08 || H->elementOfInverse(i, j) > -0.200000000000 + 1e-08 ||
                  H->elementOfInverse(j, i) < -0.200000000000 - 1e-08 || H->elementOfInverse(j, i) > -0.200000000000 + 1e-08) {
                result = 1;
              }
            }
            if (i == 4 && j == 4) {
              if (H->elementOfInverse(i, j) < 0.872000000000 - 1e-08 || H->elementOfInverse(i, j) > 0.872000000000 + 1e-08) {
                result = 1;
              }
            }
          } // end for
          reporter.printf(R_NL, R_BASIC, "\n");
        } // end for
      }   // end else

      // Reset as diagonal matrix
      H->setAsDiagonal(5, 7.0);

      // Check elements
      for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
          if (i == j) {
            if (H->element(i, j) < 7.0 - 1e-08 || H->element(i, j) > 7.0 + 1e-08 ||
                H->elementOfInverse(i, j) < 1.0 / 7.0 - 1e-08 || H->elementOfInverse(i, j) > 1.0 / 7.0 + 1e-08) {
              result = 1;
            }
          }
          else {
            if (H->element(i, j) < 0.0 - 1e-08 || H->element(i, j) > 0.0 + 1e-08 ||
                H->elementOfInverse(i, j) < 0.0 - 1e-08 || H->elementOfInverse(i, j) > 0.0 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end for

      // Print diagonal matrix
      H->print(&reporter, "Resetting to diagonal matrix... matrix should be 7*I:");

      // Set as 2x2
      H->setAsDiagonal(2, 1.0);

      // Check elements
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          if (i == j) {
            if (H->element(i, j) < 1.0 - 1e-08 || H->element(i, j) > 1.0 + 1e-08 ||
                H->elementOfInverse(i, j) < 1.0 - 1e-08 || H->elementOfInverse(i, j) > 1.0 + 1e-08) {
              result = 1;
            }
          }
          else {
            if (H->element(i, j) < 0.0 - 1e-08 || H->element(i, j) > 0.0 + 1e-08 ||
                H->elementOfInverse(i, j) < 0.0 - 1e-08 || H->elementOfInverse(i, j) > 0.0 + 1e-08) {
              result = 1;
            }
          }
        } // end for
      }   // end for

      // Print 2x2 matrix
      H->print(&reporter, "Resizing to 2-by-2 matrix... should be identity:");

    } // end for

  } // end for

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

} // end testSymmetricMatrixImplementation

#endif /* __TESTSYMMETRICMATRIX_HPP__ */
