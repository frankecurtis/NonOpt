// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Authors : Frank E. Curtis

#ifndef __TESTSYMMETRICMATRIX_HPP__
#define __TESTSYMMETRICMATRIX_HPP__

#include <iostream>
#include <iterator>

#include "NonOptOptions.hpp"
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

  // Declare pointer to symmetric matrix object
  std::shared_ptr<SymmetricMatrix> symmetric_matrix;
  symmetric_matrix = std::make_shared<SymmetricMatrixDense>();
  symmetric_matrix->addOptions(&options, &reporter);
  symmetric_matrix = std::make_shared<SymmetricMatrixLimitedMemory>();
  symmetric_matrix->addOptions(&options, &reporter);

  // Loop over strategies
  for (int symmetric_matrix_number = 0; symmetric_matrix_number < 1; symmetric_matrix_number++) {

    // Set symmetric matrix
    if (symmetric_matrix_number == 0) {
      symmetric_matrix = std::make_shared<SymmetricMatrixDense>();
      reporter.printf(R_NL, R_BASIC, "TESTING SYMMETRIC MATRIX DENSE\n");
    }  // end if
    else {
      symmetric_matrix = std::make_shared<SymmetricMatrixLimitedMemory>();
      reporter.printf(R_NL, R_BASIC, "TESTING SYMMETRIC MATRIX LIMITED-MEMORY\n");
    }  // end else

    // Set options
    symmetric_matrix->setOptions(&options, &reporter);

    // Set size
    symmetric_matrix->setAsDiagonal(5, 1.0);

    // Use I as name
    std::shared_ptr<SymmetricMatrix> I = symmetric_matrix;

    // Check elements
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        if (i == j) {
          if (I->elementHessian(i, j) < 1.0 - 1e-12 || I->elementHessian(i, j) > 1.0 + 1e-12) {
            result = 1;
          }
        }
        else {
          if (I->elementHessian(i, j) < 0.0 - 1e-12 || I->elementHessian(i, j) > 0.0 + 1e-12) {
            result = 1;
          }
        }
      }  // end for
    }    // end for

    // Print identity matrix

    //I->print(&reporter, "Testing constructor with default value... should be identity matrix:");

    // Declare ones vector
    Vector s(5, 1.0);

    // Print creation message
    //s.print(&reporter, "Creating vector... should be ones vector:");

    // Declare [1,2,3,4,5] vector and set values
    Vector y(5, 1.0);
    y.set(1, 2.0);
    y.set(2, 3.0);
    y.set(3, 4.0);
    y.set(4, 5.0);

    // Print creation message
    //y.print(&reporter, "Creating vector... should be [1,2,3,4,5]:");

    // Perform BFGS Update
    I->updateBFGS(s, y);

    // Print matrix
    //I->print(&reporter, "Testing Symmetric update... result with I, s=ones, y=1..5:");

    reporter.printf(R_NL, R_BASIC, "Testing element access... printing same matrix again:\n");
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        reporter.printf(R_NL, R_BASIC, " %23.16e", I->elementHessianInverse(i, j));
      }
      reporter.printf(R_NL, R_BASIC, "\n");
    }

    // Print element-by-element while checking values
    reporter.printf(R_NL, R_BASIC, "Testing elementHessian access... printing same matrix again:\n");
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        reporter.printf(R_NL, R_BASIC, " %23.16e", I->elementHessian(i, j));
        if (i == 0 && j == 0) {
          if (I->elementHessian(i, j) < 1.177777777777 - 1e-12 || I->elementHessian(i, j) > 1.177777777777 + 1e-12) {
            result = 1;
          }
        }
        if (i == 0 && j == 1) {
          if (I->elementHessian(i, j) < 0.111111111111 - 1e-12 || I->elementHessian(i, j) > 0.111111111111 + 1e-12) {
            result = 1;
          }
        }
        if (i == 0 && j == 2) {
          if (I->elementHessian(i, j) < 0.044444444444 - 1e-12 || I->elementHessian(i, j) > 0.044444444444 + 1e-12) {
            result = 1;
          }
        }
        if (i == 0 && j == 3) {
          if (I->elementHessian(i, j) < -0.022222222222 - 1e-12 || I->elementHessian(i, j) > -0.022222222222 + 1e-12) {
            result = 1;
          }
        }
        if (i == 0 && j == 4) {
          if (I->elementHessian(i, j) < -0.088888888888 - 1e-12 || I->elementHessian(i, j) > -0.088888888888 + 1e-12) {
            result = 1;
          }
        }
        if (i == 1 && j == 0) {
          if (I->elementHessian(i, j) < 0.111111111111 - 1e-12 || I->elementHessian(i, j) > 0.111111111111 + 1e-12) {
            result = 1;
          }
        }
        if (i == 1 && j == 1) {
          if (I->elementHessian(i, j) < 1.044444444444 - 1e-12 || I->elementHessian(i, j) > 1.044444444444 + 1e-12) {
            result = 1;
          }
        }
        if (i == 1 && j == 2) {
          if (I->elementHessian(i, j) < -0.022222222222 - 1e-12 || I->elementHessian(i, j) > -0.022222222222 + 1e-12) {
            result = 1;
          }
        }
        if (i == 1 && j == 3) {
          if (I->elementHessian(i, j) < -0.088888888888 - 1e-12 || I->elementHessian(i, j) > -0.088888888888 + 1e-12) {
            result = 1;
          }
        }
        if (i == 1 && j == 4) {
          if (I->elementHessian(i, j) < -0.155555555555 - 1e-12 || I->elementHessian(i, j) > -0.155555555555 + 1e-12) {
            result = 1;
          }
        }
        if (i == 2 && j == 0) {
          if (I->elementHessian(i, j) < 0.044444444444 - 1e-12 || I->elementHessian(i, j) > 0.044444444444 + 1e-12) {
            result = 1;
          }
        }
        if (i == 2 && j == 1) {
          if (I->elementHessian(i, j) < -0.022222222222 - 1e-12 || I->elementHessian(i, j) > -0.022222222222 + 1e-12) {
            result = 1;
          }
        }
        if (i == 2 && j == 2) {
          if (I->elementHessian(i, j) < 0.911111111111 - 1e-12 || I->elementHessian(i, j) > 0.911111111111 + 1e-12) {
            result = 1;
          }
        }
        if (i == 2 && j == 3) {
          if (I->elementHessian(i, j) < -0.155555555555 - 1e-12 || I->elementHessian(i, j) > -0.155555555555 + 1e-12) {
            result = 1;
          }
        }
        if (i == 2 && j == 4) {
          if (I->elementHessian(i, j) < -0.222222222222 - 1e-12 || I->elementHessian(i, j) > -0.222222222222 + 1e-12) {
            result = 1;
          }
        }
        if (i == 3 && j == 0) {
          if (I->elementHessian(i, j) < -0.022222222222 - 1e-12 || I->elementHessian(i, j) > -0.022222222222 + 1e-12) {
            result = 1;
          }
        }
        if (i == 3 && j == 1) {
          if (I->elementHessian(i, j) < -0.088888888888 - 1e-12 || I->elementHessian(i, j) > -0.088888888888 + 1e-12) {
            result = 1;
          }
        }
        if (i == 3 && j == 2) {
          if (I->elementHessian(i, j) < -0.155555555555 - 1e-12 || I->elementHessian(i, j) > -0.155555555555 + 1e-12) {
            result = 1;
          }
        }
        if (i == 3 && j == 3) {
          if (I->elementHessian(i, j) < 0.777777777777 - 1e-12 || I->elementHessian(i, j) > 0.777777777777 + 1e-12) {
            result = 1;
          }
        }
        if (i == 3 && j == 4) {
          if (I->elementHessian(i, j) < -0.288888888888 - 1e-12 || I->elementHessian(i, j) > -0.288888888888 + 1e-12) {
            result = 1;
          }
        }
        if (i == 4 && j == 0) {
          if (I->elementHessian(i, j) < -0.088888888888 - 1e-12 || I->elementHessian(i, j) > -0.088888888888 + 1e-12) {
            result = 1;
          }
        }
        if (i == 4 && j == 1) {
          if (I->elementHessian(i, j) < -0.155555555555 - 1e-12 || I->elementHessian(i, j) > -0.155555555555 + 1e-12) {
            result = 1;
          }
        }
        if (i == 4 && j == 2) {
          if (I->elementHessian(i, j) < -0.222222222222 - 1e-12 || I->elementHessian(i, j) > -0.222222222222 + 1e-12) {
            result = 1;
          }
        }
        if (i == 4 && j == 3) {
          if (I->elementHessian(i, j) < -0.288888888888 - 1e-12 || I->elementHessian(i, j) > -0.288888888888 + 1e-12) {
            result = 1;
          }
        }
        if (i == 4 && j == 4) {
          if (I->elementHessian(i, j) < 0.644444444444 - 1e-12 || I->elementHessian(i, j) > 0.644444444444 + 1e-12) {
            result = 1;
          }
        }
      }  // end for
      reporter.printf(R_NL, R_BASIC, "\n");
    }  // end for


    // Test Hessian multiply Inverse Hessian
    std::vector<double> res;
    for(int i=0;i<5;i++){
      Vector col(5, 0.0);
      I->column_Hessian(i, col);
      Vector prod(5, 0.0);
      I->matrixVectorProduct_HessianInverse(col,prod);

    for(int j=0;j<5;j++){
//      printf("%.4e",prod.values()[j]);
      res.push_back(prod.values()[j]);
    }
    }

    reporter.printf(R_NL, R_BASIC, "Testing Matrix Matrix multiplication... printing the product:\n");
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        reporter.printf(R_NL, R_BASIC, " %23.16e", res[i+5*j]);
        if (i == j) {
          if (res[i+5*j] < 1.0 - 1e-12 || res[i+5*j] > 1.0 + 1e-12) {
            result = 1;
          }
        }
        else {
          if (res[i+5*j] < 0.0 - 1e-12 || res[i+5*j] > 0.0 + 1e-12) {
            result = 1;
          }
        }
      }  // end for
      reporter.printf(R_NL, R_BASIC, "\n");
    }    // end for







    // Declare vector for column access
    Vector c(5, 0.0);

    // Get column
    I->column(3, c);

    // Check elements of column
    for (int i = 0; i < 5; i++) {
      if (i == 0) {
        if (c.values()[i] < -0.022222222222 - 1e-12 || c.values()[i] > -0.022222222222 + 1e-12) {
          result = 1;
        }
      }
      if (i == 1) {
        if (c.values()[i] < -0.088888888888 - 1e-12 || c.values()[i] > -0.088888888888 + 1e-12) {
          result = 1;
        }
      }
      if (i == 2) {
        if (c.values()[i] < -0.155555555555 - 1e-12 || c.values()[i] > -0.155555555555 + 1e-12) {
          result = 1;
        }
      }
      if (i == 3) {
        if (c.values()[i] < 0.777777777777 - 1e-12 || c.values()[i] > 0.777777777777 + 1e-12) {
          result = 1;
        }
      }
      if (i == 4) {
        if (c.values()[i] < -0.288888888888 - 1e-12 || c.values()[i] > -0.288888888888 + 1e-12) {
          result = 1;
        }
      }
    }  // end for

    // Print column
    c.print(&reporter, "Column... should be vector of elements of column 3 of I:");

    // Declare vector for matrix-vector product
    Vector p(5, 0.0);

    // Compute matrix-vector product
    I->matrixVectorProduct_HessianInverse(s, p);

    // Check elements of matrix-vector product
    for (int i = 0; i < 5; i++) {
      if (i == 0) {
        if (p.values()[i] < 1.222222222222 - 1e-12 || p.values()[i] > 1.222222222222 + 1e-12) {
          result = 1;
        }
      }
      if (i == 1) {
        if (p.values()[i] < 0.888888888888 - 1e-12 || p.values()[i] > 0.888888888888 + 1e-12) {
          result = 1;
        }
      }
      if (i == 2) {
        if (p.values()[i] < 0.555555555555 - 1e-12 || p.values()[i] > 0.555555555555 + 1e-12) {
          result = 1;
        }
      }
      if (i == 3) {
        if (p.values()[i] < 0.222222222222 - 1e-12 || p.values()[i] > 0.222222222222 + 1e-12) {
          result = 1;
        }
      }
      if (i == 4) {
        if (p.values()[i] < -0.111111111111 - 1e-12 || p.values()[i] > -0.111111111111 + 1e-12) {
          result = 1;
        }
      }
    }  // end for

    // Print matrix-vector product
    p.print(&reporter, "Matrix-vector product... should be vector of sums of rows of I:");

    // Declare inner product
    double inner_product = I->innerProduct_HessianInverse(s);

    // Check inner product
    if (inner_product < 2.777777777777 - 1e-12 || inner_product > 2.777777777777 + 1e-12) {
      result = 1;
    }

    // Compute inner product
    reporter.printf(R_NL, R_BASIC, "Inner product... should be 2.7777777777777772: %+23.16e\n", inner_product);

    // Reset as diagonal matrix
    I->setAsDiagonal(5, 7.0);

    // Check elements
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        if (i == j) {
          if (I->elementHessian(i, j) < 7.0 - 1e-12 || I->elementHessian(i, j) > 7.0 + 1e-12) {
            result = 1;
          }
        }
        else {
          if (I->elementHessian(i, j) < 0.0 - 1e-12 || I->elementHessian(i, j) > 0.0 + 1e-12) {
            result = 1;
          }
        }
      }  // end for
    }    // end for

    // Print diagonal matrix
    I->print(&reporter, "Resetting to diagonal matrix... should be 7*I:");

    // Set as 2x2
    I->setAsDiagonal(2, 1.0);

    // Check elements
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        if (i == j) {
          if (I->elementHessian(i, j) < 1.0 - 1e-12 || I->elementHessian(i, j) > 1.0 + 1e-12) {
            result = 1;
          }
        }
        else {
          if (I->elementHessian(i, j) < 0.0 - 1e-12 || I->elementHessian(i, j) > 0.0 + 1e-12) {
            result = 1;
          }
        }
      }  // end for
    }    // end for

    // Print 2x2 matrix
    I->print(&reporter, "Resizing to 2-by-2 matrix... should be identity:");

  }  // end for

  // Return
  return result;

}  // end testSymmetricMatrixImplementation

#endif /* __TESTSYMMETRICMATRIX_HPP__ */
