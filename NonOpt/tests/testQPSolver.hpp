// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Baoyu Zhou

#ifndef __TESTQPSOLVER_HPP__
#define __TESTQPSOLVER_HPP__

#include <iostream>

#include <cmath>
#include <random>

#include "NonOptQPSolverDualActiveSet.hpp"
#include "NonOptQPSolverInteriorPoint.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptSymmetricMatrixDense.hpp"

using namespace NonOpt;

// Implementation of test
int testQPSolverImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare test numbers
  int test_start = 0;
  int test_end = 20;
  int test_adds = 1;

  // Declare number of variables
  int numberVariables = 20;

  // Declare reporter
  Reporter reporter;
  Quantities quantities;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> s(new StreamReport("s", R_QP, R_BASIC));

    // Set stream report to standard output
    s->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(s);

  } // end if

  // Declare random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);

  // Declare options
  Options options;

  // Declare QP solver objects
  QPSolverDualActiveSet qDAS;
  QPSolverInteriorPoint qIPM;

  // Add options
  qDAS.addOptions(&options);
  qIPM.addOptions(&options);

  // Use exact solves
  options.modifyBoolValue("QPDAS_allow_inexact_termination", false);
  options.modifyBoolValue("QPIPM_allow_inexact_termination", false);
  options.modifyDoubleValue("QPDAS_kkt_tolerance", 1e-06);
  options.modifyDoubleValue("QPIPM_kkt_tolerance", 1e-06);
  options.modifyIntegerValue("QPDAS_iteration_limit", 1e+03);
  options.modifyIntegerValue("QPIPM_iteration_limit", 1e+03);

  // Set options
  qDAS.setOptions(&options);
  qIPM.setOptions(&options);

  // Initialize data
  qDAS.initializeData(numberVariables);
  qIPM.initializeData(numberVariables);

  // Loop over number of tests
  for (int test = test_start; test <= test_end; test++) {

    // Print test number
    reporter.printf(R_QP, R_BASIC, "Running test %8d... ", test);

    // Set random seed
    generator.seed(test);

    // Declare problem parameters
    int numberAffine = 2 * (test + 1);
    int numberActive = (test + 1);
    int numberPoints = 10 * (test + 1);
    int numberPointsAdd = 50;

    // Declare scaling factors
    double hess_scaling = 10;

    // Declare regularization
    double regularization = 1.0 / ((double)(test) + 1.0);

    /*
     * Set up generator function
     */

    // Initialize linear vector term
    Vector gen_maxLinearVector(numberAffine, 0.0);

    // Set random elements of linear vector term
    for (int i = numberActive; i < numberAffine; i++) {
      gen_maxLinearVector.set(i, -pow(uniform(generator), 2.0));
    }

    // Initialize and set linear matrix term
    double* gen_maxLinearMatrix = new double[numberAffine * numberVariables];
    for (int i = 0; i < numberAffine * numberVariables; i++) {
      gen_maxLinearMatrix[i] = normal(generator);
    }

    // Initialize weights
    Vector gen_weights(numberAffine, 0.0);

    // Set weights and normalize
    for (int i = 0; i < numberActive; i++) {
      gen_weights.set(i, uniform(generator));
    }
    gen_weights.scale(1.0 / gen_weights.norm1());

    // Initialize and set linear term
    Vector gen_linear(numberVariables);
    for (int i = 0; i < numberVariables; i++) {
      for (int j = 0; j < numberAffine; j++) {
        gen_linear.set(i, gen_linear.values()[i] - gen_maxLinearMatrix[j * numberVariables + i] * gen_weights.values()[j]);
      }
    } // end for

    // Initialize quadratic term
    double* gen_quadratic_init = new double[numberVariables * numberVariables];
    for (int i = 0; i < numberVariables * numberVariables; i++) {
      gen_quadratic_init[i] = normal(generator);
    }

    // Set quadratic term
    double* gen_quadratic = new double[numberVariables * numberVariables];
    for (int i = 0; i < numberVariables * numberVariables; i++) {
      gen_quadratic[i] = 0.0;
    }
    for (int i = 0; i < numberVariables; i++) {
      for (int j = 0; j < numberVariables; j++) {
        for (int k = 0; k < numberVariables; k++) {
          gen_quadratic[i * numberVariables + j] = gen_quadratic[i * numberVariables + j] + gen_quadratic_init[i * numberVariables + k] * gen_quadratic_init[j * numberVariables + k];
        }
      } // end for
    }   // end for
    for (int i = 0; i < numberVariables * numberVariables; i++) {
      gen_quadratic[i] = hess_scaling * gen_quadratic[i];
    }

    // Set center point
    Vector gen_center(numberVariables);
    for (int i = 0; i < numberVariables; i++) {
      gen_center.set(i, normal(generator));
    }

    /*
     * Set up problem data from generator function
     */

    // Initialize vector list
    std::vector<std::shared_ptr<Vector>> vector_list;
    std::vector<std::shared_ptr<Vector>> new_vector_list;

    // Initialize vector
    std::vector<double> vector;
    std::vector<double> new_vector;

    // Loop over number of points
    for (int points = 0; points < numberPoints + numberPointsAdd; points++) {

      // Declare point
      std::shared_ptr<Vector> point(new Vector(numberVariables));

      // Set point values
      for (int i = 0; i < numberVariables; i++) {
        point->set(i, gen_center.values()[i] + normal(generator));
      }

      // Declare affine vector
      Vector affine(numberAffine);

      // Set affine vector values
      for (int i = 0; i < numberAffine; i++) {
        double sum = 0.0;
        for (int j = 0; j < numberVariables; j++) {
          sum = sum + gen_maxLinearMatrix[i * numberVariables + j] * point->values()[j];
        }
        affine.set(i, sum + gen_maxLinearVector.values()[i]);
      } // end for

      // Initialize maximum of affine functions
      int index = 0;
      double maxAffine = affine.values()[0];

      // Determine maximum of affine functions
      for (int i = 1; i < numberAffine; i++) {
        if (affine.values()[i] > maxAffine) {
          index = i;
          maxAffine = affine.values()[i];
        }
      } // end for

      // Declare and set generator function value
      double gen_f = maxAffine + point->innerProduct(gen_linear);
      double* temp = new double[numberVariables];
      for (int i = 0; i < numberVariables; i++) {
        temp[i] = 0.0;
      }
      for (int i = 0; i < numberVariables; i++) {
        for (int j = 0; j < numberVariables; j++) {
          temp[i] = temp[i] + gen_quadratic[i * numberVariables + j] * point->values()[j];
        }
      } // end for
      for (int i = 0; i < numberVariables; i++) {
        gen_f = gen_f + 0.5 * point->values()[i] * temp[i];
      }

      // Declare and set generator gradient value
      std::shared_ptr<Vector> gen_g(new Vector(numberVariables));
      for (int i = 0; i < numberVariables; i++) {
        gen_g->set(i, temp[i] + gen_linear.values()[i] + gen_maxLinearMatrix[index * numberVariables + i]);
      }

      // Push vector into list
      if (points < numberPoints) {
        vector_list.push_back(gen_g);
      }
      else {
        new_vector_list.push_back(gen_g);
      }

      // Declare temporary scalar
      double junk = 0.0;
      for (int i = 0; i < numberVariables; i++) {
        junk = junk + gen_g->values()[i] * (gen_center.values()[i] - point->values()[i]);
      }
      if (points < numberPoints) {
        vector.push_back(gen_f + junk);
      }
      else {
        new_vector.push_back(gen_f + junk);
      }

      // Delete temporary vector
      delete[] temp;

    } // end for

    // Initialize inverse Hessian
    double* hessianInverse_init = new double[numberVariables * numberVariables];
    for (int i = 0; i < numberVariables * numberVariables; i++) {
      hessianInverse_init[i] = normal(generator);
    }

    // Set quadratic term
    std::shared_ptr<SymmetricMatrixDense> matrix = std::make_shared<SymmetricMatrixDense>();
    matrix->setAsDiagonal(numberVariables, 1.0);
    for (int i = 0; i < numberVariables; i++) {
      for (int j = 0; j < numberVariables; j++) {
        for (int k = 0; k < numberVariables; k++) {
          matrix->valuesOfInverseModifiable()[i * numberVariables + j] = matrix->valuesOfInverse()[i * numberVariables + j] + hessianInverse_init[i * numberVariables + k] * hessianInverse_init[j * numberVariables + k];
        }
      } // end for
    }   // end for
    for (int i = 0; i < numberVariables; i++) {
      for (int j = 0; j < numberVariables; j++) {
        matrix->valuesOfInverseModifiable()[i * numberVariables + j] = hess_scaling * matrix->valuesOfInverse()[i * numberVariables + j];
      }
    } // end for

    // Loop over cold and hot solves
    for (int solve_count = 0; solve_count < test_adds + 1; solve_count++) {

      // Loop over solvers
      for (int solver = 0; solver < 2; solver++) {

        // Check solve counter
        if (solve_count == 0 && solver == 0) {

          // Set QP data for dual active-set solver
          qDAS.setMatrix(matrix);
          qDAS.setVectorList(vector_list);
          qDAS.setVector(vector);
          qDAS.setScalar(regularization);

          // Solve QP
          qDAS.solveQP(&options, &reporter, &quantities);

        } // end if

        else if (solve_count == 1 && solver == 0) {

          // Add vectors
          qDAS.addData(new_vector_list, new_vector);

          // Solve QP hot
          qDAS.solveQPHot(&options, &reporter, &quantities);

          // Print new line
          reporter.printf(R_QP, R_BASIC, "... adding %2d vectors... ", numberPointsAdd);

        } // end else if

        else if (solve_count == 0 && solver == 1) {

          // Set QP data for interior-point solver
          qIPM.setMatrix(matrix);
          qIPM.setVectorList(vector_list);
          qIPM.setVector(vector);
          qIPM.setScalar(regularization);

          // Solve QP
          qIPM.solveQP(&options, &reporter, &quantities);

          // Print new line
          reporter.printf(R_QP, R_BASIC, "........................ ", numberPointsAdd);

        } // end else if

        else {

          // Add vectors
          qIPM.addData(new_vector_list, new_vector);

          // Solve QP (no hot start!)
          qIPM.solveQP(&options, &reporter, &quantities);

          // Print new line
          reporter.printf(R_QP, R_BASIC, "........................ ", numberPointsAdd);

        } // end else

        // Check for pass or fail
        if (solver == 0 && qDAS.status() == QP_SUCCESS) {
          reporter.printf(R_QP, R_BASIC, "DAS pass");
        }
        else if (solver == 0 && qDAS.status() != QP_SUCCESS) {
          reporter.printf(R_QP, R_BASIC, "DAS fail");
        }
        else if (solver == 1 && qIPM.status() == QP_SUCCESS) {
          reporter.printf(R_QP, R_BASIC, "IPM pass");
        }
        else {
          reporter.printf(R_QP, R_BASIC, "IPM fail");
        }

        // Declare primal solution
        Vector primal_solution(numberVariables);

        // Get primal solution
        if (solver == 0) {
          qDAS.primalSolution(primal_solution.valuesModifiable());
        }
        else {
          qIPM.primalSolution(primal_solution.valuesModifiable());
        }

        // Check results
        if (solver == 0 && (qDAS.status() != QP_SUCCESS || qDAS.KKTError() > 1e-03 || qDAS.KKTErrorDual() > 1e-03 || primal_solution.normInf() > 1.0 / ((double)(test) + 1.0) + 1e-03)) {
          result = 1;
        }
        else if (solver == 1 && (qIPM.status() != QP_SUCCESS || qIPM.KKTError() > 1e-03 || qIPM.KKTErrorDual() > 1e-03 || primal_solution.normInf() > 1.0 / ((double)(test) + 1.0) + 1e-03)) {
          result = 1;
        }

        // Print solve status information
        if (solver == 0) {
          reporter.printf(R_QP, R_BASIC, "  status: %d  iters: %6d  kkt error: %+.4e  kkt error (dual): %+.4e  ||step||_inf: %+.4e\n", qDAS.status(), qDAS.numberOfIterations(), qDAS.KKTError(), qDAS.KKTErrorDual(), primal_solution.normInf());
        }
        else {
          reporter.printf(R_QP, R_BASIC, "  status: %d  iters: %6d  kkt error: %+.4e  kkt error (dual): %+.4e  ||step||_inf: %+.4e\n", qIPM.status(), qIPM.numberOfIterations(), qIPM.KKTError(), qIPM.KKTErrorDual(), primal_solution.normInf());
        }

      } // end for

    } // end for

    // Delete matrix
    delete[] gen_maxLinearMatrix;
    delete[] gen_quadratic_init;
    delete[] gen_quadratic;
    delete[] hessianInverse_init;

  } // end for

  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      reporter.printf(R_QP, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    }
    else {
      reporter.printf(R_QP, R_BASIC, "TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testQPSolverImplementation

#endif /* __TESTQPSOLVER_HPP__ */
