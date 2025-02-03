// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Lara Zebiane

#ifndef __TESTQPSOLVER_HPP__
#define __TESTQPSOLVER_HPP__

#include <iomanip>
#include <iostream>

#include <cmath>
#include <random>

#include "NonOptDefinitions.hpp"
#include "NonOptQPSolverDualActiveSet.hpp"
#include "NonOptQPSolverInteriorPoint.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptSymmetricMatrixDense.hpp"

using namespace NonOpt;

// Implementation of test
int testQPSolverImplementation(int option)
{

  // Declare test numbers
  int test_start = 1;
  int test_end = 1;

  // Declare n sizes
  int n_start = 0;
  int n_end = 0;
  int n_increment = 200;

  // Declare m sizes
  int m_sizes = 3;
  std::array<double, 3> m_factors = {1.0, 1.5, 2.0};
  std::array<int, 3> m_additions = {1, 0, 0};

  // Declare time array
  int p1 = 3; // d_type values
  int p2 = 1; // n values
  int p3 = 3; // m values
  int p4 = 2; // algorithms
  int p5 = 1; // runs
  double time_array[3][1][3][2][1];
  for (int i = 0; i < p1; i++) {
    for (int j = 0; j < p2; j++) {
      for (int k = 0; k < p3; k++) {
        for (int l = 0; l < p4; l++) {
          for (int z = 0; z < p5; z++) {
            time_array[i][j][k][l][z] = 0.0;
          }
        }
      }
    }
  }

  // Declare reporter
  Reporter reporter;
  Quantities quantities;

  // Check option
  if (option == 1)
  {

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
  qIPM.addOptions(&options);
  qDAS.addOptions(&options);

  // Set options
  options.modifyBoolValue("QPDAS_allow_inexact_termination", false);
  options.modifyBoolValue("QPIPM_allow_inexact_termination", false);
  options.modifyDoubleValue("QPDAS_kkt_tolerance", 1e-06);
  options.modifyDoubleValue("QPIPM_kkt_tolerance", 1e-06);

  // Set options
  qIPM.setOptions(&options);
  qDAS.setOptions(&options);

  // Loop over solution types (d = 0, half at bounds, or all at bounds)
  for (int d_type = 0; d_type <= 2; d_type++) {

    // Loop over n
    for (int n_index = n_start; n_index <= n_end; n_index++) {

      // Set n
      int n = (n_index + 1) * n_increment;

      // Loop over m
      for (int m_index = 0; m_index < m_sizes; m_index++) {

        // Set m
        int m = (int)(m_factors[m_index] * n) + m_additions[m_index];

        // Print header line
        printf("d_type = %6d\n", d_type);
        printf("n      = %6d\n", n);
        printf("m      = %6d\n", m);
        printf("%11s  %11s  %11s  %6s  %11s  %11s  %11s  %6s\n",
               "KKT DAS",
               "d err DAS",
               "time DAS",
               "k DAS",
               "KKT IPM",
               "d err IPM",
               "time IPM",
               "k IPM");

        // Loop over number of tests
        for (int test = test_start; test <= test_end; test++) {

          // Set random seed
          generator.seed(test);

          // Initialize data
          qIPM.initializeData(n);
          qDAS.initializeData(n);

          // Set trust region radius
          double scalar = 1.0;

          // Initialize solution, unconstrained solution, and ones vector
          Vector d_star(n, 0.0);
          Vector d_unc(n, 0.0);
          Vector half(n, 0.0);
          Vector ones(n, 1.0);
          for (int i = 0; i < n / 2; i++) {
            half.set(i, 1.0);
          }

          // Check solution type
          if (d_type == 0) {
            // do nothing, d_star and d_unc already set
          }
          else if (d_type == 1) {
            d_star.addScaledVector(scalar, half);
            d_unc.addScaledVector(2.0 * scalar, half);
          }
          else {
            d_star.addScaledVector(scalar, ones);
            d_unc.addScaledVector(2.0 * scalar, ones);
          }

          // Set t
          Vector t(n, 0.0);
          t.addScaledVector(-1.0, d_unc);

          // Set matrix
          std::shared_ptr<SymmetricMatrixDense> matrix = std::make_shared<SymmetricMatrixDense>();
          matrix->setAsDiagonal(n, 1.0);

          // Initialize G
          std::vector<std::shared_ptr<Vector>> G;

          // Randomly generate G of size n by n-1
          for (int i = 0; i < n - 1; i++) {

            // Create random vector
            std::shared_ptr<Vector> G_col = std::make_shared<Vector>(n);
            for (int j = 0; j < n; j++) {
              G_col->set(j, normal(generator));
            }

            // Add it to G
            G.push_back(G_col);

          } // end for

          // Generate "omega hat"
          Vector omega_hat(n - 1, 0.0);
          for (int i = 0; i < n - 1; i++) {
            omega_hat.set(i, uniform(generator));
          }

          // Compute G*omega_hat
          Vector G_omega_hat(n, 0.0);
          for (int i = 0; i < n - 1; i++) {
            G_omega_hat.addScaledVector(omega_hat.values()[i], *G[i].get());
          }

          // Compute g = t - G*omega_hat
          std::shared_ptr<Vector> g = std::make_shared<Vector>(n, 0.0);
          g->linearCombination(1.0, t, -1.0, G_omega_hat);

          // Add g to G
          G.push_back(g);

          // Set "omega bar"
          Vector omega_bar(n, 0.0);
          for (int i = 0; i < n - 1; i++) {
            omega_bar.set(i, omega_hat.values()[i]);
          }
          omega_bar.set(n - 1, 1.0);

          // Set scaling
          double scaling = omega_bar.norm1();

          // Scale omega_bar
          omega_bar.scale(1.0 / scaling);

          // Scale all vectors in G
          for (auto &vec_ptr : G) {
            vec_ptr->scale(scaling);
          }

          // Extend G with m - n random vectors
          for (int i = 0; i < m - n; i++) {

            // Create random vector
            std::shared_ptr<Vector> G_col = std::make_shared<Vector>(n);
            for (int j = 0; j < n; j++) {
              G_col->set(j, normal(generator));
            }

            // Add it to G
            G.push_back(G_col);

          } // end for

          // Extend omega_bar to omega_star (by padding with zeros)
          Vector omega_star(m, 0.0);
          for (int i = 0; i < n; i++) {
            omega_star.set(i, omega_bar.values()[i]);
          }

          // Compute G*omega_star
          Vector G_omega_star(n, 0.0);
          for (int i = 0; i < m; i++) {
            G_omega_star.addScaledVector(omega_star.values()[i], *G[i].get());
          }

          // Compute q
          Vector q(n, 0.0);
          q.linearCombination(-1.0, d_star, -1.0, G_omega_star);

          // Initialize sigma and rho
          Vector sigma(n, 0.0);
          Vector rho(n, 0.0);

          // Set sigma and rho
          for (int i = 0; i < n; i++) {
            if (q.values()[i] > 0.0) {
              sigma.set(i, q.values()[i]);
            }
            else if (q.values()[i] < 0.0) {
              rho.set(i, -q.values()[i]);
            }
          } // end for

          // Set u
          double u = 5.0;

          // Generate v_omega
          Vector v_omega(m, 0.0);
          for (int i = 0; i < m; i++) {
            if (omega_star.values()[i] == 0.0) {
              v_omega.set(i, abs(uniform(generator)));
            }
          }

          // Compute d
          Vector d(n, 0.0);
          d.linearCombination(1.0, G_omega_star, 1.0, sigma);
          d.addScaledVector(-1.0, rho);

          // Compute G'*d
          Vector Gtd(m, 0.0);
          for (int i = 0; i < m; i++) {
            Gtd.values()[i] = G[i]->innerProduct(d);
          }

          // Compute b
          std::vector<double> b(m, 0.0);
          for (int i = 0; i < m; i++) {
            b[i] += Gtd.values()[i] + u - v_omega.values()[i];
          }

          // Set dual active-set solver data
          qDAS.setMatrix(matrix);
          qDAS.setVectorList(G);
          qDAS.setVector(b);
          qDAS.setScalar(scalar);

          // Set interior-point solver data
          qIPM.setMatrix(matrix);
          qIPM.setVectorList(G);
          qIPM.setVector(b);
          qIPM.setScalar(scalar);

          // Solve with dual active-set method
          std::clock_t start_DAS = std::clock();
          qDAS.solveQP(&options, &reporter, &quantities);
          std::clock_t end_DAS = std::clock();
          double time_DAS = static_cast<double>(end_DAS - start_DAS) / CLOCKS_PER_SEC;

          // Solve with interior-point method
          std::clock_t start_IPM = std::clock();
          qIPM.solveQP(&options, &reporter, &quantities);
          std::clock_t end_IPM = std::clock();
          double time_IPM = static_cast<double>(end_IPM - start_IPM) / CLOCKS_PER_SEC;

          // Get primal solution
          Vector primal_solution_DAS(n);
          qDAS.primalSolution(primal_solution_DAS.valuesModifiable());
          Vector vector_d_dstar_diff_DAS(n);
          vector_d_dstar_diff_DAS.linearCombination(1, d_star, -1, primal_solution_DAS);

          // Get primal solution
          Vector primal_solution_IPM(n);
          qIPM.primalSolution(primal_solution_IPM.valuesModifiable());
          Vector vector_d_dstar_diff_IPM(n);
          vector_d_dstar_diff_IPM.linearCombination(1, d_star, -1, primal_solution_IPM);

          // Increment time matrix
          time_array[d_type][n_index][m_index][0][test] = time_array[d_type][n_index][m_index][0][test] + time_DAS;
          time_array[d_type][n_index][m_index][1][test] = time_array[d_type][n_index][m_index][1][test] + time_IPM;

          // Print results
          printf("%+.4e  %+.4e  %+.4e  %6d  %+.4e  %+.4e  %+.4e  %6d\n",
                 qDAS.KKTError(),
                 vector_d_dstar_diff_DAS.normInf(),
                 time_DAS,
                 qDAS.numberOfIterations(),
                 qIPM.KKTError(),
                 vector_d_dstar_diff_IPM.normInf(),
                 time_IPM,
                 qIPM.numberOfIterations());

        } // end for (over test number)

        // Print new line
        printf("\n");

      } // end for (over m)

    } // end for (over n)

  } // end for (over d type)

  // Print times
  for (int i = 0; i < p1; i++) {
    for (int j = 0; j < p2; j++) {
      for (int k = 0; k < p3; k++) {
        for (int l = 0; l < p4; l++) {
          for (int z = 0; z < p5; z++) {
            printf(" %1d, %4d, %4d, %1d, %1d, %.4e\n", i, (j+1)*n_increment, (int)(m_factors[k] * (j+1)*n_increment) + m_additions[k], l, z, time_array[i][j][k][l][z]);
          }
        }
      }
    }
  }

  // Return
  return 0;

} // end testQPSolverImplementation

#endif /* __TESTQPSOLVER_HPP__ */
