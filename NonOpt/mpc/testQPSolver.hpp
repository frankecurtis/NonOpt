// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Baoyu Zhou

#ifndef __TESTQPSOLVER_HPP__
#define __TESTQPSOLVER_HPP__

#include <iostream>

#include <cmath>
#include <random>
#include <iomanip>
#include <iostream>

#include "NonOptQPSolverInteriorPoint.hpp"
#include "NonOptQPSolverDualActiveSet.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptSymmetricMatrixDense.hpp"

#include "NonOptDeclarations.hpp"  // Lara
#include "NonOptDefinitions.hpp"   // Lara

using namespace NonOpt;

// Implementation of test
int testQPSolverImplementation(int option)
{
  std::clock_t start = std::clock();

  // Initialize output
  int result_DAS = 0;
  int result_IPM = 0;

  // Declare reporter
  Reporter reporter;
  Quantities quantities;

  // Check option
  if (option == 1) {

    // Declare stream report
    // std::shared_ptr<StreamReport> s(new StreamReport("s", R_QP, R_PER_INNER_ITERATION));
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
  
  std::cout << "n,"
          << "m,"
          << "IPM_KKT_Error,"
          << "IPM_d_Error,"
          << "IPM_Time,"
          << "IPM_Iters,"
          << "DAS_KKT_Error,"
          << "DAS_d_Error,"
          << "DAS_Time,"
          << "DAS_Iters"
          << std::endl;


  int test_number_total = 10;
  for (int k = 0; k < test_number_total; k++){
  // int k = 5;{
    // std::cout << "Test Number: "<< k + 1 << std::endl;

  // std::array<double, 4> mu_aff_factors = {1, 1.5, 2.0, 2.5};
  std::array<double, 1> mu_aff_factors = {1};
  // std::array<double, 3> init_params = {1e-3, 1e-2, 1e-1};
  std::array<double, 1> init_params = {1e-1};

  for (double mu_aff_factor : mu_aff_factors) {
    for (double init_param : init_params) {

      // std::cout << "mu_aff_factor = " << mu_aff_factor << std::endl;
      // std::cout << "init_param = " << init_param << std::endl;

    
      // std::cout << std::setw(10) << "n" 
      //                             << std::setw(10) << "m" 
      //                             << std::setw(25) << "IPM_KKT_Error"
      //                             << std::setw(25) << "IPM_d_Error"
      //                             << std::setw(15) << "IPM_Time"
      //                             << std::setw(15) << "IPM_Iters"
      //                             << std::setw(20) << "DAS_KKT_Error"
      //                             << std::setw(25) << "DAS_d_Error"
      //                             << std::setw(15) << "DAS_Time"
      //                             << std::setw(15) << "DAS_Iters" 
      //                             << std::endl;
      // std::cout << "--------------------------------------------------------------------------------------------" 
      //           << "--------------------------------------------------------------------------------------------"<< std::endl;

      for (int size_n : {300}) {
      // for (int size_n = 100; size_n <= 2500; size_n += 200) {
        for (int size_m : { size_n + 1}) {
        // for (int size_m : {size_n + 1, 3 * size_n / 2 , 2 * size_n }) {



          // Declare options
          Options options;

          // Declare QP solver object
          QPSolverInteriorPoint q_IPM; //add here
          QPSolverDualActiveSet q_DAS; //add here

          // Add options
          q_IPM.addOptions(&options);
          q_DAS.addOptions(&options);

          // Use exact solves
          options.modifyIntegerValue("QPDAS_iteration_limit_minimum", 1e+6);  // Check this!!
          options.modifyDoubleValue("QPDAS_kkt_tolerance", 1e-5);  // Check this!!
          options.modifyBoolValue("QPDAS_allow_inexact_termination", false);  // Check this!!

          // Set options
          q_IPM.setOptions(&options);
          q_DAS.setOptions(&options);

          // Initialize data
          q_IPM.initializeData(size_n);
          q_DAS.initializeData(size_n);


            // Print test number
            // reporter.printf(R_QP, R_BASIC, "Running test %8d... ", test);

            // Set random seeds
            generator.seed(k); // Instead of test that was 0...

        
            /**
             * Lara 
             */
            
            // double scalar_ = NONOPT_DOUBLE_INFINITY;  // this is my delta!!
            double scalar_ = 10;  // this is my delta!!
            Vector vector_ones_n(size_n, 1.0);
            Vector vector_ones_m(size_m, 1.0);

            Vector vector_d_star(size_n, 0.0); 

            if (scalar_ == NONOPT_DOUBLE_INFINITY){
              vector_d_star.addScaledVector(0.0, vector_ones_n);
            }
            else{
              // vector_d_star.addScaledVector(scalar_, vector_ones_n); //for now, d=delta (it can also be equal to 0 in this code)
              vector_d_star.addScaledVector(0.0, vector_ones_n); 
            }

            Vector vector_d_unc(size_n, 0);
            if (scalar_ == NONOPT_DOUBLE_INFINITY){
              vector_d_unc.addScaledVector(1, vector_d_star);
            }
            else{
            if (vector_d_star.normInf() == scalar_){
              vector_d_unc.addScaledVector(2*scalar_, vector_ones_n);
            }
            }

            
            Vector vector_t(size_n, 0.0);
            vector_t.addScaledVector(-1,vector_d_unc);


            std::shared_ptr<SymmetricMatrixDense> matrix = std::make_shared<SymmetricMatrixDense>();
            matrix->setAsDiagonal(size_n, 1.0);  // W is the identity matrix

            std::vector<std::shared_ptr<Vector>> vector_list_G_hat;

            // Generate the outer vector with 'numberVectors' inner vectors
            for (int i = 0; i < size_n-1; i++) {
                // Create a new Vector with random components
                std::shared_ptr<Vector> new_vector = std::make_shared<Vector>(size_n);

                for (int j = 0; j < size_n; j++) {
                    new_vector->set(j, normal(generator));  // Set random component
                }

                // Add it to the vector list
                vector_list_G_hat.push_back(new_vector);
            }


            // // Initialize weights
            Vector omega_hat(size_n-1, 0.0);

            // Set weights and normalize
            for (int i = 0; i < size_n-1; i++) {
              omega_hat.set(i, uniform(generator));
            }


            Vector vector_g_hat(size_n, 0.0);
            
            Vector vector_Ghat_omegaHat(size_n, 0.0);
            
            for (int i = 0; i < size_n-1; i++) {
              vector_Ghat_omegaHat.addScaledVector(omega_hat.values()[i], *vector_list_G_hat[i].get());
            }

            vector_g_hat.linearCombination(1, vector_t, -1, vector_Ghat_omegaHat);

            // Create a new vector with the same size as vector_g_hat
            std::shared_ptr<Vector> g_hat_ptr = std::make_shared<Vector>(size_n);

            // Manually copy the data from vector_g_hat to the new vector
            for (int i = 0; i < size_n; ++i) {
                g_hat_ptr->set(i, vector_g_hat.values()[i]);
            }

            // Add the shared pointer to the list
            vector_list_G_hat.push_back(g_hat_ptr);



            Vector omega_bar(size_n, 0.0);
            for (int i=0; i<size_n-1; i++){
              omega_bar.set(i, omega_hat.values()[i]);
            }
            omega_bar.set(size_n-1, 1);
            // omega_bar.scale(1.0 / omega_bar.norm1());

            // Scale all vectors in vector_list_G_hat
            for (auto& vec_ptr : vector_list_G_hat) {
                vec_ptr->scale(omega_bar.norm1());
            }

            
            // Extend vector_list_G_hat by adding m - n random vectors of size n
            for (int i = 0; i < size_m - size_n; i++) {
                std::shared_ptr<Vector> new_vector = std::make_shared<Vector>(size_n);

                for (int j = 0; j < size_n; j++) {
                    new_vector->set(j, normal(generator)); 
                }

                vector_list_G_hat.push_back(new_vector);
            }

            Vector vector_omega_star(size_m, 0.0);
            for (int i = 0; i < size_n; i++){
              vector_omega_star.set(i, omega_bar.values()[i] / omega_bar.norm1());
            }


            Vector vector_q(size_n, 0.0);

            Vector vector_G_omegaStar(size_n, 0.0);
            
            for (int i = 0; i < size_m ; i++) {
              vector_G_omegaStar.addScaledVector(vector_omega_star.values()[i], *vector_list_G_hat[i].get());
            }

            vector_q.linearCombination(-1, vector_d_star, -1, vector_G_omegaStar);

            Vector vector_sigma(size_n, 0.0);
            Vector vector_rho(size_n, 0.0);

            for (int i = 0; i < size_n; i++){
              if (vector_q.values()[i] > 0){
                vector_sigma.set(i, vector_q.values()[i]);
              }
              else{
                vector_rho.set(i, -vector_q.values()[i]);
              }
            }



            double my_scalar_u = 5;


            Vector vector_d(size_n, 0.0);  // Gomega + sigma - rho
            
            vector_d.linearCombination(1, vector_G_omegaStar, 1, vector_sigma);
            vector_d.addScaledVector(-1, vector_rho);


            Vector vector_v_sigma(size_n, 0.0);
            Vector vector_v_rho(size_n, 0.0);


            vector_v_sigma.linearCombination(1, vector_d, scalar_, vector_ones_n);
            vector_v_rho.linearCombination(-1, vector_d, scalar_, vector_ones_n);



            Vector vector_v_omega(size_m, 0.0); 

            for (int i = 0; i < size_m; i++){
              if (vector_omega_star.values()[i] == 0){
                vector_v_omega.set(i, abs(uniform(generator)) );
              }
            }

            Vector vector_Gtd(size_m, 0.0);
            for (int i = 0; i < size_m; i++) {
              vector_Gtd.values()[i] = vector_list_G_hat[i]->innerProduct(vector_d);
            }

            std::vector<double> vector_bk(size_m, 0.0);  // Initialize with zeroes


            for (int i = 0; i < size_m; i++) {
                vector_bk[i] += my_scalar_u * 1 - vector_v_omega.values()[i];
                vector_bk[i] += vector_Gtd.values()[i];  // Add vector_Gtd directly (scaled by 1)
            }

            if (scalar_ == NONOPT_DOUBLE_INFINITY){
              for (int i = 0; i < size_m; i++) {
                vector_bk[i] *= 0;
                vector_bk[i] += my_scalar_u * 1 - vector_v_omega.values()[i];
                vector_bk[i] += vector_list_G_hat[i]->innerProduct(vector_G_omegaStar);
                }
            }

            q_IPM.setMatrix(matrix);
            q_IPM.setVectorList(vector_list_G_hat);
            q_IPM.setVector(vector_bk);
            q_IPM.setScalar(scalar_);
            q_IPM.setMuAffFactor(mu_aff_factor); //lara
            q_IPM.setInitParam(init_param); //lara

            q_DAS.setMatrix(matrix);
            q_DAS.setVectorList(vector_list_G_hat);
            q_DAS.setVector(vector_bk);
            q_DAS.setScalar(scalar_);


            // Solve QP
            // Measure time for q_IPM.solveQP
            std::clock_t start_IPM = std::clock();
            q_IPM.solveQP(&options, &reporter, &quantities);
            std::clock_t end_IPM = std::clock();
            double time_IPM = static_cast<double>(end_IPM - start_IPM) / CLOCKS_PER_SEC;


            // Measure time for q_DAS.solveQP
            std::clock_t start_DAS = std::clock();
            q_DAS.solveQP(&options, &reporter, &quantities);
            std::clock_t end_DAS = std::clock();
            double time_DAS = static_cast<double>(end_DAS - start_DAS) / CLOCKS_PER_SEC;


            // Check for pass or fail
            if (q_IPM.status() == QP_SUCCESS) {
              reporter.printf(R_QP, R_PER_INNER_ITERATION, "pass");  //Lara changed to R_PER...
            }
            else {
              reporter.printf(R_QP, R_PER_INNER_ITERATION, "fail");
            }

            // Check for pass or fail
            if (q_DAS.status() == QP_SUCCESS) {
              reporter.printf(R_QP, R_PER_INNER_ITERATION, "pass");
            }
            else {
              reporter.printf(R_QP, R_PER_INNER_ITERATION, "fail");
            }

            // Get primal solution
            Vector primal_solution_IPM(size_n);
            q_IPM.primalSolution(primal_solution_IPM.valuesModifiable());   // This is -W(G*omega + gamma) right?

            Vector vector_d_dstar_diff_IPM(size_n);
            vector_d_dstar_diff_IPM.linearCombination(1, vector_d_star, -1, primal_solution_IPM);  // Difference between d and d_star
          

            // Check results CHECK LAST COMPONENT
            if (q_IPM.status() != QP_SUCCESS || q_IPM.KKTError() > 1e-03 || q_IPM.KKTErrorDual() > 1e-03 || vector_d_dstar_diff_IPM.normInf() > 1e-03) {
              result_IPM = 1;
            }


            // Get primal solution
            Vector primal_solution_DAS(size_n);
            q_DAS.primalSolution(primal_solution_DAS.valuesModifiable());

            Vector vector_d_dstar_diff_DAS(size_n);
            vector_d_dstar_diff_DAS.linearCombination(1, vector_d_star, -1, primal_solution_DAS);  // Difference between d and d_star
          

            // Check results
            if (q_DAS.status() != QP_SUCCESS || q_DAS.KKTError() > 1e-03 || q_DAS.KKTErrorDual() > 1e-03 || primal_solution_DAS.normInf() > 1.0 / ((double)(0) + 1.0) + 1e-04) {
              result_DAS = 1;
            }
            
            // std::cout << std::setw(10) << size_n 
            //                   << std::setw(10) << size_m
            //                   << std::setw(25) << std::scientific << q_IPM.KKTError()
            //                   << std::setw(25) << vector_d_dstar_diff_IPM.normInf()
            //                   << std::setw(15) << std::fixed << std::setprecision(4) << time_IPM
            //                   << std::setw(10) << q_IPM.numberOfIterations()
            //                   << std::setw(25) << std::scientific << q_DAS.KKTError()
            //                   << std::setw(25) << vector_d_dstar_diff_DAS.normInf()
            //                   << std::setw(15) << std::fixed << std::setprecision(4) << time_DAS
            //                   << std::setw(10) << q_DAS.numberOfIterations()
            //                   << std::endl;

            std::cout << size_n << ","
              << size_m << ","
              << std::scientific << q_IPM.KKTError() << ","
              << std::scientific << vector_d_dstar_diff_IPM.normInf() << ","
              << std::fixed << std::setprecision(4) << time_IPM << ","
              << q_IPM.numberOfIterations() << ","
              << std::scientific << q_DAS.KKTError() << ","
              << std::scientific << vector_d_dstar_diff_DAS.normInf() << ","
              << std::fixed << std::setprecision(4) << time_DAS << ","
              << q_DAS.numberOfIterations()
              << std::endl;

            // // Print solve status information
            // reporter.printf(R_QP, R_BASIC, "  status: %d  iters: %6d  kkt error: %+.4e  kkt error (dual): %+.4e  ||step||_inf: %+.4e\n", q.status(), q.numberOfIterations(), q.KKTError(), q.KKTErrorDual(), primal_solution.normInf());
            // reporter.printf(R_QP, R_BASIC, "  status: %d  iters: %6d  kkt error: %+.4e  kkt error (dual): %+.4e  ||step||_inf: %+.4e\n", q.status(), q.numberOfIterations(), q.KKTError(), q.KKTErrorDual(), primal_solution.normInf());
            // std::cout << "Hi Lara Hi 2"<< std::endl;
            // // std::cout << "time = " << q.elapsedCPUtime() << std::endl;
            // reporter.printf(R_QP, R_BASIC, "time = %.14f\n", q.elapsedCPUtime());
            // reporter.printf(R_QP, R_BASIC, "kkt_error = %.14f\n", q.KKTError());
            // reporter.printf(R_QP, R_BASIC, "d_error = %.14f\n", vector_d_dstar_diff.normInf() );

            //   // Print solve status information
            // reporter.printf(R_QP, R_PER_INNER_ITERATION, "  status: %d  iters: %6d  kkt error: %+.4e  kkt error (dual): %+.4e  ||step||_inf: %+.4e\n", q.status(), q.numberOfIterations(), q.KKTError(), q.KKTErrorDual(), primal_solution.normInf());
            // // std::cout << "Hi Lara Hi 2"<< std::endl;
            // // std::cout << "time = " << q.elapsedCPUtime() << std::endl;
            // reporter.printf(R_QP, R_PER_INNER_ITERATION, "time = %.14f\n", q.elapsedCPUtime());
            // reporter.printf(R_QP, R_PER_INNER_ITERATION, "kkt_error = %.14f\n", q.KKTError());
            // reporter.printf(R_QP, R_PER_INNER_ITERATION, "d_error = %.14f\n", vector_d_dstar_diff.normInf() );

          } // end for

        }
        // // Check option
        // if (option == 1) {
        //     // Print final messages
        //     if (result_IPM == 0) {
        //         reporter.printf(R_QP, R_BASIC, "IPM TEST WAS SUCCESSFUL.\n");
        //     } else {
        //         reporter.printf(R_QP, R_BASIC, "IPM TEST FAILED.\n");
        //     }

        //     if (result_DAS == 0) {
        //         reporter.printf(R_QP, R_BASIC, "DAS TEST WAS SUCCESSFUL.\n");
        //     } else {
        //         reporter.printf(R_QP, R_BASIC, "DAS TEST FAILED.\n");
        //     }
        // } // end if
        // std::cout << " " << std::endl;
      }
      }
      }
  

  
  std::clock_t end = std::clock();

  // Calculate the elapsed time in seconds
  double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;

  // std::cout << "Execution total time: " << elapsed_time << " seconds" << std::endl;

  // Return
  return result_IPM;

} // end testQPSolverImplementation

#endif /* __TESTQPSOLVER_HPP__ */
