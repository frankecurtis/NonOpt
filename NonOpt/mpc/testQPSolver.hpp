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

#include "NonOptQPSolverInteriorPoint.hpp"
#include "NonOptQPSolverDualActiveSet.hpp"
#include "NonOptSymmetricMatrix.hpp"
#include "NonOptSymmetricMatrixDense.hpp"

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

  // Set random seeds
  // generator.seed(42); // Instead of test that was 0...
  int test_number_total = 5;
  // for (int k = 0; k < test_number_total; k++){
  int k = 4;{
    std::cout << "Test Number: "<< k + 1 << std::endl;


    
    std::cout << std::setw(10) << "n" 
                                << std::setw(10) << "m" 
                                << std::setw(25) << "IPM KKT Error"
                                << std::setw(25) << "IPM d Error"
                                << std::setw(15) << "IPM Time"
                                << std::setw(15) << "IPM Iters"
                                << std::setw(20) << "DAS KKT Error"
                                << std::setw(25) << "DAS d Error"
                                << std::setw(15) << "DAS Time"
                                << std::setw(15) << "DAS Iters" 
                                << std::endl;
    std::cout << "--------------------------------------------------------------------------------------------" 
              << "--------------------------------------------------------------------------------------------"<< std::endl;

    for (int size_n : {500}) {
    // for (int size_n : {10, 100, 1000, 1200}) {
      // std::cout << "\n" << "\n"<< "\n" << "\n" << "Outer loop (n = " << size_n << "):\n";
      
      // Inner loop for m = n + 1 and m = 2n
      for (int size_m : { 2 * size_n }) {
      // for (int size_m : {size_n + 1, 3 * size_n / 2 , 2 * size_n }) {
        // std::cout << "\n" << "    m = " << size_m << std::endl;



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
        
        double scalar_ = 1;  // this is my delta!!
        Vector vector_ones_n(size_n, 1.0);
        Vector vector_ones_m(size_m, 1.0);



        Vector vector_d_star(size_n, scalar_); //for now, d=delta

        Vector vector_d_unc(size_n, 0);
        if (vector_d_star.normInf() == scalar_){    //not good because it only work for our case only
          vector_d_unc.addScaledVector(2*scalar_, vector_ones_n);
        }

        Vector vector_t(size_n, 0.0);
        vector_t.addScaledVector(-1,vector_d_unc);

        // std::cout << "Elements of RHS_matrix:" << std::endl;
        // for (int i = 0; i < size_n; i++) {
        //     // Access element LHS(i, total_length)
            
        //     std::cout << "dunc(" << i << ", " << ") = " << vector_d_unc.values()[i] << std::endl;
        // }

        // std::cout << "Elements of RHS_matrix:" << std::endl;
        // for (int i = 0; i < size_n; i++) {
        //     // Access element LHS(i, total_length)
            
        //     std::cout << "t(" << i << ", " << ") = " << vector_t.values()[i] << std::endl;
        // }
        // std::exit(0); 

        std::shared_ptr<SymmetricMatrixDense> matrix = std::make_shared<SymmetricMatrixDense>();
        matrix->setAsDiagonal(size_n, 1.0);  // W is the identity matrix

        std::vector<std::shared_ptr<Vector>> vector_list_G_hat;

        // Generate the outer vector with 'numberVectors' inner vectors
        for (int i = 0; i < size_n-1; i++) {
            // Create a new Vector with random components directly
            std::shared_ptr<Vector> new_vector = std::make_shared<Vector>(size_n);

            for (int j = 0; j < size_n; j++) {
                new_vector->set(j, normal(generator));  // Set random component directly
            }

            // Add it to the vector list
            vector_list_G_hat.push_back(new_vector);
        }
            // std::cout << "dunc("<< vector_list_G_hat[0]->values()[0]<< std::endl;




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

        // std::vector<std::shared_ptr<Vector>> vector_list_G_bar;
        

        // Create a new vector with the same size as vector_g_hat
        std::shared_ptr<Vector> g_hat_ptr = std::make_shared<Vector>(size_n);

        // Manually copy the data from vector_g_hat to the new vector
        for (int i = 0; i < size_n; ++i) {
            g_hat_ptr->set(i, vector_g_hat.values()[i]);  // Assuming 'set' and 'get' methods exist
        }

        // Add the shared pointer to the list
        vector_list_G_hat.push_back(g_hat_ptr);



        Vector omega_bar(size_n, 0.0);
        for (int i=0; i<size_n-1; i++){
          omega_bar.set(i, omega_hat.values()[i]);
        }
        omega_bar.set(size_n-1, 1);
        // omega_bar.scale(1.0 / omega_bar.norm1());

        // Scale all vectors in vector_list_G_hat by 2
        for (auto& vec_ptr : vector_list_G_hat) {
            vec_ptr->scale(omega_bar.norm1());  // Scale each vector by a factor of 2
        }

        
        // Extend vector_list_G_hat by adding m - n random vectors of size n
        for (int i = 0; i < size_m - size_n; i++) {
            // Step 1: Create a shared pointer to a new vector with the specified size
            std::shared_ptr<Vector> new_vector = std::make_shared<Vector>(size_n);

            // Step 2: Generate random components and populate the vector
            for (int j = 0; j < size_n; j++) {
                // Assuming `normal(generator)` generates a random number
                new_vector->set(j, normal(generator));  // Set the random value at position 'j'
            }

            // Step 3: Add the shared pointer to the vector list
            vector_list_G_hat.push_back(new_vector);
        }

        Vector vector_omega_star(size_m, 0.0);
        for (int i = 0; i < size_n; i++){
          vector_omega_star.set(i, omega_bar.values()[i] / omega_bar.norm1());
        }

        Vector lara(size_n, 0.0);

        for (int i = 0; i < size_n; i++) {
          lara.addScaledVector(vector_omega_star.values()[i], *vector_list_G_hat[i].get());
        }


        // for (int i = 0; i <size_n; i++) {
        //     std::cout << "Gomega*(" << i << ", " << ") = " << lara.values()[i] << std::endl;
        // }











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

        // for (int i = 0; i <size_n; i++) {
        //     std::cout << "rho(" << i << ", " << ") = " << vector_rho.values()[i] << std::endl;
        // }
        // // std::exit(0); 




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


        // std::vector<double> vector_b;

      

        // vector_b.linearCombination(my_scalar_u, vector_ones_m, -1, vector_v_omega);
        // vector_b.addScaledVector(1, vector_Gtd);

        // Step 1: Initialize vector_b as a std::vector<double> of size size_n
        std::vector<double> vector_bk(size_m, 0.0);  // Initialize with zeroes


        // Step 2: Perform linearCombination(my_scalar_u, vector_ones_m, -1, vector_v_omega)
        for (int i = 0; i < size_m; i++) {
            vector_bk[i] += my_scalar_u * 1 - vector_v_omega.values()[i];
        // }

        // // Step 3: Perform addScaledVector(1, vector_Gtd)
        // for (int i = 0; i < size_m; i++) {
            vector_bk[i] += vector_Gtd.values()[i];  // Add vector_Gtd directly (scaled by 1)
        }

        // std::cout << "Elements of RHS_matrix:" << std::endl;
        // for (int i = 0; i < size_n; i++) {
        //     // Access element LHS(i, total_length)
            
        //     std::cout << "v rho(" << i << ", " << ") = " << vector_v_rho.values()[i] << std::endl;
        // }

        // std::cout << "Elements of vector b_k:" << std::endl;
        // for (int i = 0; i < size_m; i++) {
        //     // Access element LHS(i, total_length)
            
        //     std::cout << "bk(" << i << ", " << ") = " << vector_bk[i] << std::endl;
        // }

        // std::exit(0); 

        // // Set QP data
        // q.setMatrix(matrix);
        // q.setVectorList(vector_list);
        // q.setVector(vector);
        // q.setScalar(regularization);
        // std::cout << " " << vector_list_G_hat[0]->values()[0]
        // << " " << vector_list_G_hat[0]->values()[1]
        // << " " << vector_list_G_hat[0]->values()[2]
        // << " " << vector_list_G_hat[1]->values()[0]
        // << " " << vector_list_G_hat[1]->values()[1]
        // << " " << vector_list_G_hat[1]->values()[2]
        // << " " << vector_list_G_hat[2]->values()[0]
        // << " " << vector_list_G_hat[2]->values()[1]
        // << " " << vector_list_G_hat[2]->values()[2]
        // << " " << vector_list_G_hat[3]->values()[0]
        // << " " << vector_list_G_hat[3]->values()[1]
        // << " " << vector_list_G_hat[3]->values()[2]<< std::endl;
        // std::exit(0); 


        q_IPM.setMatrix(matrix);
        q_IPM.setVectorList(vector_list_G_hat);
        q_IPM.setVector(vector_bk);
        q_IPM.setScalar(scalar_);

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
        // std::exit(0); 


        // Measure time for q_DAS.solveQP
        std::clock_t start_DAS = std::clock();
        q_DAS.solveQP(&options, &reporter, &quantities);
        std::clock_t end_DAS = std::clock();
        double time_DAS = static_cast<double>(end_DAS - start_DAS) / CLOCKS_PER_SEC;


          // std::cout << "Stopping program here for debugging purposes." << std::endl;
          // std::exit(0); 
        // Print all elements in the last column of LHS_matrix (column index = total_length)
        // std::cout << "Elements of RHS_matrix:" << std::endl;
        // for (int i = 0; i < size_m; i++) {
        //     // Access element LHS(i, total_length)
            
        //     std::cout << "sigma(" << i << ", " << ") = " << vector_bk[i] << std::endl;
        // }
        // std::cout << "G value(" << vector_list_G_hat[184]->values()[193] << std::endl;

        // for (int i = 0; i < size_n; i++) {
        //   std::cout << "d star*(" << i << ", " << ") = " << vector_d_star.values()[i] << std::endl;
        // }

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
        if (q_IPM.status() != QP_SUCCESS || q_IPM.KKTError() > 1e-03 || q_IPM.KKTErrorDual() > 1e-03 || primal_solution_IPM.normInf() > 1.0 / ((double)(0) + 1.0) + 1e-04) {
          result_IPM = 1;
        }


        // Get primal solution
        Vector primal_solution_DAS(size_n);
        q_DAS.primalSolution(primal_solution_DAS.valuesModifiable());   // This is -W(G*omega + gamma) right?

        Vector vector_d_dstar_diff_DAS(size_n);
        vector_d_dstar_diff_DAS.linearCombination(1, vector_d_star, -1, primal_solution_DAS);  // Difference between d and d_star
      

        // Check results CHECK LAST COMPONENT
        if (q_DAS.status() != QP_SUCCESS || q_DAS.KKTError() > 1e-03 || q_DAS.KKTErrorDual() > 1e-03 || primal_solution_DAS.normInf() > 1.0 / ((double)(0) + 1.0) + 1e-04) {
          result_DAS = 1;
        }

        // reporter.printf(R_QP, R_BASIC, 
        //     "IPM:  KKT error = %+.4e  d error = %+.4e  Total time = %.7f  iters = %6d\n", 
        //     q_IPM.KKTError(), vector_d_dstar_diff_IPM.normInf(),time_IPM, q_IPM.numberOfIterations());

        // reporter.printf(R_QP, R_BASIC, 
        //     "DAS:  KKT error = %+.16e  d error = %+.4e  Total time = %.7f  iters = %6d\n", 
        //     q_DAS.KKTError(), vector_d_dstar_diff_DAS.normInf(),time_DAS, q_DAS.numberOfIterations());



        // std::cout << std::setw(10) << "n" 
        //                           << std::setw(10) << "m" 
        //                           << std::setw(25) << "IPM KKT Error"
        //                           << std::setw(25) << "IPM d Error"
        //                           << std::setw(15) << "IPM Time"
        //                           << std::setw(10) << "IPM Iterations"
        //                           << std::setw(25) << "DAS KKT Error"
        //                           << std::setw(25) << "DAS d Error"
        //                           << std::setw(15) << "DAS Time"
        //                           << std::setw(10) << "DAS Iterations" 
        //                           << std::endl;
        // std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
        
        std::cout << std::setw(10) << size_n 
                          << std::setw(10) << size_m
                          << std::setw(25) << std::scientific << q_IPM.KKTError()
                          << std::setw(25) << vector_d_dstar_diff_IPM.normInf()
                          << std::setw(15) << std::fixed << std::setprecision(4) << time_IPM
                          << std::setw(10) << q_IPM.numberOfIterations()
                          << std::setw(25) << std::scientific << q_DAS.KKTError()
                          << std::setw(25) << vector_d_dstar_diff_DAS.normInf()
                          << std::setw(15) << std::fixed << std::setprecision(4) << time_DAS
                          << std::setw(10) << q_DAS.numberOfIterations()
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
    // Check option
    if (option == 1) {
        // Print final messages
        if (result_IPM == 0) {
            reporter.printf(R_QP, R_BASIC, "IPM TEST WAS SUCCESSFUL.\n");
        } else {
            reporter.printf(R_QP, R_BASIC, "IPM TEST FAILED.\n");
        }

        if (result_DAS == 0) {
            reporter.printf(R_QP, R_BASIC, "DAS TEST WAS SUCCESSFUL.\n");
        } else {
            reporter.printf(R_QP, R_BASIC, "DAS TEST FAILED.\n");
        }
    } // end if
    std::cout << " " << std::endl;
}
  

  
  std::clock_t end = std::clock();

  // Calculate the elapsed time in seconds
  double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;

  // Print the result
  std::cout << "Execution total time: " << elapsed_time << " seconds" << std::endl;

  // Return
  return result_IPM;

} // end testQPSolverImplementation

#endif /* __TESTQPSOLVER_HPP__ */
