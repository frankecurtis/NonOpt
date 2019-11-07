// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTNONOPT_HPP__
#define __TESTNONOPT_HPP__

#include <cmath>
#include <cstring>

#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"

#include "ActiveFaces.hpp"
#include "BrownFunction_2.hpp"
#include "ChainedCB3_1.hpp"
#include "ChainedCB3_2.hpp"
#include "ChainedCrescent_1.hpp"
#include "ChainedCrescent_2.hpp"
#include "ChainedLQ.hpp"
#include "ChainedMifflin_2.hpp"
#include "MaxQ.hpp"
#include "MxHilb.hpp"
#include "QuadPoly.hpp"

#include "Test29_2.hpp"
#include "Test29_2.hpp"
#include "Test29_5.hpp"
#include "Test29_5.hpp"
#include "Test29_6.hpp"
#include "Test29_6.hpp"
#include "Test29_11.hpp"
#include "Test29_11.hpp"
#include "Test29_13.hpp"
#include "Test29_13.hpp"
#include "Test29_17.hpp"
#include "Test29_17.hpp"
#include "Test29_19.hpp"
#include "Test29_19.hpp"
#include "Test29_20.hpp"
#include "Test29_20.hpp"
#include "Test29_22.hpp"
#include "Test29_22.hpp"
#include "Test29_24.hpp"
#include "Test29_24.hpp"

using namespace NonOpt;

// Implementation of test
int runExperiments()
{

  // Declare NonOptSolver
  NonOptSolver nonopt;

  // Delete reports
  nonopt.reporter()->deleteReports();

  // Declare problem pointer
  std::shared_ptr<Problem> problem;

//  // Declare vectors of strings for strategies
//  std::vector<std::string> direction_computation_names;
//  std::vector<std::string> inverse_hessian_update_names;
//  std::vector<std::string> line_search_names;
//  std::vector<std::string> point_set_update_names;
//  std::vector<std::string> qp_solver_names;
//  std::vector<std::string> symmetric_matrix_names;
//
//  // Push names of strategies
//  direction_computation_names.push_back("CuttingPlane");
//  direction_computation_names.push_back("GradientCombination");
//  direction_computation_names.push_back("Gradient");
//  inverse_hessian_update_names.push_back("BFGS");
//  line_search_names.push_back("WeakWolfe");
//  line_search_names.push_back("Backtracking");
//  point_set_update_names.push_back("Proximity");
//  qp_solver_names.push_back("ActiveSet");
//  symmetric_matrix_names.push_back("Dense");
//  symmetric_matrix_names.push_back("LimitedMemory");



    std::vector<std::string> settings_names;
    //settings_names.push_back("original");
	settings_names.push_back("only_inexact");

	for (int settings_number = 0; settings_number < settings_names.size(); settings_number++) {

	            // Set qp solver strategy
	     nonopt.options()->modifyOptionsFromFile(nonopt.reporter(), settings_names[settings_number]+".txt");

	     std::string settings_name=settings_names[settings_number];



//              // Declare strategy names
//              std::string direction_computation_name;
//              std::string inverse_hessian_update_name;
//              std::string line_search_name;
//              std::string point_set_update_name;
//              std::string qp_solver_name;
//              std::string symmetric_matrix_name;
//
//              // Read integer options
//              nonopt.options()->valueAsString(nonopt.reporter(), "direction_computation", direction_computation_name);
//              nonopt.options()->valueAsString(nonopt.reporter(), "inverse_hessian_update", inverse_hessian_update_name);
//              nonopt.options()->valueAsString(nonopt.reporter(), "line_search", line_search_name);
//              nonopt.options()->valueAsString(nonopt.reporter(), "point_set_update", point_set_update_name);
//              nonopt.options()->valueAsString(nonopt.reporter(), "qp_solver", qp_solver_name);
//              nonopt.options()->valueAsString(nonopt.reporter(), "symmetric_matrix", symmetric_matrix_name);
//
//              // Print header
//              printf("%s, %s, %s, %s, %s, %s\n",
//                     direction_computation_name.c_str(),
//                     inverse_hessian_update_name.c_str(),
//                     line_search_name.c_str(),
//                     point_set_update_name.c_str(),
//                     qp_solver_name.c_str(),
//                     symmetric_matrix_name.c_str());
	     	 printf("Start testing %s\n",settings_name.c_str());
              printf("================================================================================================================================================================\n");
              printf("Problem            Status   Stat. Rad.    Optimal     Objective    Iterations    Inn. Iters.    QP Iters.    Func. Eval.   Grad. Eval.      Time     NonOpt Time\n");
              printf("================================================================================================================================================================\n");

              // Declare optimal value
              double optimal_value;

              // Loop through test problems
              for (int problem_count = 0; problem_count < 11; problem_count++) {

                // Delete reports
                nonopt.reporter()->deleteReports();

                // Declare file report
                std::shared_ptr<FileReport> r(new FileReport("f", R_NL, R_PER_INNER_ITERATION));
                std::shared_ptr<FileReport> r_QP(new FileReport("f", R_QP, R_PER_INNER_ITERATION_IN));

                // Declare output file name
                char out_file[100];

                // Switch on problems
                switch (problem_count) {
                  case 0:
                    problem = std::make_shared<ActiveFaces>();
                    printf("ActiveFaces        ");
                    strcpy(out_file, (char*)"output/ActiveFaces");
                    optimal_value = 0.0;
                    break;
                  case 1:
                    problem = std::make_shared<BrownFunction_2>();
                    printf("BrownFunction_2    ");
                    strcpy(out_file, (char*)"output/BrownFunction_2");
                    optimal_value = 0.0;
                    break;
                  case 2:
                    problem = std::make_shared<ChainedCB3_1>();
                    printf("ChainedCB3_1       ");
                    strcpy(out_file, (char*)"output/ChainedCB3_1");
                    optimal_value = 49.0 * 2.0;
                    break;
                  case 3:
                    problem = std::make_shared<ChainedCB3_2>();
                    printf("ChainedCB3_2       ");
                    strcpy(out_file, (char*)"output/ChainedCB3_2");
                    optimal_value = 49.0 * 2.0;
                    break;
                  case 4:
                    problem = std::make_shared<ChainedCrescent_1>();
                    printf("ChainedCrescent_1  ");
                    strcpy(out_file, (char*)"output/ChainedCrescent_1");
                    optimal_value = 0.0;
                    break;
                  case 5:
                    problem = std::make_shared<ChainedCrescent_2>();
                    printf("ChainedCrescent_2  ");
                    strcpy(out_file, (char*)"output/ChainedCrescent_2");
                    optimal_value = 0.0;
                    break;
                  case 6:
                    problem = std::make_shared<ChainedLQ>();
                    printf("ChainedLQ          ");
                    strcpy(out_file, (char*)"output/ChainedLQ");
                    optimal_value = -49.0 * sqrt(2.0);
                    break;
                  case 7:
                    problem = std::make_shared<ChainedMifflin_2>();
                    printf("ChainedMifflin_2   ");
                    strcpy(out_file, (char*)"output/ChainedMifflin_2");
                    optimal_value = -34.795;
                    break;
                  case 8:
                    problem = std::make_shared<MaxQ>();
                    printf("MaxQ               ");
                    strcpy(out_file, (char*)"output/MaxQ");
                    optimal_value = 0.0;
                    break;
                  case 9:
                    problem = std::make_shared<MxHilb>();
                    printf("MxHilb             ");
                    strcpy(out_file, (char*)"output/MxHilb");
                    optimal_value = 0.0;
                    break;
                  case 10:
                    problem = std::make_shared<QuadPoly>(50, 10, 5, 10.0);
                    printf("QuadPoly           ");
                    strcpy(out_file, (char*)"output/QuadPoly");
                    optimal_value = 0.0;
                    break;
                  case 11:
                    problem = std::make_shared<Test29_2>();
                    printf("Test29_2           ");
                    strcpy(out_file, (char*)"output/Test29_2");
                    optimal_value = 0.0;
                    break;
                  case 12:
                    problem = std::make_shared<Test29_5>();
                    printf("Test29_5           ");
                    strcpy(out_file, (char*)"output/Test29_5");
                    optimal_value = 0.0;
                    break;
                  case 13:
                    problem = std::make_shared<Test29_6>();
                    printf("Test29_6           ");
                    strcpy(out_file, (char*)"output/Test29_6");
                    optimal_value = 0.0;
                    break;
                  case 14:
                    problem = std::make_shared<Test29_11>();
                    printf("Test29_11           ");
                    strcpy(out_file, (char*)"output/Test29_11");
                    optimal_value = 20.0;
                    break;
                  case 15:
                    problem = std::make_shared<Test29_13>();
                    printf("Test29_13           ");
                    strcpy(out_file, (char*)"output/Test29_13");
                    optimal_value = 0.94;
                    break;
                  case 16:
                    problem = std::make_shared<Test29_17>();
                    printf("Test29_17           ");
                    strcpy(out_file, (char*)"output/Test29_17");
                    optimal_value = 0.0;
                    break;
                  case 17:
                    problem = std::make_shared<Test29_19>();
                    printf("Test29_19           ");
                    strcpy(out_file, (char*)"output/Test29_19");
                    optimal_value = 0.0;
                    break;
                  case 18:
                    problem = std::make_shared<Test29_20>();
                    printf("Test29_20           ");
                    strcpy(out_file, (char*)"output/Test29_20");
                    optimal_value = 0.0;
                    break;
                  case 19:
                    problem = std::make_shared<Test29_22>();
                    printf("Test29_22           ");
                    strcpy(out_file, (char*)"output/Test29_22");
                    optimal_value = 0.0;
                    break;
                  case 20:
                    problem = std::make_shared<Test29_24>();
                    printf("Test29_24           ");
                    strcpy(out_file, (char*)"output/Test29_24");
                    optimal_value = 0.4;
                    break;
                }  // end switch

                // Set rest of output file name
                //sprintf(out_file, "%s_%d_%d_%d_%d_%d_%d.out", out_file, direction_computation_number, inverse_hessian_update_number, line_search_number, point_set_update_number, qp_solver_number, symmetric_matrix_number);
                sprintf(out_file, "%s_%s.out", out_file, settings_name.c_str());
                char out_file_QP[100];
                std::copy(out_file, out_file+100, out_file_QP);
                sprintf(out_file_QP, "%s_%s_qp.out", out_file_QP, settings_name.c_str());
                // Open output file
                r->open(out_file);
                r_QP->open(out_file_QP);

                // Add to reporter
                nonopt.reporter()->addReport(r);
                nonopt.reporter()->addReport(r_QP);

                // Optimize
                nonopt.optimize(problem);

                // Print results
                printf("%6d  %+.4e  %+.4e  %+.4e  %12d  %12d %12d %12d  %12d  %+.4e  %+.4e\n",
                       nonopt.status(),
                       nonopt.stationarityRadius(),
                       optimal_value,
                       nonopt.objective(),
                       nonopt.iterations(),
                       nonopt.totalInnerIterations(),
					   nonopt.totalQPIterations(),
                       nonopt.functionEvaluations(),
                       nonopt.gradientEvaluations(),
                       nonopt.time(),
                       nonopt.timeNonOpt());

              }  // end for

              // Print new line
              printf("\n");

	}

  // Return
  return 0;

}  // end testNonOptImplementation

#endif /* __TESTNONOPT_HPP__ */
