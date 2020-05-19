// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <string>

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
#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"
#include "QuadPoly.hpp"
#include "Test29_11.hpp"
#include "Test29_13.hpp"
#include "Test29_17.hpp"
#include "Test29_19.hpp"
#include "Test29_2.hpp"
#include "Test29_20.hpp"
#include "Test29_22.hpp"
#include "Test29_24.hpp"
#include "Test29_5.hpp"
#include "Test29_6.hpp"

// Main function
int main()
{

  // Declare NonOptSolver
  NonOptSolver nonopt;

  // Modify options from file
  nonopt.options()->modifyOptionsFromFile(nonopt.reporter(), "nonopt.opt");

  // Delete reports
  nonopt.reporter()->deleteReports();

  // Declare problem pointer
  std::shared_ptr<Problem> problem;

  // Declare vectors of strings for strategies
  std::vector<std::string> direction_computation_names;
  std::vector<std::string> inverse_hessian_update_names;
  std::vector<std::string> line_search_names;
  std::vector<std::string> point_set_update_names;
  std::vector<std::string> qp_solver_names;
  std::vector<std::string> symmetric_matrix_names;

  // Push names of strategies
  direction_computation_names.push_back("CuttingPlane");
  direction_computation_names.push_back("GradientCombination");
  direction_computation_names.push_back("Gradient");
  inverse_hessian_update_names.push_back("BFGS");
  line_search_names.push_back("WeakWolfe");
  line_search_names.push_back("Backtracking");
  point_set_update_names.push_back("Proximity");
  qp_solver_names.push_back("DualActiveSet");
  symmetric_matrix_names.push_back("Dense");
  symmetric_matrix_names.push_back("LimitedMemory");

  // Loop through direction computations
  for (int direction_computation_number = 0; direction_computation_number < direction_computation_names.size(); direction_computation_number++) {

    // Set direction computation strategy
    nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", direction_computation_names[direction_computation_number]);

    // Loop through inverse Hessian update strategies
    for (int inverse_hessian_update_number = 0; inverse_hessian_update_number < inverse_hessian_update_names.size(); inverse_hessian_update_number++) {

      // Set inverse Hessian update strategy
      nonopt.options()->modifyStringValue(nonopt.reporter(), "inverse_hessian_update", inverse_hessian_update_names[inverse_hessian_update_number]);

      // Loop through line searches
      for (int line_search_number = 0; line_search_number < line_search_names.size(); line_search_number++) {

        // Set line search strategy
        nonopt.options()->modifyStringValue(nonopt.reporter(), "line_search", line_search_names[line_search_number]);

        // Loop through point set updates
        for (int point_set_update_number = 0; point_set_update_number < point_set_update_names.size(); point_set_update_number++) {

          // Set point set update strategy
          nonopt.options()->modifyStringValue(nonopt.reporter(), "point_set_update", point_set_update_names[point_set_update_number]);

          // Loop through qp solvers
          for (int qp_solver_number = 0; qp_solver_number < qp_solver_names.size(); qp_solver_number++) {

            // Set qp solver strategy
            nonopt.options()->modifyStringValue(nonopt.reporter(), "qp_solver", qp_solver_names[qp_solver_number]);

            // Loop through symmetric matrices
            for (int symmetric_matrix_number = 0; symmetric_matrix_number < symmetric_matrix_names.size(); symmetric_matrix_number++) {

              // Set symmetric matrix strategy
              nonopt.options()->modifyStringValue(nonopt.reporter(), "symmetric_matrix", symmetric_matrix_names[symmetric_matrix_number]);

              // Declare strategy names
              std::string direction_computation_name;
              std::string inverse_hessian_update_name;
              std::string line_search_name;
              std::string point_set_update_name;
              std::string qp_solver_name;
              std::string symmetric_matrix_name;

              // Read integer options
              nonopt.options()->valueAsString(nonopt.reporter(), "direction_computation", direction_computation_name);
              nonopt.options()->valueAsString(nonopt.reporter(), "inverse_hessian_update", inverse_hessian_update_name);
              nonopt.options()->valueAsString(nonopt.reporter(), "line_search", line_search_name);
              nonopt.options()->valueAsString(nonopt.reporter(), "point_set_update", point_set_update_name);
              nonopt.options()->valueAsString(nonopt.reporter(), "qp_solver", qp_solver_name);
              nonopt.options()->valueAsString(nonopt.reporter(), "symmetric_matrix", symmetric_matrix_name);

              // Print header
              printf("%s, %s, %s, %s, %s, %s\n",
                     direction_computation_name.c_str(),
                     inverse_hessian_update_name.c_str(),
                     line_search_name.c_str(),
                     point_set_update_name.c_str(),
                     qp_solver_name.c_str(),
                     symmetric_matrix_name.c_str());
              printf("==================================================================================================================================================\n");
              printf("Problem            Status   Stat. Rad.    Optimal     Objective    Iterations    Inn. Iters.   Func. Eval.   Grad. Eval.      Time     NonOpt Time\n");
              printf("==================================================================================================================================================\n");

              // Declare optimal value
              double optimal_value;

              // Loop through test problems
              for (int problem_count = 0; problem_count < 21; problem_count++) {

                // Delete reports
                nonopt.reporter()->deleteReports();

                // Declare file report
                std::shared_ptr<FileReport> r(new FileReport("f", R_NL, R_PER_INNER_ITERATION));

                // Declare output file name
                char out_file[100];

                // Declare problem dimension
                int const dimension = 10;

                // Switch on problems
                switch (problem_count) {
                case 0:
                  problem = std::make_shared<ActiveFaces>(dimension);
                  printf("ActiveFaces        ");
                  strcpy(out_file, (char*)"output/ActiveFaces");
                  optimal_value = 0.0;
                  break;
                case 1:
                  problem = std::make_shared<BrownFunction_2>(dimension);
                  printf("BrownFunction_2    ");
                  strcpy(out_file, (char*)"output/BrownFunction_2");
                  optimal_value = 0.0;
                  break;
                case 2:
                  problem = std::make_shared<ChainedCB3_1>(dimension);
                  printf("ChainedCB3_1       ");
                  strcpy(out_file, (char*)"output/ChainedCB3_1");
                  optimal_value = ((double)dimension - 1.0) * 2.0;
                  break;
                case 3:
                  problem = std::make_shared<ChainedCB3_2>(dimension);
                  printf("ChainedCB3_2       ");
                  strcpy(out_file, (char*)"output/ChainedCB3_2");
                  optimal_value = ((double)dimension - 1.0) * 2.0;
                  break;
                case 4:
                  problem = std::make_shared<ChainedCrescent_1>(dimension);
                  printf("ChainedCrescent_1  ");
                  strcpy(out_file, (char*)"output/ChainedCrescent_1");
                  optimal_value = 0.0;
                  break;
                case 5:
                  problem = std::make_shared<ChainedCrescent_2>(dimension);
                  printf("ChainedCrescent_2  ");
                  strcpy(out_file, (char*)"output/ChainedCrescent_2");
                  optimal_value = 0.0;
                  break;
                case 6:
                  problem = std::make_shared<ChainedLQ>(dimension);
                  printf("ChainedLQ          ");
                  strcpy(out_file, (char*)"output/ChainedLQ");
                  optimal_value = -((double)dimension - 1.0) * sqrt(2.0);
                  break;
                case 7:
                  problem = std::make_shared<ChainedMifflin_2>(dimension);
                  printf("ChainedMifflin_2   ");
                  strcpy(out_file, (char*)"output/ChainedMifflin_2");
                  optimal_value = -6.51;
                  break;
                case 8:
                  problem = std::make_shared<MaxQ>(dimension);
                  printf("MaxQ               ");
                  strcpy(out_file, (char*)"output/MaxQ");
                  optimal_value = 0.0;
                  break;
                case 9:
                  problem = std::make_shared<MxHilb>(dimension);
                  printf("MxHilb             ");
                  strcpy(out_file, (char*)"output/MxHilb");
                  optimal_value = 0.0;
                  break;
                case 10:
                  problem = std::make_shared<QuadPoly>(dimension, 10, 5, 10.0);
                  printf("QuadPoly           ");
                  strcpy(out_file, (char*)"output/QuadPoly");
                  optimal_value = 0.0;
                  break;
                case 11:
                  problem = std::make_shared<Test29_11>(dimension);
                  printf("Test29_11          ");
                  strcpy(out_file, (char*)"output/Test29_11");
                  optimal_value = 49.838;
                  break;
                case 12:
                  problem = std::make_shared<Test29_13>(dimension);
                  printf("Test29_13          ");
                  strcpy(out_file, (char*)"output/Test29_13");
                  optimal_value = 4.538;
                  break;
                case 13:
                  problem = std::make_shared<Test29_17>(dimension);
                  printf("Test29_17          ");
                  strcpy(out_file, (char*)"output/Test29_17");
                  optimal_value = 0.0;
                  break;
                case 14:
                  problem = std::make_shared<Test29_19>(dimension);
                  printf("Test29_19          ");
                  strcpy(out_file, (char*)"output/Test29_19");
                  optimal_value = 0.0;
                  break;
                case 15:
                  problem = std::make_shared<Test29_2>(dimension);
                  printf("Test29_2           ");
                  strcpy(out_file, (char*)"output/Test29_2");
                  optimal_value = 0.0;
                  break;
                case 16:
                  problem = std::make_shared<Test29_20>(dimension);
                  printf("Test29_20          ");
                  strcpy(out_file, (char*)"output/Test29_20");
                  optimal_value = 0.0;
                  break;
                case 17:
                  problem = std::make_shared<Test29_22>(dimension);
                  printf("Test29_22          ");
                  strcpy(out_file, (char*)"output/Test29_22");
                  optimal_value = 0.0;
                  break;
                case 18:
                  problem = std::make_shared<Test29_24>(dimension);
                  printf("Test29_24          ");
                  strcpy(out_file, (char*)"output/Test29_24");
                  optimal_value = 0.0;
                  break;
                case 19:
                  problem = std::make_shared<Test29_5>(dimension);
                  printf("Test29_5           ");
                  strcpy(out_file, (char*)"output/Test29_5");
                  optimal_value = 0.0;
                  break;
                case 20:
                  problem = std::make_shared<Test29_6>(dimension);
                  printf("Test29_6           ");
                  strcpy(out_file, (char*)"output/Test29_6");
                  optimal_value = 0.0;
                  break;
                } // end switch

                // Set rest of output file name
                sprintf(out_file, "%s_%d_%d_%d_%d_%d_%d.out", out_file, direction_computation_number, inverse_hessian_update_number, line_search_number, point_set_update_number, qp_solver_number, symmetric_matrix_number);

                // Open output file
                r->open(out_file);

                // Add to reporter
                nonopt.reporter()->addReport(r);

                // Optimize
                nonopt.optimize(problem);

                // Print results
                printf("%6d  %+.4e  %+.4e  %+.4e  %12d  %12d  %12d  %12d  %+.4e  %+.4e\n",
                       nonopt.status(),
                       nonopt.stationarityRadius(),
                       optimal_value,
                       nonopt.objective(),
                       nonopt.iterations(),
                       nonopt.totalInnerIterations(),
                       nonopt.functionEvaluations(),
                       nonopt.gradientEvaluations(),
                       nonopt.time(),
                       nonopt.timeNonOpt());

              } // end for

              // Print new line
              printf("\n");

            } // end for

          } // end for

        } // end for

      } // end for

    } // end for

  } // end for

  // Return
  return 0;

} // end main
