// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cstring>
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
int main() {

  // Declare solver object
  NonOptSolver nonopt;

  // Modify options from file
  nonopt.options()->modifyOptionsFromFile("nonopt.opt");

  // Set print level to 0
  nonopt.options()->modifyIntegerValue("print_level", 0);

  // Declare problem pointer
  std::shared_ptr<Problem> problem;

  // Declare vectors of strings for strategies
  std::vector<std::string> direction_computation_names;
  std::vector<std::string> approximate_hessian_update_names;
  std::vector<std::string> line_search_names;
  std::vector<std::string> point_set_update_names;
  std::vector<std::string> qp_solver_names;
  std::vector<std::string> symmetric_matrix_names;
  std::vector<std::string> termination_names;

  // Push names of strategies
  direction_computation_names.push_back("CuttingPlane");
  direction_computation_names.push_back("GradientCombination");
  direction_computation_names.push_back("Gradient");
  approximate_hessian_update_names.push_back("BFGS");
  approximate_hessian_update_names.push_back("DFP");
  line_search_names.push_back("WeakWolfe");
  line_search_names.push_back("Backtracking");
  point_set_update_names.push_back("Proximity");
  qp_solver_names.push_back("DualActiveSet");
  qp_solver_names.push_back("InteriorPoint");
  symmetric_matrix_names.push_back("Dense");
  symmetric_matrix_names.push_back("LimitedMemory");
  termination_names.push_back("Basic");
  termination_names.push_back("SecondQP");

  // Loop through direction computations
  for (int direction_computation_number = 0; direction_computation_number < (int)direction_computation_names.size(); direction_computation_number++) {

    // Set direction computation strategy
    nonopt.options()->modifyStringValue("direction_computation", direction_computation_names[direction_computation_number]);

    // Loop through approximate Hessian update strategies
    for (int approximate_hessian_update_number = 0; approximate_hessian_update_number < (int)approximate_hessian_update_names.size(); approximate_hessian_update_number++) {

      // Set approximate Hessian update strategy
      nonopt.options()->modifyStringValue("approximate_hessian_update", approximate_hessian_update_names[approximate_hessian_update_number]);

      // Loop through line searches
      for (int line_search_number = 0; line_search_number < (int)line_search_names.size(); line_search_number++) {

        // Set line search strategy
        nonopt.options()->modifyStringValue("line_search", line_search_names[line_search_number]);

        // Loop through point set updates
        for (int point_set_update_number = 0; point_set_update_number < (int)point_set_update_names.size(); point_set_update_number++) {

          // Set point set update strategy
          nonopt.options()->modifyStringValue("point_set_update", point_set_update_names[point_set_update_number]);

          // Loop through small scale qp solvers
          for (int qp_solver_small_scale_number = 0; qp_solver_small_scale_number < (int)qp_solver_names.size(); qp_solver_small_scale_number++) {

            // Set small scale qp solver strategy
            nonopt.options()->modifyStringValue("qp_solver_small_scale", qp_solver_names[qp_solver_small_scale_number]);

            // Loop through small scale qp solvers
            for (int qp_solver_large_scale_number = 0; qp_solver_large_scale_number < (int)qp_solver_names.size(); qp_solver_large_scale_number++) {

              // Set small scale qp solver strategy
              nonopt.options()->modifyStringValue("qp_solver_large_scale", qp_solver_names[qp_solver_large_scale_number]);

              // Loop through symmetric matrices
              for (int symmetric_matrix_number = 0; symmetric_matrix_number < (int)symmetric_matrix_names.size(); symmetric_matrix_number++) {

                // Set symmetric matrix strategy
                nonopt.options()->modifyStringValue("symmetric_matrix", symmetric_matrix_names[symmetric_matrix_number]);

                // Loop through termination
                for (int termination_number = 0; termination_number < (int)termination_names.size(); termination_number++) {

                  // Set termination strategy
                  nonopt.options()->modifyStringValue("termination", termination_names[termination_number]);

                  // Declare strategy names
                  std::string direction_computation_name;
                  std::string approximate_hessian_update_name;
                  std::string line_search_name;
                  std::string point_set_update_name;
                  std::string qp_solver_small_scale_name;
                  std::string qp_solver_large_scale_name;
                  std::string symmetric_matrix_name;
                  std::string termination_name;

                  // Read integer options
                  nonopt.options()->valueAsString("direction_computation", direction_computation_name);
                  nonopt.options()->valueAsString("approximate_hessian_update", approximate_hessian_update_name);
                  nonopt.options()->valueAsString("line_search", line_search_name);
                  nonopt.options()->valueAsString("point_set_update", point_set_update_name);
                  nonopt.options()->valueAsString("qp_solver_small_scale", qp_solver_small_scale_name);
                  nonopt.options()->valueAsString("qp_solver_large_scale", qp_solver_large_scale_name);
                  nonopt.options()->valueAsString("symmetric_matrix", symmetric_matrix_name);
                  nonopt.options()->valueAsString("termination", termination_name);

                  // Print header
                  printf("%s, %s, %s, %s, %s, %s, %s, %s\n",
                         direction_computation_name.c_str(),
                         approximate_hessian_update_name.c_str(),
                         line_search_name.c_str(),
                         point_set_update_name.c_str(),
                         qp_solver_small_scale_name.c_str(),
                         qp_solver_large_scale_name.c_str(),
                         symmetric_matrix_name.c_str(),
                         termination_name.c_str());
                  printf("=====================================================================================================================================\n");
                  printf("Problem            Status   Stat. Rad.   Objective    Iterations    Inn. Iters.   Func. Eval.   Grad. Eval.      Time     NonOpt Time\n");
                  printf("=====================================================================================================================================\n");

                  // Loop through test problems
                  for (int problem_count = 0; problem_count < 21; problem_count++) {

                    // Delete reports
                    nonopt.reporter()->deleteReports();

                    // Declare file report
                    std::shared_ptr<FileReport> r(new FileReport("f", R_NL, R_BASIC));

                    // Declare output file name
                    char out_file[200];

                    // Declare problem dimension
                    int const dimension = 100;

                    // Switch on problems
                    switch (problem_count) {
                    case 0:
                      problem = std::make_shared<ActiveFaces>(dimension);
                      printf("ActiveFaces        ");
                      strcpy(out_file, (char *)"output/ActiveFaces");
                      break;
                    case 1:
                      problem = std::make_shared<BrownFunction_2>(dimension);
                      printf("BrownFunction_2    ");
                      strcpy(out_file, (char *)"output/BrownFunction_2");
                      break;
                    case 2:
                      problem = std::make_shared<ChainedCB3_1>(dimension);
                      printf("ChainedCB3_1       ");
                      strcpy(out_file, (char *)"output/ChainedCB3_1");
                      break;
                    case 3:
                      problem = std::make_shared<ChainedCB3_2>(dimension);
                      printf("ChainedCB3_2       ");
                      strcpy(out_file, (char *)"output/ChainedCB3_2");
                      break;
                    case 4:
                      problem = std::make_shared<ChainedCrescent_1>(dimension);
                      printf("ChainedCrescent_1  ");
                      strcpy(out_file, (char *)"output/ChainedCrescent_1");
                      break;
                    case 5:
                      problem = std::make_shared<ChainedCrescent_2>(dimension);
                      printf("ChainedCrescent_2  ");
                      strcpy(out_file, (char *)"output/ChainedCrescent_2");
                      break;
                    case 6:
                      problem = std::make_shared<ChainedLQ>(dimension);
                      printf("ChainedLQ          ");
                      strcpy(out_file, (char *)"output/ChainedLQ");
                      break;
                    case 7:
                      problem = std::make_shared<ChainedMifflin_2>(dimension);
                      printf("ChainedMifflin_2   ");
                      strcpy(out_file, (char *)"output/ChainedMifflin_2");
                      break;
                    case 8:
                      problem = std::make_shared<MaxQ>(dimension);
                      printf("MaxQ               ");
                      strcpy(out_file, (char *)"output/MaxQ");
                      break;
                    case 9:
                      problem = std::make_shared<MxHilb>(dimension);
                      printf("MxHilb             ");
                      strcpy(out_file, (char *)"output/MxHilb");
                      break;
                    case 10:
                      problem = std::make_shared<QuadPoly>(dimension, 10, 5, 10.0, 0);
                      printf("QuadPoly           ");
                      strcpy(out_file, (char *)"output/QuadPoly");
                      break;
                    case 11:
                      problem = std::make_shared<Test29_11>(dimension);
                      printf("Test29_11          ");
                      strcpy(out_file, (char *)"output/Test29_11");
                      break;
                    case 12:
                      problem = std::make_shared<Test29_13>(dimension);
                      printf("Test29_13          ");
                      strcpy(out_file, (char *)"output/Test29_13");
                      break;
                    case 13:
                      problem = std::make_shared<Test29_17>(dimension);
                      printf("Test29_17          ");
                      strcpy(out_file, (char *)"output/Test29_17");
                      break;
                    case 14:
                      problem = std::make_shared<Test29_19>(dimension);
                      printf("Test29_19          ");
                      strcpy(out_file, (char *)"output/Test29_19");
                      break;
                    case 15:
                      problem = std::make_shared<Test29_2>(dimension);
                      printf("Test29_2           ");
                      strcpy(out_file, (char *)"output/Test29_2");
                      break;
                    case 16:
                      problem = std::make_shared<Test29_20>(dimension);
                      printf("Test29_20          ");
                      strcpy(out_file, (char *)"output/Test29_20");
                      break;
                    case 17:
                      problem = std::make_shared<Test29_22>(dimension);
                      printf("Test29_22          ");
                      strcpy(out_file, (char *)"output/Test29_22");
                      break;
                    case 18:
                      problem = std::make_shared<Test29_24>(dimension);
                      printf("Test29_24          ");
                      strcpy(out_file, (char *)"output/Test29_24");
                      break;
                    case 19:
                      problem = std::make_shared<Test29_5>(dimension);
                      printf("Test29_5           ");
                      strcpy(out_file, (char *)"output/Test29_5");
                      break;
                    case 20:
                      problem = std::make_shared<Test29_6>(dimension);
                      printf("Test29_6           ");
                      strcpy(out_file, (char *)"output/Test29_6");
                      break;
                    } // end switch

                    // Set rest of output file name
                    snprintf(out_file, 200, "%s_%s_%s_%s_%s_%s_%s_%s_%s.out", out_file, direction_computation_name.c_str(), approximate_hessian_update_name.c_str(), line_search_name.c_str(), point_set_update_name.c_str(), qp_solver_small_scale_name.c_str(), qp_solver_large_scale_name.c_str(), symmetric_matrix_name.c_str(), termination_name.c_str());

                    // Open output file
                    r->open(out_file);

                    // Add to reporter
                    nonopt.reporter()->addReport(r);

                    // Optimize
                    nonopt.optimize(problem);

                    // Print results
                    printf("%6d  %+.4e  %+.4e  %12d  %12d  %12d  %12d  %+.4e  %+.4e\n",
                           nonopt.status(),
                           nonopt.stationarityRadius(),
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

    } // end for

  } // end for

  // Return
  return 0;

} // end main
