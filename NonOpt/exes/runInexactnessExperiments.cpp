// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis and Minhan Li

#include <cmath>
#include <cstring>

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

using namespace NonOpt;

int main()
{

  // Declare NonOptSolver
  NonOptSolver nonopt;

  // Delete reports
  nonopt.reporter()->deleteReports();

  // Declare problem pointer
  std::shared_ptr<Problem> problem;

  // Declare settings to run
  std::vector<std::string> settings_names;
  settings_names.push_back("exact");
  settings_names.push_back("inexact");

  // Loop over settings
  for (int settings_number = 0; settings_number < settings_names.size(); settings_number++) {

    // Set settings name
    std::string settings_name = settings_names[settings_number];

    // Print header
    printf("Start testing %s\n", settings_name.c_str());
    printf("================================================================================================================================================================\n");
    printf("Problem            Status   Stat. Rad.    Optimal     Objective    Iterations    Inn. Iters.    QP Iters.    Func. Eval.   Grad. Eval.      Time     NonOpt Time\n");
    printf("================================================================================================================================================================\n");

    // Declare optimal value
    double optimal_value;

    // Loop through test problems
    for (int problem_count = 0; problem_count < 21; problem_count++) {

      // Delete reports
      nonopt.reporter()->deleteReports();

      // Declare file report
      std::shared_ptr<FileReport> r_NL(new FileReport("f", R_NL, R_PER_INNER_ITERATION));
      std::shared_ptr<FileReport> r_QP(new FileReport("f", R_QP, R_PER_INNER_ITERATION));

      // Declare output file name stub
      char out_file_stub[100];

      // Switch on problems
      switch (problem_count) {
        case 0:
          problem = std::make_shared<ActiveFaces>(50);
          printf("ActiveFaces        ");
          strcpy(out_file_stub, (char*)"output/ActiveFaces");
          optimal_value = 0.0;
          break;
        case 1:
          problem = std::make_shared<BrownFunction_2>(50);
          printf("BrownFunction_2    ");
          strcpy(out_file_stub, (char*)"output/BrownFunction_2");
          optimal_value = 0.0;
          break;
        case 2:
          problem = std::make_shared<ChainedCB3_1>(50);
          printf("ChainedCB3_1       ");
          strcpy(out_file_stub, (char*)"output/ChainedCB3_1");
          optimal_value = 49.0 * 2.0;
          break;
        case 3:
          problem = std::make_shared<ChainedCB3_2>(50);
          printf("ChainedCB3_2       ");
          strcpy(out_file_stub, (char*)"output/ChainedCB3_2");
          optimal_value = 49.0 * 2.0;
          break;
        case 4:
          problem = std::make_shared<ChainedCrescent_1>(50);
          printf("ChainedCrescent_1  ");
          strcpy(out_file_stub, (char*)"output/ChainedCrescent_1");
          optimal_value = 0.0;
          break;
        case 5:
          problem = std::make_shared<ChainedCrescent_2>(50);
          printf("ChainedCrescent_2  ");
          strcpy(out_file_stub, (char*)"output/ChainedCrescent_2");
          optimal_value = 0.0;
          break;
        case 6:
          problem = std::make_shared<ChainedLQ>(50);
          printf("ChainedLQ          ");
          strcpy(out_file_stub, (char*)"output/ChainedLQ");
          optimal_value = -49.0 * sqrt(2.0);
          break;
        case 7:
          problem = std::make_shared<ChainedMifflin_2>(50);
          printf("ChainedMifflin_2   ");
          strcpy(out_file_stub, (char*)"output/ChainedMifflin_2");
          optimal_value = -34.795;
          break;
        case 8:
          problem = std::make_shared<MaxQ>(50);
          printf("MaxQ               ");
          strcpy(out_file_stub, (char*)"output/MaxQ");
          optimal_value = 0.0;
          break;
        case 9:
          problem = std::make_shared<MxHilb>(50);
          printf("MxHilb             ");
          strcpy(out_file_stub, (char*)"output/MxHilb");
          optimal_value = 0.0;
          break;
        case 10:
          problem = std::make_shared<QuadPoly>(50, 10, 5, 10.0);
          printf("QuadPoly           ");
          strcpy(out_file_stub, (char*)"output/QuadPoly");
          optimal_value = 0.0;
          break;
        case 11:
          problem = std::make_shared<Test29_2>(50);
          printf("Test29_2           ");
          strcpy(out_file_stub, (char*)"output/Test29_2");
          optimal_value = 0.0;
          break;
        case 12:
          problem = std::make_shared<Test29_5>(50);
          printf("Test29_5           ");
          strcpy(out_file_stub, (char*)"output/Test29_5");
          optimal_value = 0.0;
          break;
        case 13:
          problem = std::make_shared<Test29_6>(50);
          printf("Test29_6           ");
          strcpy(out_file_stub, (char*)"output/Test29_6");
          optimal_value = 0.0;
          break;
        case 14:
          problem = std::make_shared<Test29_11>(50);
          printf("Test29_11          ");
          strcpy(out_file_stub, (char*)"output/Test29_11");
          optimal_value = 0.0;
          break;
        case 15:
          problem = std::make_shared<Test29_13>(50);
          printf("Test29_13          ");
          strcpy(out_file_stub, (char*)"output/Test29_13");
          optimal_value = 0.0;
          break;
        case 16:
          problem = std::make_shared<Test29_17>(50);
          printf("Test29_17          ");
          strcpy(out_file_stub, (char*)"output/Test29_17");
          optimal_value = 0.0;
          break;
        case 17:
          problem = std::make_shared<Test29_19>(50);
          printf("Test29_19          ");
          strcpy(out_file_stub, (char*)"output/Test29_19");
          optimal_value = 0.0;
          break;
        case 18:
          problem = std::make_shared<Test29_20>(50);
          printf("Test29_20          ");
          strcpy(out_file_stub, (char*)"output/Test29_20");
          optimal_value = 0.0;
          break;
        case 19:
          problem = std::make_shared<Test29_22>(50);
          printf("Test29_22          ");
          strcpy(out_file_stub, (char*)"output/Test29_22");
          optimal_value = 0.0;
          break;
        case 20:
          problem = std::make_shared<Test29_24>(50);
          printf("Test29_24          ");
          strcpy(out_file_stub, (char*)"output/Test29_24");
          optimal_value = 0.0;
          break;
      }  // end switch

      // Declare strings for output file handles
      char out_file_NL[100];
      char out_file_QP[100];

      // Set rest of output file name
      sprintf(out_file_NL, "%s_%s.out", out_file_stub, settings_name.c_str());
      sprintf(out_file_QP, "%s_%s_qp.out", out_file_stub, settings_name.c_str());

      // Open output file
      r_NL->open(out_file_NL);
      r_QP->open(out_file_QP);

      // Add to reporter
      nonopt.reporter()->addReport(r_NL);
      nonopt.reporter()->addReport(r_QP);

      // Read options
      nonopt.options()->modifyOptionsFromFile(nonopt.reporter(), settings_names[settings_number] + ".txt");

      // Optimize
      nonopt.optimize(problem);

      // Print results
      printf("%6d  %+.4e  %+.4e  %+.4e  %12d  %12d  %12d  %12d  %12d  %+.4e  %+.4e\n",
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

    }  // end for (over problems)

    // Print new line
    printf("\n");

  }  // end for (over settings)

  // Return
  return 0;

}  // end main
