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

int main(int argc, char* argv[])
{
  // Declare NonOptSolver
  NonOptSolver nonopt;

  // Delete reports
  nonopt.reporter()->deleteReports();

  // Declare problem pointer
  std::shared_ptr<Problem> problem;


  std::vector<std::string> problem_types;
//    problem_types.push_back("exact");
//     problem_types.push_back("inexact");
    problem_types.push_back("inexact_agg");


//    double sizes[]={1e-10};//,1e-11,1e-12
//    double ls_size;
//
//for(int i=0;i<3;i++){
//	ls_size=sizes[i];

  // Loop over settings
  for (int settings_number = 0; settings_number < problem_types.size(); settings_number++) {

    // Set settings name
    std::string problem_type = problem_types[settings_number];
//    if( (settings_number==0||settings_number==1)&&(Agg_size!=10.0) ){
//    	continue;
//    }

    // Print header
    printf("Starting %s problem\n", problem_type.c_str());
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
      //std::shared_ptr<FileReport> r_QP(new FileReport("f", R_QP, R_PER_INNER_ITERATION));


      //int problem_count=std::stoi(argv[1]);
      // Declare output file name stub
      char out_file_stub[100];
      int const p_size=0;



      switch (problem_count) {
        case 0:
          problem = std::make_shared<ActiveFaces>(2830-p_size);
          printf("ActiveFaces        ");
          strcpy(out_file_stub, (char*)"output/ActiveFaces");
          optimal_value = 0.0;
          break;
        case 1:
          problem = std::make_shared<BrownFunction_2>(370-p_size);
          printf("BrownFunction_2    ");
          strcpy(out_file_stub, (char*)"output/BrownFunction_2");
          optimal_value = 0.0;
          break;
        case 2:
          problem = std::make_shared<ChainedCB3_1>(340-p_size);
          printf("ChainedCB3_1       ");
          strcpy(out_file_stub, (char*)"output/ChainedCB3_1");
          optimal_value = 49.0 * 2.0;
          break;
        case 3:
          problem = std::make_shared<ChainedCB3_2>(2860-p_size);
          printf("ChainedCB3_2       ");
          strcpy(out_file_stub, (char*)"output/ChainedCB3_2");
          optimal_value = 49.0 * 2.0;
          break;
        case 4:
          problem = std::make_shared<ChainedCrescent_1>(4300-p_size);
          printf("ChainedCrescent_1  ");
          strcpy(out_file_stub, (char*)"output/ChainedCrescent_1");
          optimal_value = 0.0;
          break;
        case 5:
          problem = std::make_shared<ChainedCrescent_2>(280-p_size);
          printf("ChainedCrescent_2  ");
          strcpy(out_file_stub, (char*)"output/ChainedCrescent_2");
          optimal_value = 0.0;
          break;
        case 6:
          problem = std::make_shared<ChainedLQ>(280-p_size);
          printf("ChainedLQ          ");
          strcpy(out_file_stub, (char*)"output/ChainedLQ");
          optimal_value = -49.0 * sqrt(2.0);
          break;
        case 7:
          problem = std::make_shared<ChainedMifflin_2>(190-p_size);
          printf("ChainedMifflin_2   ");
          strcpy(out_file_stub, (char*)"output/ChainedMifflin_2");
          optimal_value = -34.795;
          break;
        case 8:
          problem = std::make_shared<MaxQ>(700-p_size);
          printf("MaxQ               ");
          strcpy(out_file_stub, (char*)"output/MaxQ");
          optimal_value = 0.0;
          break;
        case 9:
          problem = std::make_shared<MxHilb>(940-p_size);
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
          problem = std::make_shared<Test29_2>(610-p_size);
          printf("Test29_2           ");
          strcpy(out_file_stub, (char*)"output/Test29_2");
          optimal_value = 0.0;
          break;
        case 12:
          problem = std::make_shared<Test29_5>(970-p_size);
          printf("Test29_5           ");
          strcpy(out_file_stub, (char*)"output/Test29_5");
          optimal_value = 0.0;
          break;
        case 13:
          problem = std::make_shared<Test29_6>(400-p_size);
          printf("Test29_6           ");
          strcpy(out_file_stub, (char*)"output/Test29_6");
          optimal_value = 0.0;
          break;
        case 14:
          problem = std::make_shared<Test29_11>(160-p_size);
          printf("Test29_11          ");
          strcpy(out_file_stub, (char*)"output/Test29_11");
          optimal_value = 0.0;
          break;
        case 15:
          problem = std::make_shared<Test29_13>(280-p_size);
          printf("Test29_13          ");
          strcpy(out_file_stub, (char*)"output/Test29_13");
          optimal_value = 0.0;
          break;
        case 16:
          problem = std::make_shared<Test29_17>(640-p_size);
          printf("Test29_17          ");
          strcpy(out_file_stub, (char*)"output/Test29_17");
          optimal_value = 0.0;
          break;
        case 17:
          problem = std::make_shared<Test29_19>(430-p_size);
          printf("Test29_19          ");
          strcpy(out_file_stub, (char*)"output/Test29_19");
          optimal_value = 0.0;
          break;
        case 18:
          problem = std::make_shared<Test29_20>(220-p_size);
          printf("Test29_20          ");
          strcpy(out_file_stub, (char*)"output/Test29_20");
          optimal_value = 0.0;
          break;
        case 19:
          problem = std::make_shared<Test29_22>(1120-p_size);
          printf("Test29_22          ");
          strcpy(out_file_stub, (char*)"output/Test29_22");
          optimal_value = 0.0;
          break;
        case 20:
          problem = std::make_shared<Test29_24>(100-p_size);
          printf("Test29_24          ");
          strcpy(out_file_stub, (char*)"output/Test29_24");
          optimal_value = 0.0;
          break;
      }  // end switch


       //Declare strings for output file handles
      char out_file_NL[100];
      //char out_file_QP[100];

      // Set rest of output file name
      sprintf(out_file_NL, "%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.out", out_file_stub,problem_type.c_str(),argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9]);
      //sprintf(out_file_QP, "%s_%s_qp.out", out_file_stub, settings_name.c_str());

      // Open output file
      r_NL->open(out_file_NL);
      //r_QP->open(out_file_QP);

      // Add to reporter
      nonopt.reporter()->addReport(r_NL);
      //nonopt.reporter()->addReport(r_QP);



	  nonopt.options()->modifyBoolValue(nonopt.reporter(), "QPAS_allow_skip_inexact", true);
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "QPAS_inexact_termination_ratio_min", 0.01);
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "cpu_time_limit", 300.0);
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "PSP_size_factor", std::stod(argv[5]));
      nonopt.options()->modifyIntegerValue(nonopt.reporter(), "gradient_evaluation_limit", 2e+05);

      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "QPAS_skip_factor", std::stod(argv[1]));
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "PSP_envelope_factor", std::stod(argv[2]));
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCGC_random_sample_fraction", std::stod(argv[3]));
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCAgg_random_sample_fraction", std::stod(argv[3]));
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCGC_shortened_stepsize", std::stod(argv[4]));
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCAgg_shortened_stepsize", std::stod(argv[4]));

      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "inexact_termination_factor_initial", std::stod(argv[6]));
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "inexact_termination_update_factor", std::stod(argv[7]));

      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCAgg_size_factor", std::stod(argv[9]));


	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "LSWW_stepsize_sufficient_decrease_threshold", std::stod(argv[8])*1e-10);

      if(strcmp(problem_type.c_str(), "exact") == 0){
    	  nonopt.options()->modifyBoolValue(nonopt.reporter(), "QPAS_allow_inexact_termination", false);
    	  nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", "GradientCombination");
      }
      else if(strcmp(problem_type.c_str(), "inexact") == 0){
    	  nonopt.options()->modifyBoolValue(nonopt.reporter(), "QPAS_allow_inexact_termination", true);
    	  nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", "GradientCombination");
      }
      else{
    	  nonopt.options()->modifyBoolValue(nonopt.reporter(), "QPAS_allow_inexact_termination", true);
    	  nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", "Aggregation");
      }


      // Optimize
      nonopt.optimize(problem);

      // Print results
      printf("%6d  %+.4e  %+.4e  %+.4e  %12d  %12d  %12d  %12d  %12d  %+.4e  %+12d\n",
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
			 nonopt.numberOfVariables());
             //nonopt.timeNonOpt());

    }  // end for (over problems)

    // Print new line
    printf("\n");
}
//}


  // Return
  return 0;

}  // end main
