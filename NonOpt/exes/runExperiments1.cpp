// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis


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

// Main function
int main(int argc, char* argv[])
{

	double optimal_value=0.0;

	  // Declare NonOptSolver
	  NonOptSolver nonopt;

	  // Delete reports
	  nonopt.reporter()->deleteReports();

  std::shared_ptr<Problem> problem;

  std::vector<std::string> problem_types;
    problem_types.push_back("exact");
     problem_types.push_back("inexact");
    problem_types.push_back("inexact_agg");


//    double sizes[]={1.0,2.0,3.0,5.0,8.0,20.0};
//    double Agg_size;

//for(int i=0;i<6;i++){
//	Agg_size=sizes[i];

  // Loop over settings
  for (int settings_number = 0; settings_number < problem_types.size(); settings_number++) {

	  std::string problem_type = problem_types[settings_number];
	    // Print header
	    printf("Starting %s problem\n", problem_type.c_str());
	    printf("================================================================================================================================================================\n");
	    printf("Problem            Status   Stat. Rad.    Optimal     Objective    Iterations    Inn. Iters.    QP Iters.    Func. Eval.   Grad. Eval.      Time     NonOpt Time\n");
	    printf("================================================================================================================================================================\n");

	      // Delete reports
	      nonopt.reporter()->deleteReports();

	      // Declare file report
	      std::shared_ptr<FileReport> r_NL(new FileReport("f", R_NL, R_PER_INNER_ITERATION));

	      char out_file_stub[100];
	      int const p_size=std::stoi(argv[2]);

	  if (strcmp(argv[1], "ActiveFaces") == 0) {
		problem = std::make_shared<ActiveFaces>(p_size);
		strcpy(out_file_stub, (char*)"output/ActiveFaces");
	  }
	  else if (strcmp(argv[1], "BrownFunction_2") == 0) {
		problem = std::make_shared<BrownFunction_2>(p_size);
		strcpy(out_file_stub, (char*)"output/BrownFunction_2");
	  }
	  else if (strcmp(argv[1], "ChainedCB3_1") == 0) {
		problem = std::make_shared<ChainedCB3_1>(p_size);
		 strcpy(out_file_stub, (char*)"output/ChainedCB3_1");
	  }
	  else if (strcmp(argv[1], "ChainedCB3_2") == 0) {
		problem = std::make_shared<ChainedCB3_2>(p_size);
		strcpy(out_file_stub, (char*)"output/ChainedCB3_2");
	  }
	  else if (strcmp(argv[1], "ChainedCrescent_1") == 0) {
		problem = std::make_shared<ChainedCrescent_1>(p_size);
		strcpy(out_file_stub, (char*)"output/ChainedCrescent_1");
	  }
	  else if (strcmp(argv[1], "ChainedCrescent_2") == 0) {
		problem = std::make_shared<ChainedCrescent_2>(p_size);
		strcpy(out_file_stub, (char*)"output/ChainedCrescent_2");
	  }
	  else if (strcmp(argv[1], "ChainedLQ") == 0) {
		problem = std::make_shared<ChainedLQ>(p_size);
		strcpy(out_file_stub, (char*)"output/ChainedLQ");
	  }
	  else if (strcmp(argv[1], "ChainedMifflin_2") == 0) {
		problem = std::make_shared<ChainedMifflin_2>(p_size);
		strcpy(out_file_stub, (char*)"output/ChainedMifflin_2");
	  }
	  else if (strcmp(argv[1], "MaxQ") == 0) {
		problem = std::make_shared<MaxQ>(p_size);
		strcpy(out_file_stub, (char*)"output/MaxQ");
	  }
	  else if (strcmp(argv[1], "MxHilb") == 0) {
		problem = std::make_shared<MxHilb>(p_size);
		strcpy(out_file_stub, (char*)"output/MxHilb");
	  }
	  else if (strcmp(argv[1], "Test29_2") == 0) {
		problem = std::make_shared<Test29_2>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_2");
	  }
	  else if (strcmp(argv[1], "Test29_5") == 0) {
		problem = std::make_shared<Test29_5>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_5");
	  }
	  else if (strcmp(argv[1], "Test29_6") == 0) {
		problem = std::make_shared<Test29_6>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_6");
	  }
	  else if (strcmp(argv[1], "Test29_11") == 0) {
		problem = std::make_shared<Test29_11>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_11");
	  }
	  else if (strcmp(argv[1], "Test29_13") == 0) {
		problem = std::make_shared<Test29_13>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_13");
	  }
	  else if (strcmp(argv[1], "Test29_17") == 0) {
		problem = std::make_shared<Test29_17>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_17");
	  }
	  else if (strcmp(argv[1], "Test29_19") == 0) {
		problem = std::make_shared<Test29_19>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_19");
	  }
	  else if (strcmp(argv[1], "Test29_20") == 0) {
		problem = std::make_shared<Test29_20>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_20");
	  }
	  else if (strcmp(argv[1], "Test29_22") == 0) {
		problem = std::make_shared<Test29_22>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_22");
	  }
	  else if (strcmp(argv[1], "Test29_24") == 0) {
		problem = std::make_shared<Test29_24>(p_size);
		strcpy(out_file_stub, (char*)"output/Test29_24");
	  }
	  else {
		problem = std::make_shared<QuadPoly>(50, 10, 5, 10.0);
		strcpy(out_file_stub, (char*)"output/QuadPoly");
	  }


	  // Modify direction computation option
	//  nonopt.options()->modifyStringValue(nonopt.reporter(), "direction_computation", argv[2]);
	//  nonopt.options()->modifyStringValue(nonopt.reporter(), "inverse_hessian_update", argv[3]);
	//  nonopt.options()->modifyStringValue(nonopt.reporter(), "line_search", argv[4]);
	//  nonopt.options()->modifyStringValue(nonopt.reporter(), "point_set_update", argv[5]);
	//  nonopt.options()->modifyStringValue(nonopt.reporter(), "qp_solver", argv[6]);
	//  nonopt.options()->modifyStringValue(nonopt.reporter(), "symmetric_matrix", argv[7]);


	  nonopt.options()->modifyBoolValue(nonopt.reporter(), "QPAS_allow_skip_inexact", true);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "QPAS_inexact_termination_ratio_min", 0.01);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "cpu_time_limit", 1000.0);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "PSP_size_factor", 10.0);
	  nonopt.options()->modifyIntegerValue(nonopt.reporter(), "gradient_evaluation_limit", 3e+05);

	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "QPAS_skip_factor", 0.25);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "PSP_envelope_factor", 2000);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCGC_random_sample_fraction", 5.0);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCAgg_random_sample_fraction", 3.0);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCGC_shortened_stepsize", 0.2);
	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCAgg_shortened_stepsize", 0.2);
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "inexact_termination_factor_initial", 1.0);
      nonopt.options()->modifyDoubleValue(nonopt.reporter(), "inexact_termination_update_factor", std::stod(argv[3]));

	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "LSWW_stepsize_sufficient_decrease_threshold", 1e-12);

	  nonopt.options()->modifyDoubleValue(nonopt.reporter(), "DCAgg_size_factor", 10.0);

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



      char out_file_NL[100];
      //char out_file_QP[100];

      // Set rest of output file name
      sprintf(out_file_NL, "%s_%s_0.25_2000_5.0_0.2_10.0_1.0_1.0_0.01_%s.out", out_file_stub,problem_type.c_str(),argv[3]);
      //sprintf(out_file_QP, "%s_%s_qp.out", out_file_stub, settings_name.c_str());

      // Open output file
      r_NL->open(out_file_NL);
      //r_QP->open(out_file_QP);

      // Add to reporter
      nonopt.reporter()->addReport(r_NL);
      //nonopt.reporter()->addReport(r_QP);

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

  }
	  // Return

  return 0;
}  // end main
