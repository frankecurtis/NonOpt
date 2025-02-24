// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Lara Zebiane

#include <fstream>
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

  // Declare problem dimension
  int const dimension = 1000;

  // Declare solver object
  NonOptSolver nonopt;

  // Set print level to 0
  nonopt.options()->modifyIntegerValue("print_level", 0);

  // Declare option name
  std::string option_name;

  // Declare problem pointer
  std::shared_ptr<Problem> problem;

  // Declare vectors of strings for strategies
  std::vector<std::string> direction_computation_names;
  std::vector<std::string> direction_computation_short_names;

  // Push names of strategies
  direction_computation_names.push_back("CuttingPlane");
  direction_computation_names.push_back("GradientCombination");
  direction_computation_names.push_back("Gradient");
  direction_computation_short_names.push_back("CP");
  direction_computation_short_names.push_back("GC");
  direction_computation_short_names.push_back("G ");

  // Loop through option files (speed vs. accuracy)
  for (int option_number = 0; option_number <= 1; option_number++) {

    // Initialize total time
    double total_time = 0.0;

    // Modify options from file
    if (option_number == 0) {
      nonopt.options()->modifyOptionsFromFile("options_speed.opt");
      option_name = "Speed";
    }
    else {
      nonopt.options()->modifyOptionsFromFile("options_accuracy.opt");
      option_name = "Accuracy";
    }

    // Declare output stream
    std::ofstream outputFile(option_name + ".tex");
    if (!outputFile.is_open()) {
      std::cerr << "Could not open output file." << std::endl;
      return 1;
    }

    // Print header
    printf("%s\n", option_name.c_str());
    printf("==================================================================================================================================================\n");
    printf("Problem            Status   Stat. Rad.    Optimal     Objective    Iterations    Inn. Iters.   Func. Eval.   Grad. Eval.      Time     NonOpt Time\n");
    printf("==================================================================================================================================================\n");

    // Declare problem temporary values
    std::string stub;
    std::string stubLaTeX;
    double optimal_value;
    bool optimal_known;

    // Loop through test problems
    for (int problem_count = 0; problem_count < 20; problem_count++) {

      // Switch on problems
      switch (problem_count) {
      case 0:
        problem = std::make_shared<ActiveFaces>(dimension);
        stub = "ActiveFaces";
        stubLaTeX = "ActiveFaces";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 1:
        problem = std::make_shared<BrownFunction_2>(dimension);
        stub = "BrownFunction_2";
        stubLaTeX = "BrownFunction\\_2";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 2:
        problem = std::make_shared<ChainedCB3_1>(dimension);
        stub = "ChainedCB3_1";
        stubLaTeX = "ChainedCB3\\_1";
        optimal_value = ((double)dimension - 1.0) * 2.0;
        optimal_known = true;
        break;
      case 3:
        problem = std::make_shared<ChainedCB3_2>(dimension);
        stub = "ChainedCB3_2";
        stubLaTeX = "ChainedCB3\\_2";
        optimal_value = ((double)dimension - 1.0) * 2.0;
        optimal_known = true;
        break;
      case 4:
        problem = std::make_shared<ChainedCrescent_1>(dimension);
        stub = "ChainedCrescent_1";
        stubLaTeX = "ChainedCrescent\\_1";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 5:
        problem = std::make_shared<ChainedCrescent_2>(dimension);
        stub = "ChainedCrescent_2";
        stubLaTeX = "ChainedCrescent\\_2";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 6:
        problem = std::make_shared<ChainedLQ>(dimension);
        stub = "ChainedLQ";
        stubLaTeX = "ChainedLQ";
        optimal_value = -((double)dimension - 1.0) * sqrt(2.0);
        optimal_known = true;
        break;
      case 7:
        problem = std::make_shared<ChainedMifflin_2>(dimension);
        stub = "ChainedMifflin_2";
        stubLaTeX = "ChainedMifflin\\_2";
        optimal_value = -706.55;
        optimal_known = true;
        break;
      case 8:
        problem = std::make_shared<MaxQ>(dimension);
        stub = "MaxQ";
        stubLaTeX = "MaxQ";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 9:
        problem = std::make_shared<MxHilb>(dimension);
        stub = "MxHilb";
        stubLaTeX = "MxHilb";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 10:
        problem = std::make_shared<Test29_2>(dimension);
        stub = "Test29_2";
        stubLaTeX = "Test29\\_2";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 11:
        problem = std::make_shared<Test29_5>(dimension);
        stub = "Test29_5";
        stubLaTeX = "Test29\\_5";
        optimal_value = 0.0;
        optimal_known = true;
        break;
      case 12:
        problem = std::make_shared<Test29_6>(dimension);
        stub = "Test29_6";
        stubLaTeX = "Test29\\_6";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 13:
        problem = std::make_shared<Test29_11>(dimension);
        stub = "Test29_11";
        stubLaTeX = "Test29\\_11";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 14:
        problem = std::make_shared<Test29_13>(dimension);
        stub = "Test29_13";
        stubLaTeX = "Test29\\_13";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 15:
        problem = std::make_shared<Test29_17>(dimension);
        stub = "Test29_17";
        stubLaTeX = "Test29\\_17";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 16:
        problem = std::make_shared<Test29_19>(dimension);
        stub = "Test29_19";
        stubLaTeX = "Test29\\_19";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 17:
        problem = std::make_shared<Test29_20>(dimension);
        stub = "Test29_20";
        stubLaTeX = "Test29\\_20";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 18:
        problem = std::make_shared<Test29_22>(dimension);
        stub = "Test29_22";
        stubLaTeX = "Test29\\_22";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      case 19:
        problem = std::make_shared<Test29_24>(dimension);
        stub = "Test29_24";
        stubLaTeX = "Test29\\_24";
        optimal_value = 0.0;
        optimal_known = false;
        break;
      } // end switch

      // Declare direction computation temporary value
      std::string direction_computation_name;

      // Loop through direction computations
      for (int direction_computation_number = 0; direction_computation_number < (int)direction_computation_names.size(); direction_computation_number++) {

        // Set direction computation strategy
        nonopt.options()->modifyStringValue("direction_computation", direction_computation_names[direction_computation_number]);

        // Read options
        nonopt.options()->valueAsString("direction_computation", direction_computation_name);

        // Delete reports
        nonopt.reporter()->deleteReports();

        // Declare file report
        std::shared_ptr<FileReport> r(new FileReport("f", R_NL, R_PER_ITERATION));

        // Set output file name
        std::string out_name = "output/" + stub + "_" + direction_computation_name + "_" + option_name + ".out";
        char out_file[200];
        strcpy(out_file, out_name.c_str());

        // Open output file
        r->open(out_file);

        // Add to reporter
        nonopt.reporter()->addReport(r);

        // Optimize
        nonopt.optimize(problem);

        // Print results
        printf("%17s  %6d  %+.4e",
               stub.c_str(),
               nonopt.status(),
               nonopt.stationarityRadius());

        // Check if optimal value known
        if (optimal_known) {
          printf("  %+.4e", optimal_value);
        }
        else {
          printf("  -----------");
        }

        // Print remaining results
        printf("  %+.4e  %12d  %12d  %12d  %12d  %+.4e  %+.4e\n",
               nonopt.objective(),
               nonopt.iterations(),
               nonopt.totalInnerIterations(),
               nonopt.functionEvaluations(),
               nonopt.gradientEvaluations(),
               nonopt.time(),
               nonopt.timeNonOpt());
        
        // Increment total time
        total_time += nonopt.timeNonOpt();
        
        // Print to LaTeX
        if (direction_computation_number == 0) {
          outputFile << stubLaTeX << " & ";
        }
        else {
          outputFile << " & ";
        }
        outputFile << direction_computation_short_names[direction_computation_number];
        outputFile << " & ";
        outputFile << nonopt.iterations();
        outputFile << " & ";
        outputFile << nonopt.functionEvaluations();
        outputFile << " & ";
        outputFile << nonopt.gradientEvaluations();
        outputFile << " & ";
        outputFile << std::showpos << std::scientific << nonopt.objective() << std::noshowpos;
        outputFile << " & ";
        outputFile << std::fixed << nonopt.time();
        outputFile << " \\\\" << std::endl;

      } // end for

      // Print line for LaTeX
      outputFile << "\\hline" << std::endl;

    } // end for

    // Print total time
    printf("\nTotal time = %+.4e \n\n", total_time);

    // Close output file
    outputFile.close();

  } // end for

  // Return
  return 0;

} // end main
