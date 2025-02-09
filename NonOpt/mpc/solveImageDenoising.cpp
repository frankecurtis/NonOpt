// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <iostream>
#include <fstream>

#include "ImageDenoising.hpp"
#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"

using namespace NonOpt;

// Main function
int main(int argc, char* argv[])
{

  // Declare image dimensions
  int rows = 300;
  int cols = 400;

  // Declare parameters
  double regularization = 6e-02;
  double smoothing      = 1e-08;

  // Declare input stream
  std::ifstream inputFile("croissant_matrix.txt");
  if (!inputFile.is_open()) {
    std::cerr << "Could not open the file." << std::endl;
    return 1;
  }
  
  // Read matrix
  int temp;
  double* image = new double[rows*cols];
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      inputFile >> temp;
      image[r*cols + c] = ((double)temp) / 256.0;
    }
  }
  inputFile.close();

  // Declare problem
  std::shared_ptr<Problem> problem = std::make_shared<ImageDenoising>(rows,cols,image,regularization,smoothing);

  // Declare solver object
  NonOptSolver nonopt;

  // Modify options from file
  nonopt.options()->modifyOptionsFromFile("nonopt.opt");

  // Optimize
  nonopt.optimize(problem);

  // Get solution
  nonopt.solution(image);

  // Declare output stream
  std::ofstream outputFile("croissant_recovered_matrix.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Could not open output file." << std::endl;
    return 1;
  }

  // Print matrix
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      outputFile << (int)(image[r*cols + c] * 256.0) << " ";
    }
    outputFile << std::endl;
  }
  outputFile.close();

  // Delete array
  if (image != nullptr) {
    delete[] image;
    image = nullptr;
  } // end if

  // Return
  return 0;

} // end main
