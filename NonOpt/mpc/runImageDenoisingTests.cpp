// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Lara Zebiane

#include <iostream>
#include <fstream>
#include <string>

#include "ImageDenoising.hpp"
#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"

using namespace NonOpt;

// Main function
int main(int argc, char *argv[]) {

  // Declare image dimensions
  int rows = 605;
  int cols = 807;

  // Store MSEs for final printing
  double mses[4];
  double times[4];

  // Loop over regularizers
  for (int regularizer = 0; regularizer <= 3; regularizer++) {

    // Declare parameters
    double regularization;
    double smoothing;

    // Set tuned parameters that minimize MSE
    switch (regularizer) {
    case 0:
      regularization = pow(2.0, 6.0); // TUNED!
      smoothing = pow(2.0, 0.0);
      break;
    case 1:
      regularization = pow(2.0, 21.0); // TUNED!
      smoothing = pow(2.0, -7.0);
      break;
    case 2:
      regularization = pow(2.0, 25.0); // TUNED!
      smoothing = pow(2.0, -19.0);
      break;
    case 3:
      regularization = pow(2.0, 6.0); // TUNED!
      smoothing = pow(2.0, 18.0);
    }

    // Declare input streams
    std::ifstream inputFile("croissant_matrix.txt");
    std::ifstream inputOriginalFile("croissant_original_matrix.txt");
    if (!inputFile.is_open() || !inputOriginalFile.is_open()) {
      std::cerr << "Could not open input file." << std::endl;
      return 1;
    }

    // Declare image and recovered_image
    double* image = new double[rows * cols];
    double* image_original = new double[rows * cols];
    double* image_recovered = new double[rows * cols];

    // Read matrix
    int temp;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        inputFile >> temp;
        image[r * cols + c] = (double)temp;
        inputOriginalFile >> temp;
        image_original[r * cols + c] = (double)temp;
      }
    }
    inputFile.close();

    // Declare problem
    std::shared_ptr<Problem> problem = std::make_shared<ImageDenoising>(rows, cols, image, regularizer, regularization, smoothing);

    // Declare solver object
    NonOptSolver nonopt;

    // Modify options from file
    nonopt.options()->modifyOptionsFromFile("options_image.opt");

    // Optimize
    nonopt.optimize(problem);

    // Get solution
    nonopt.solution(image_recovered);

    // Declare output stream
    std::ofstream outputFile("croissant_recovered_" + std::to_string(regularizer) + ".txt");
    if (!outputFile.is_open()) {
      std::cerr << "Could not open output file." << std::endl;
      return 1;
    }

    // Initialize mean-squared error
    double mse = 0.0;

    // Print matrix
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) {
        mse += pow(image_original[r * cols + c] - image_recovered[r * cols + c], 2.0);
        outputFile << (int)(image_recovered[r * cols + c] * 1.0) << " ";
      }
      outputFile << std::endl;
    }
    outputFile.close();

    // Print mean-squared error
    std::cout << std::endl << "MSE = " << mse << std::endl;

    // Store MSE and time
    mses[regularizer] = mse;
    times[regularizer] = nonopt.time();

    // Delete arrays
    if (image != nullptr) {
      delete[] image;
      image = nullptr;
    } // end if
    if (image_original != nullptr) {
      delete[] image_original;
      image_original = nullptr;
    } // end if
    if (image_recovered != nullptr) {
      delete[] image_recovered;
      image_recovered = nullptr;
    } // end if

  } // end for

  // Print MSEs
  std::cout << std::endl;
  for (int regularizer = 0; regularizer <= 3; regularizer++) {
    std::cout << "MSE [" << regularizer << "] = " << mses[regularizer] << std::endl;
    std::cout << "time[" << regularizer << "] = " << times[regularizer] << std::endl;
  }

  // Return
  return 0;

} // end main
