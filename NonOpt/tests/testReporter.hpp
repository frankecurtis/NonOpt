// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTREPORTER_HPP__
#define __TESTREPORTER_HPP__

#include <fstream>
#include <iostream>

#include "NonOptReporter.hpp"

using namespace NonOpt;

// Implementation of test
int testReporterImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter r;

  // Declare stream report
  std::shared_ptr<StreamReport> rs(new StreamReport("s", R_NL, R_BASIC));

  // Set stream
  rs->setStream(&std::cout);

  // Check option
  if (option == 1) {
    // Add to reporter
    r.addReport(rs);
  }

  // Declare file report
  std::shared_ptr<FileReport> rf(new FileReport("f", R_NL, R_BASIC));

  // Open file
  rf->open("NonOpt_filereport_NL.txt");

  // Add to reporter
  r.addReport(rf);

  // Add file report directly
  r.addFileReport("g", "NonOpt_filereport_QP.txt", R_QP, R_BASIC);

  // Get report that was added directly
  std::shared_ptr<Report> rg = r.report("g");

  // Print something
  r.printf(R_NL, R_BASIC, "This line should appear in all reports.\n");
  r.printf(R_QP, R_BASIC, "This line should appear in all reports.\n");
  r.printf(R_NL, R_BASIC, "This is a string: %s\n", "NonOpt");
  r.printf(R_QP, R_BASIC, "This is a string: %s\n", "NonOpt");
  r.printf(R_NL, R_BASIC, "This is an integer: %d\n", 1);
  r.printf(R_QP, R_BASIC, "This is an integer: %d\n", 1);
  r.printf(R_NL, R_BASIC, "This is a float: %f\n", 2.3);
  r.printf(R_QP, R_BASIC, "This is a float: %f\n", 2.3);
  r.printf(R_NL, R_BASIC, "This is scientific notation: %e\n", 4.56);
  r.printf(R_QP, R_BASIC, "This is scientific notation: %e\n", 4.56);
  r.printf(R_NL, R_BASIC, "Here are all again: %s, %d, %f, %e\n", "NonOpt", 1, 2.3, 4.56);
  r.printf(R_QP, R_BASIC, "Here are all again: %s, %d, %f, %e\n", "NonOpt", 1, 2.3, 4.56);

  // Set types and levels
  rs->setTypeAndLevel(R_NL, R_PER_ITERATION);
  rf->setTypeAndLevel(R_NL, R_PER_INNER_ITERATION);
  rg->setTypeAndLevel(R_QP, R_PER_ITERATION);

  // Print stuff
  r.printf(R_NL, R_BASIC, "NL ALWAYS\n");
  r.printf(R_NL, R_PER_ITERATION, "NL PER ITERATION\n");
  r.printf(R_NL, R_PER_INNER_ITERATION, "NL PER INNER ITERATION\n");
  r.printf(R_QP, R_BASIC, "QP ALWAYS\n");
  r.printf(R_QP, R_PER_ITERATION, "QP PER ITERATION\n");
  r.printf(R_QP, R_PER_INNER_ITERATION, "QP PER INNER ITERATION\n");

  // Delete reports
  r.deleteReports();

  // Read NonOpt_filereport_NL.txt and check values
  std::ifstream infile("NonOpt_filereport_NL.txt");
  std::string line;
  std::getline(infile, line);
  if (line.compare("This line should appear in all reports.") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is a string: NonOpt") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is an integer: 1") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is a float: 2.300000") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("This is scientific notation: 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("Here are all again: NonOpt, 1, 2.300000, 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("NL ALWAYS") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("NL PER ITERATION") != 0) {
    result = 1;
  }
  std::getline(infile, line);
  if (line.compare("NL PER INNER ITERATION") != 0) {
    result = 1;
  }

  // Read NonOpt_filereport_QP.txt and check values
  std::ifstream infile2("NonOpt_filereport_QP.txt");
  std::string line2;
  std::getline(infile2, line2);
  if (line2.compare("This line should appear in all reports.") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is a string: NonOpt") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is an integer: 1") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is a float: 2.300000") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("This is scientific notation: 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("Here are all again: NonOpt, 1, 2.300000, 4.560000e+00") != 0) {
    result = 1;
  }
  std::getline(infile2, line2);
  if (line2.compare("QP PER ITERATION") != 0) {
    result = 1;
  }

  // Delete files
  remove("NonOpt_filereport_NL.txt");
  remove("NonOpt_filereport_QP.txt");

  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      printf("TEST WAS SUCCESSFUL.\n");
    }
    else {
      printf("TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testReporterImplementation

#endif /* __TESTREPORTER_HPP__ */
