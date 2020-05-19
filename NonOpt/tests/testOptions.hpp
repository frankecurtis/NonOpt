// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __TESTOPTIONS_HPP__
#define __TESTOPTIONS_HPP__

#include <iostream>

#include "NonOptOptions.hpp"
#include "NonOptReporter.hpp"

using namespace NonOpt;

// Implementation of test
int testOptionsImplementation(int option)
{

  // Initialize output
  int result = 0;

  // Declare reporter
  Reporter reporter;

  // Check option
  if (option == 1) {

    // Declare stream report
    std::shared_ptr<StreamReport> sr(new StreamReport("s", R_NL, R_BASIC));

    // Set stream report to standard output
    sr->setStream(&std::cout);

    // Add stream report to reporter
    reporter.addReport(sr);

  } // end if

  // Declare option
  Options o;

  // Declare bool for testing
  bool temp = true;

  // Print empty option list
  reporter.printf(R_NL, R_BASIC, "Printing options list... should be empty:\n");
  o.print(&reporter);

  // Add bool option (correctly)
  reporter.printf(R_NL, R_BASIC, "Adding bool option... should be no error message:\n");
  temp = o.addBoolOption(&reporter, "b", true, "Correctly added bool option");

  // Check result
  if (temp == false) {
    result = 1;
  }

  // Add double option (correctly)
  reporter.printf(R_NL, R_BASIC, "Adding double option... should be no error message:\n");
  temp = o.addDoubleOption(&reporter, "d", 1e-2, 0.0, 1.0, "Correctly added double option");

  // Check result
  if (temp == false) {
    result = 1;
  }

  // Add integer option (correctly)
  reporter.printf(R_NL, R_BASIC, "Adding integer option... should be no error message:\n");
  temp = o.addIntegerOption(&reporter, "i", 1, 0, 2, "Correctly added integer option");

  // Check result
  if (temp == false) {
    result = 1;
  }

  // Add string option (correctly)
  reporter.printf(R_NL, R_BASIC, "Adding string option... should be no error message:\n");
  temp = o.addStringOption(&reporter, "s", "string", "Correctly added string option");

  // Check result
  if (temp == false) {
    result = 1;
  }

  // Add bool option (repeated name)
  reporter.printf(R_NL, R_BASIC, "Adding bool option... should be error (duplicate name):\n");
  temp = o.addBoolOption(&reporter, "b", false, "Incorrectly added option, repeated name");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add double option (repeated name)
  reporter.printf(R_NL, R_BASIC, "Adding double option... should be error (duplicate name):\n");
  temp = o.addDoubleOption(&reporter, "i", 1.0, 0.0, 2.0, "Incorrectly added option, repeated name");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add integer option (repeated name)
  reporter.printf(R_NL, R_BASIC, "Adding integer option... should be error (duplicate name):\n");
  temp = o.addIntegerOption(&reporter, "i", 1, 0, 2, "Incorrectly added option, repeated name");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add bool option (repeated name)
  reporter.printf(R_NL, R_BASIC, "Adding string option... should be error (duplicate name):\n");
  temp = o.addStringOption(&reporter, "b", "true", "Incorrectly added option, repeated name");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add integer option (bad bounds)
  reporter.printf(R_NL, R_BASIC, "Adding integer option... should be error (bad bounds):\n");
  temp = o.addIntegerOption(&reporter, "j", 1, 2, 0, "Incorrectly added option, bad bounds");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add double option (bad bounds)
  reporter.printf(R_NL, R_BASIC, "Adding double option... should be error (bad bounds):\n");
  temp = o.addDoubleOption(&reporter, "j", 1.0, 2.0, 0.0, "Incorrectly added option, bad bounds");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add integer option (value out of bounds)
  reporter.printf(R_NL, R_BASIC, "Adding integer option... should be error (value out of bounds):\n");
  temp = o.addIntegerOption(&reporter, "j", -1, 0, 2, "Incorrectly added option, value out of bounds");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add integer option (value out of bounds)
  reporter.printf(R_NL, R_BASIC, "Adding integer option... should be error (value out of bounds):\n");
  temp = o.addIntegerOption(&reporter, "j", 3, 0, 2, "Incorrectly added option, value out of bounds");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add double option (value out of bounds)
  reporter.printf(R_NL, R_BASIC, "Adding double option... should be error (value out of bounds):\n");
  temp = o.addDoubleOption(&reporter, "j", -1.0, 0.0, 2.0, "Incorrectly added option, value out of bounds");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Add double option (value out of bounds)
  reporter.printf(R_NL, R_BASIC, "Adding double option... should be error (value out of bounds):\n");
  temp = o.addDoubleOption(&reporter, "j", 3.0, 0.0, 2.0, "Incorrectly added option, value out of bounds");

  // Check result
  if (temp == true) {
    result = 1;
  }

  // Print option list
  reporter.printf(R_NL, R_BASIC, "Printing options list...:\n");
  o.print(&reporter);

  // Declare values
  bool b;
  double d;
  int i;
  std::string s;

  // Get values
  o.valueAsBool(&reporter, "b", b);
  o.valueAsDouble(&reporter, "d", d);
  o.valueAsInteger(&reporter, "i", i);
  o.valueAsString(&reporter, "s", s);

  // Check result
  if (b != true) {
    result = 1;
  }
  if (d <= 1e-2 - 1e-12 || d >= 1e-2 + 1e-12) {
    result = 1;
  }
  if (i != 1) {
    result = 1;
  }
  if (s.compare("string") != 0) {
    result = 1;
  }

  // Print values
  reporter.printf(R_NL, R_BASIC, "Printing values only...:\n");
  reporter.printf(R_NL, R_BASIC, "b = %d\nd = %f\ni = %d\ns = %s\n", b, d, i, s.c_str());

  // Modify bool value
  reporter.printf(R_NL, R_BASIC, "Modifying \'b\'... should be no error message (value now false):\n");
  temp = o.modifyBoolValue(&reporter, "b", false);
  o.valueAsBool(&reporter, "b", b);

  // Check value and result
  if (temp == false) {
    result = 1;
  }
  if (b != false) {
    result = 1;
  }

  // Modify double value
  reporter.printf(R_NL, R_BASIC, "Modifying \'d\'... should be error (value out of bounds):\n");
  temp = o.modifyDoubleValue(&reporter, "d", -1.0);
  o.valueAsDouble(&reporter, "d", d);

  // Check value and result
  if (temp == true) {
    result = 1;
  }
  if (d <= 1e-2 - 1e-12 || d >= 1e-2 + 1e-12) {
    result = 1;
  }

  // Modify double value
  reporter.printf(R_NL, R_BASIC, "Modifying \'d\'... should be no error message (value now 0.5):\n");
  temp = o.modifyDoubleValue(&reporter, "d", 0.5);
  o.valueAsDouble(&reporter, "d", d);

  // Check value and result
  if (temp == false) {
    result = 1;
  }
  if (d <= 0.5 - 1e-12 || d >= 0.5 + 1e-12) {
    result = 1;
  }

  // Modify double value
  reporter.printf(R_NL, R_BASIC, "Modifying \'d\'... should be error (called integer by mistake):\n");
  temp = o.modifyIntegerValue(&reporter, "d", 0);
  o.valueAsDouble(&reporter, "d", d);

  // Check value and result
  if (temp == true) {
    result = 1;
  }
  if (d <= 0.5 - 1e-12 || d >= 0.5 + 1e-12) {
    result = 1;
  }

  // Modify integer value
  reporter.printf(R_NL, R_BASIC, "Modifying \'i\'... should be error (value out of bounds):\n");
  temp = o.modifyIntegerValue(&reporter, "i", 3);
  o.valueAsInteger(&reporter, "i", i);

  // Check value and result
  if (temp == true) {
    result = 1;
  }
  if (i != 1) {
    result = 1;
  }

  // Modify integer value
  reporter.printf(R_NL, R_BASIC, "Modifying \'i\'... should be no error message (value now 2):\n");
  temp = o.modifyIntegerValue(&reporter, "i", 2);
  o.valueAsInteger(&reporter, "i", i);

  // Check value and result
  if (temp == false) {
    result = 1;
  }
  if (i != 2) {
    result = 1;
  }

  // Modify integer value
  reporter.printf(R_NL, R_BASIC, "Modifying \'i\'... should be error (called double by mistake):\n");
  temp = o.modifyDoubleValue(&reporter, "i", 1.0);
  o.valueAsInteger(&reporter, "i", i);

  // Check value and result
  if (temp == true) {
    result = 1;
  }
  if (i != 2) {
    result = 1;
  }

  // Modify string value
  reporter.printf(R_NL, R_BASIC, "Modifying \'s\'... should be no error message (value now \'characters\'):\n");
  temp = o.modifyStringValue(&reporter, "s", "characters");
  o.valueAsString(&reporter, "s", s);

  // Check value and result
  if (temp == false) {
    result = 1;
  }
  if (s.compare("characters") != 0) {
    result = 1;
  }

  // Print option list
  reporter.printf(R_NL, R_BASIC, "Printing modified options list...:\n");
  o.print(&reporter);

  // Modify options from file
  o.modifyOptionsFromFile(&reporter, "testOptions.opt");
  o.valueAsBool(&reporter, "b", b);
  o.valueAsDouble(&reporter, "d", d);
  o.valueAsInteger(&reporter, "i", i);
  o.valueAsString(&reporter, "s", s);

  // Check values
  if (b != true) {
    result = 1;
  }
  if (d != 0.9) {
    result = 1;
  }
  if (i != 0) {
    result = 1;
  }
  if (s.compare("letters") != 0) {
    result = 1;
  }

  // Check option
  if (option == 1) {
    // Print final message
    if (result == 0) {
      reporter.printf(R_NL, R_BASIC, "TEST WAS SUCCESSFUL.\n");
    }
    else {
      reporter.printf(R_NL, R_BASIC, "TEST FAILED.\n");
    }
  } // end if

  // Return
  return result;

} // end testOptionsImplementation

#endif /* __TESTOPTIONS_HPP__ */
