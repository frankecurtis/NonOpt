// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <fstream>
#include <sstream>

#include "NonOptOptions.hpp"

namespace NonOpt
{

/////////////
// Options //
/////////////

// Options: Print
void Options::print(const Reporter* reporter) const
{

  // Print message if list is empty
  if (option_list_.size() == 0) {
    reporter->printf(R_NL, R_BASIC, "Option list is empty.\n");
    return;
  }

  // Print all options in list
  for (int i = 0; i < (int)option_list_.size(); i++) {
    option_list_[i]->print(reporter);
    if (i < (int)option_list_.size() - 1) {
      reporter->printf(R_NL, R_BASIC, "\n");
    }
  }

} // end print

// Options: Add bool option
bool Options::addBoolOption(const Reporter* reporter,
                            std::string name,
                            bool value,
                            std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      reporter->printf(R_NL, R_BASIC, "Option with name \'%s\' already exists.  Ignoring addition request.\n", name.c_str());
      return false;
    } // end if
  }   // end for

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "bool", value, description));

  // Add to list
  option_list_.push_back(option);

  // Return true
  return true;

} // end addBoolOption

// Options: Add double option
bool Options::addDoubleOption(const Reporter* reporter,
                              std::string name,
                              double value,
                              double lower_bound,
                              double upper_bound,
                              std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      reporter->printf(R_NL, R_BASIC, "Option with name \'%s\' already exists.  Ignoring addition request.\n", name.c_str());
      return false;
    } // end if
  }   // end for

  // Check that lower bound is lower than upper bound
  if (lower_bound > upper_bound) {
    reporter->printf(R_NL, R_BASIC, "Attempted to add double option \'%s\', but lower bound \'%+e\' greater than upper bound \'%+e\'.  Ignoring addition request.\n", name.c_str(), lower_bound, upper_bound);
    return false;
  } // end if

  // Check that value is within bounds
  if (value < lower_bound || value > upper_bound) {
    reporter->printf(R_NL, R_BASIC, "Attempted to add double option \'%s\', but value \'%+e\' outside of bound interval \'[%+e,%+e]\'.  Ignoring addition request.\n", name.c_str(), value, lower_bound, upper_bound);
    return false;
  } // end if

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "double", value, lower_bound, upper_bound, description));

  // Add to list
  option_list_.push_back(option);

  // Return true
  return true;

} // end addDoubleOption

// Options: Add integer option
bool Options::addIntegerOption(const Reporter* reporter,
                               std::string name,
                               int value,
                               int lower_bound,
                               int upper_bound,
                               std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      reporter->printf(R_NL, R_BASIC, "Option with name \'%s\' already exists.  Ignoring addition request.\n", name.c_str());
      return false;
    } // end if
  }   // end for

  // Check that lower bound is lower than upper bound
  if (lower_bound > upper_bound) {
    reporter->printf(R_NL, R_BASIC, "Attempted to add integer option \'%s\', but lower bound \'%d\' greater than upper bound \'%d\'.  Ignoring addition request.\n", name.c_str(), lower_bound, upper_bound);
    return false;
  } // end if

  // Check that value is within bounds
  if (value < lower_bound || value > upper_bound) {
    reporter->printf(R_NL, R_BASIC, "Attempted to add integer option \'%s\', but value \'%d\' outside of bound interval \'[%d,%d]\'.  Ignoring addition request.\n", name.c_str(), value, lower_bound, upper_bound);
    return false;
  } // end if

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "integer", value, lower_bound, upper_bound, description));

  // Add to list
  option_list_.push_back(option);

  // Return true
  return true;

} // end addIntegerOption

// Options: Add string option
bool Options::addStringOption(const Reporter* reporter,
                              std::string name,
                              std::string value,
                              std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      reporter->printf(R_NL, R_BASIC, "Option with name \'%s\' already exists.  Ignoring addition request.\n", name.c_str());
      return false;
    } // end if
  }   // end for

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "string", value, description));

  // Add to list
  option_list_.push_back(option);

  // Return true
  return true;

} // end addStringOption

// Options: Get lower bound as a double
bool Options::lowerBoundAsDouble(const Reporter* reporter,
                                 std::string name,
                                 double& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("double") == 0) {
        value = option_list_[i]->lowerBoundAsDouble();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access lower bound for option \'%s\' as double, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end lowerBoundAsDouble

// Options: Get lower bound as an integer
bool Options::lowerBoundAsInteger(const Reporter* reporter,
                                  std::string name,
                                  int& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("integer") == 0) {
        value = option_list_[i]->lowerBoundAsInteger();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access lower bound for option \'%s\' as integer, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end lowerBoundAsInteger

// Options: Get upper bound as a double
bool Options::upperBoundAsDouble(const Reporter* reporter,
                                 std::string name,
                                 double& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("double") == 0) {
        value = option_list_[i]->upperBoundAsDouble();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access upper bound for option \'%s\' as double, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end upperBoundAsDouble

// Options: Get upper bound as an integer
bool Options::upperBoundAsInteger(const Reporter* reporter,
                                  std::string name,
                                  int& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("integer") == 0) {
        value = option_list_[i]->upperBoundAsInteger();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access upper bound for option \'%s\' as integer, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end upperBoundAsInteger

// Options: Get value as a bool
bool Options::valueAsBool(const Reporter* reporter,
                          std::string name,
                          bool& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("bool") == 0) {
        value = option_list_[i]->valueAsBool();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access value for option \'%s\' as string, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsBool

// Options: Get value as a double
bool Options::valueAsDouble(const Reporter* reporter,
                            std::string name,
                            double& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("double") == 0) {
        value = option_list_[i]->valueAsDouble();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access value for option \'%s\' as double, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsDouble

// Options: Get value as an integer
bool Options::valueAsInteger(const Reporter* reporter,
                             std::string name,
                             int& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("integer") == 0) {
        value = option_list_[i]->valueAsInteger();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access value for option \'%s\' as integer, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsInteger

// Options: Get value as a string
bool Options::valueAsString(const Reporter* reporter,
                            std::string name,
                            std::string& value) const
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("string") == 0) {
        value = option_list_[i]->valueAsString();
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to access value for option \'%s\' as string, but type is %s.  Returning false.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsString

// Options: Modify from file nonopt.opt
void Options::modifyOptionsFromFile(const Reporter* reporter,
                                    std::string file_name)
{

  // Declare input file stream
  std::ifstream infile(file_name);

  // Declare line
  std::string line;

  // Loop through lines
  while (std::getline(infile, line)) {

    // Declare stream
    std::istringstream iss(line);

    // Declare substrings
    std::string s1, s2;

    // Grab first and second words in line
    iss >> s1 >> s2;

    // Declare bool for name found
    bool name_found = false;

    // Loop through to find value
    for (int i = 0; i < (int)option_list_.size(); i++) {
      if (option_list_[i]->name().compare(s1) == 0) {
        name_found = true;
        if (option_list_[i]->type().compare("bool") == 0) {
          try {
            modifyBoolValue(reporter, s1, (s2.compare("true") == 0 || s2.compare("1") == 0));
          } catch (...) {
            reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\', but cannot convert \'%s\' to bool.  Ignoring request.\n", s1.c_str(), s2.c_str());
          }
        }
        else if (option_list_[i]->type().compare("double") == 0) {
          try {
            modifyDoubleValue(reporter, s1, std::stod(s2));
          } catch (...) {
            reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\', but cannot convert \'%s\' to double.  Ignoring request.\n", s1.c_str(), s2.c_str());
          }
        } // end if
        else if (option_list_[i]->type().compare("integer") == 0) {
          try {
            modifyIntegerValue(reporter, s1, (int)std::stod(s2));
          } catch (...) {
            reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\', but cannot convert \'%s\' to int.  Ignoring request.\n", s1.c_str(), s2.c_str());
          }
        } // end else if
        else {
          modifyStringValue(reporter, s1, s2);
        }
      } // end if
    }   // end for

    // Print message if name not found
    if (!name_found) {
      reporter->printf(R_NL, R_BASIC, "Option with name \'%s\' does not exist.  Ignoring request.\n", s1.c_str());
    }

  } // end while

} // end modifyOptionsFromFile

// Options: Modify bool value
bool Options::modifyBoolValue(const Reporter* reporter,
                              std::string name,
                              bool value)
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("bool") == 0) {
        option_list_[i]->modifyBoolValue(value);
        reporter->printf(R_NL, R_BASIC, "Set value for option \'%s\' as %s.\n", name.c_str(), (value) ? "true" : "false");
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\' as bool, but type is %s.  Ignoring request.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyBoolValue

// Options: Modify double value
bool Options::modifyDoubleValue(const Reporter* reporter,
                                std::string name,
                                double value)
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("double") == 0) {
        if (value >= option_list_[i]->lowerBoundAsDouble() &&
            value <= option_list_[i]->upperBoundAsDouble()) {
          option_list_[i]->modifyDoubleValue(value);
          reporter->printf(R_NL, R_BASIC, "Set value for option \'%s\' as %+e.\n", name.c_str(), value);
          return true;
        } // end if
        else {
          reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\', but value %+e outside of bound interval \'[%+e,%+e]\'.  Ignoring request.\n", name.c_str(), value, option_list_[i]->lowerBoundAsDouble(), option_list_[i]->upperBoundAsDouble());
          return false;
        } // end else
      }   // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\' as double, but type is %s.  Ignoring request.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyDoubleValue

// Options: Modify integer value
bool Options::modifyIntegerValue(const Reporter* reporter,
                                 std::string name,
                                 int value)
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("integer") == 0) {
        if (value >= option_list_[i]->lowerBoundAsInteger() &&
            value <= option_list_[i]->upperBoundAsInteger()) {
          option_list_[i]->modifyIntegerValue(value);
          reporter->printf(R_NL, R_BASIC, "Set value for option \'%s\' as %d.\n", name.c_str(), value);
          return true;
        } // end if
        else {
          reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\', but value %d outside of bound interval \'[%d,%d]\'.  Ignoring request.\n", name.c_str(), value, option_list_[i]->lowerBoundAsInteger(), option_list_[i]->upperBoundAsInteger());
          return false;
        } // end else
      }   // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\' as integer, but type is %s.  Ignoring request.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyIntegerValue

// Options: Modify string value
bool Options::modifyStringValue(const Reporter* reporter,
                                std::string name,
                                std::string value)
{

  // Loop through to find value
  for (int i = 0; i < (int)option_list_.size(); i++) {
    if (option_list_[i]->name().compare(name) == 0) {
      if (option_list_[i]->type().compare("string") == 0) {
        option_list_[i]->modifyStringValue(value);
        reporter->printf(R_NL, R_BASIC, "Set value for option \'%s\' as %s.\n", name.c_str(), value.c_str());
        return true;
      } // end if
      else {
        reporter->printf(R_NL, R_BASIC, "Attempted to set value for option \'%s\' as string, but type is %s.  Ignoring request.\n", name.c_str(), option_list_[i]->type().c_str());
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyStringValue

////////////
// Option //
////////////

// Print
void Option::print(const Reporter* reporter) const
{

  // Print information
  reporter->printf(R_NL, R_BASIC, "Name        : %s\n", name_.c_str());
  reporter->printf(R_NL, R_BASIC, "Type        : %s\n", type_.c_str());
  if (type_.compare("bool") == 0) {
    reporter->printf(R_NL, R_BASIC, "Value       : %s\n", (value_bool_) ? "true" : "false");
  }
  else if (type_.compare("double") == 0) {
    reporter->printf(R_NL, R_BASIC, "Value       : %+e\n", value_double_);
    reporter->printf(R_NL, R_BASIC, "Lower bound : %+e\n", lower_bound_double_);
    reporter->printf(R_NL, R_BASIC, "Upper bound : %+e\n", upper_bound_double_);
  } // end if
  else if (type_.compare("integer") == 0) {
    reporter->printf(R_NL, R_BASIC, "Value       : %d\n", value_int_);
    reporter->printf(R_NL, R_BASIC, "Lower bound : %d\n", lower_bound_int_);
    reporter->printf(R_NL, R_BASIC, "Upper bound : %d\n", upper_bound_int_);
  } // end else if
  else {
    reporter->printf(R_NL, R_BASIC, "Value       : %s\n", value_string_.c_str());
  }
  reporter->printf(R_NL, R_BASIC, "Description : %s\n", description_.c_str());

} // end print

} // namespace NonOpt
