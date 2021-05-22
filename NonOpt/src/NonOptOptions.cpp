// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <fstream>
#include <sstream>

#include "NonOptEnumerations.hpp"
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
  if (list_.size() == 0) {
    reporter->printf(R_NL, R_BASIC, "Option list is empty.\n");
    return;
  }

  // Print all options in list
  for (int i = 0; i < (int)list_.size(); i++) {
    list_[i]->print(reporter);
    if (i < (int)list_.size() - 1) {
      reporter->printf(R_NL, R_BASIC, "\n");
    }
  }

} // end print

// Options: Add bool option
bool Options::addBoolOption(std::string name,
                            bool value,
                            std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      message_ += "Option with name \'" + name + "\' already exists.  Ignoring addition request.\n";
      return false;
    } // end if
  }   // end for

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "bool", value, description));

  // Add to list
  list_.push_back(option);

  // Return true
  return true;

} // end addBoolOption

// Options: Add double option
bool Options::addDoubleOption(std::string name,
                              double value,
                              double lower_bound,
                              double upper_bound,
                              std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      message_ += "Option with name \'" + name + "\' already exists.  Ignoring addition request.\n";
      return false;
    } // end if
  }   // end for

  // Check that lower bound is lower than upper bound
  if (lower_bound > upper_bound) {
    message_ += "Attempted to add double option \'" + name + "\', but lower bound \'" + std::to_string(lower_bound) + "\' greater than upper bound \'" + std::to_string(upper_bound) + "\'.  Ignoring addition request.\n";
    return false;
  } // end if

  // Check that value is within bounds
  if (value < lower_bound || value > upper_bound) {
    message_ += "Attempted to add double option \'" + name + "\', but value \'" + std::to_string(value) + "\' outside of bound interval \'" + std::to_string(lower_bound) + "," + std::to_string(upper_bound) + "\'.  Ignoring addition request.\n";
    return false;
  } // end if

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "double", value, lower_bound, upper_bound, description));

  // Add to list
  list_.push_back(option);

  // Return true
  return true;

} // end addDoubleOption

// Options: Add integer option
bool Options::addIntegerOption(std::string name,
                               int value,
                               int lower_bound,
                               int upper_bound,
                               std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      message_ += "Option with name \'" + name + "\' already exists.  Ignoring addition request.\n";
      return false;
    } // end if
  }   // end for

  // Check that lower bound is lower than upper bound
  if (lower_bound > upper_bound) {
    message_ += "Attempted to add integer option \'" + name + "\', but lower bound \'" + std::to_string(lower_bound) + "\' greater than upper bound \'" + std::to_string(upper_bound) + "\'.  Ignoring addition request.\n";
    return false;
  } // end if

  // Check that value is within bounds
  if (value < lower_bound || value > upper_bound) {
    message_ += "Attempted to add integer option \'" + name + "\', but value \'" + std::to_string(value) + "\' outside of bound interval \'" + std::to_string(lower_bound) + "," + std::to_string(upper_bound) + "\'.  Ignoring addition request.\n";
    return false;
  } // end if

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "integer", value, lower_bound, upper_bound, description));

  // Add to list
  list_.push_back(option);

  // Return true
  return true;

} // end addIntegerOption

// Options: Add string option
bool Options::addStringOption(std::string name,
                              std::string value,
                              std::string description)
{

  // Check that option with given name doesn't already exist
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      message_ += "Option with name \'" + name + "\' already exists.  Ignoring addition request.\n";
      return false;
    } // end if
  }   // end for

  // Declare new pointer
  std::shared_ptr<Option> option(new Option(name, "string", value, description));

  // Add to list
  list_.push_back(option);

  // Return true
  return true;

} // end addStringOption

// Options: Get value as a bool
bool Options::valueAsBool(std::string name,
                          bool& value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("bool") == 0) {
        value = list_[i]->valueAsBool();
        return true;
      } // end if
      else {
        message_ += "Attempted to access value for option \'" + name + "\' as bool, but type is " + list_[i]->type() + ".  Returning false.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsBool

// Options: Get value as a double
bool Options::valueAsDouble(std::string name,
                            double& value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("double") == 0) {
        value = list_[i]->valueAsDouble();
        return true;
      } // end if
      else {
        message_ += "Attempted to access value for option \'" + name + "\' as double, but type is " + list_[i]->type() + ".  Returning false.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsDouble

// Options: Get value as an integer
bool Options::valueAsInteger(std::string name,
                             int& value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("integer") == 0) {
        value = list_[i]->valueAsInteger();
        return true;
      } // end if
      else {
        message_ += "Attempted to access value for option \'" + name + "\' as integer, but type is " + list_[i]->type() + ".  Returning false.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsInteger

// Options: Get value as a string
bool Options::valueAsString(std::string name,
                            std::string& value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("string") == 0) {
        value = list_[i]->valueAsString();
        return true;
      } // end if
      else {
        message_ += "Attempted to access value for option \'" + name + "\' as string, but type is " + list_[i]->type() + ".  Returning false.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false if option not found
  return false;

} // end valueAsString

// Options: Modify from file nonopt.opt
void Options::modifyOptionsFromFile(std::string file_name)
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
    for (int i = 0; i < (int)list_.size(); i++) {
      if (list_[i]->name().compare(s1) == 0) {
        name_found = true;
        if (list_[i]->type().compare("bool") == 0) {
          try {
            modifyBoolValue(s1, (s2.compare("true") == 0 || s2.compare("1") == 0));
          } catch (...) {
            message_ += "Attempted to set value for option \'" + s1 + "\', but cannot convert \'" + s2 + "\' to bool.  Ignoring request.\n";
          }
        }
        else if (list_[i]->type().compare("double") == 0) {
          try {
            modifyDoubleValue(s1, std::stod(s2));
          } catch (...) {
            message_ += "Attempted to set value for option \'" + s1 + "\', but cannot convert \'" + s2 + "\' to double.  Ignoring request.\n";
          }
        } // end if
        else if (list_[i]->type().compare("integer") == 0) {
          try {
            modifyIntegerValue(s1, (int)std::stod(s2));
          } catch (...) {
            message_ += "Attempted to set value for option \'" + s1 + "\', but cannot convert \'" + s2 + "\' to integer.  Ignoring request.\n";
          }
        } // end else if
        else {
          modifyStringValue(s1, s2);
        }
      } // end if
    }   // end for

    // Print message if name not found
    if (!name_found) {
      message_ += "Option with name \'" + s1 + "\' does not exist.  Ignoring request.\n";
    }

  } // end while

} // end modifyOptionsFromFile

// Options: Modify bool value
bool Options::modifyBoolValue(std::string name,
                              bool value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("bool") == 0) {
        list_[i]->modifyBoolValue(value);
        message_ += "Set value for option \'" + name + "\' as " + ((value) ? "true" : "false") + ".\n";
        return true;
      } // end if
      else {
        message_ += "Attempted to set value for option \'" + name + "\' as bool, but type is \'" + list_[i]->type() + "\'.  Ignoring request.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyBoolValue

// Options: Modify double value
bool Options::modifyDoubleValue(std::string name,
                                double value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("double") == 0) {
        if (value >= list_[i]->lowerBoundAsDouble() &&
            value <= list_[i]->upperBoundAsDouble()) {
          list_[i]->modifyDoubleValue(value);
          message_ += "Set value for option \'" + name + "\' as " + std::to_string(value) + ".\n";
          return true;
        } // end if
        else {
          message_ += "Attempted to set value for option \'" + name + "\', but value " + std::to_string(value) + " outside of bound interval " + std::to_string(list_[i]->lowerBoundAsDouble()) + "," + std::to_string(list_[i]->upperBoundAsDouble()) + ".  Ignoring request.\n";
          return false;
        } // end else
      }   // end if
      else {
        message_ += "Attempted to set value for option \'" + name + "\' as double, but type is \'" + list_[i]->type() + "\'.  Ignoring request.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyDoubleValue

// Options: Modify integer value
bool Options::modifyIntegerValue(std::string name,
                                 int value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("integer") == 0) {
        if (value >= list_[i]->lowerBoundAsInteger() &&
            value <= list_[i]->upperBoundAsInteger()) {
          list_[i]->modifyIntegerValue(value);
          message_ += "Set value for option \'" + name + "\' as " + std::to_string(value) + ".\n";
          return true;
        } // end if
        else {
          message_ += "Attempted to set value for option \'" + name + "\', but value " + std::to_string(value) + " outside of bound interval " + std::to_string(list_[i]->lowerBoundAsInteger()) + "," + std::to_string(list_[i]->upperBoundAsInteger()) + ".  Ignoring request.\n";
          return false;
        } // end else
      }   // end if
      else {
        message_ += "Attempted to set value for option \'" + name + "\' as integer, but type is \'" + list_[i]->type() + "\'.  Ignoring request.\n";
        return false;
      } // end else
    }   // end if
  }     // end for

  // Return false (not found)
  return false;

} // end modifyIntegerValue

// Options: Modify string value
bool Options::modifyStringValue(std::string name,
                                std::string value)
{

  // Loop through to find value
  for (int i = 0; i < (int)list_.size(); i++) {
    if (list_[i]->name().compare(name) == 0) {
      if (list_[i]->type().compare("string") == 0) {
        list_[i]->modifyStringValue(value);
        message_ += "Set value for option \'" + name + "\' as " + value + ".\n";
        return true;
      } // end if
      else {
        message_ += "Attempted to set value for option \'" + name + "\' as string, but type is \'" + list_[i]->type() + "\'.  Ignoring request.\n";
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
