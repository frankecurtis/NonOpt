// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTOPTIONS_HPP__
#define __NONOPTOPTIONS_HPP__

#include <string>
#include <vector>

#include "NonOptReporter.hpp"

namespace NonOpt
{

/**
 * Option class
 */
class Option
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor of bool option
   * \param[in] name is name of option
   * \param[in] type is type of option
   * \param[in] value is default value of option
   * \param[in] description is description of option
   */
  Option(std::string name,
         std::string type,
         bool value,
         std::string description)
    : description_(description),
      name_(name),
      type_(type),
      value_bool_(value){};
  /**
   * Constructor of double option
   * \param[in] name is name of option
   * \param[in] type is type of option
   * \param[in] value is default value of option
   * \param[in] lower_bound is lower bound for option values
   * \param[in] upper_bound is upper bound for option values
   * \param[in] description is description of option
   */
  Option(std::string name,
         std::string type,
         double value,
         double lower_bound,
         double upper_bound,
         std::string description)
    : description_(description),
      lower_bound_double_(lower_bound),
      name_(name),
      type_(type),
      upper_bound_double_(upper_bound),
      value_double_(value){};
  /**
   * Constructor of integer option
   * \param[in] name is name of option
   * \param[in] type is type of option
   * \param[in] value is default value of option
   * \param[in] lower_bound is lower bound for option values
   * \param[in] upper_bound is upper bound for option values
   * \param[in] description is description of option
   */
  Option(std::string name,
         std::string type,
         int value,
         int lower_bound,
         int upper_bound,
         std::string description)
    : description_(description),
      lower_bound_int_(lower_bound),
      name_(name),
      type_(type),
      upper_bound_int_(upper_bound),
      value_int_(value){};
  /**
   * Constructor of string option
   * \param[in] name is name of option
   * \param[in] type is type of option
   * \param[in] value is default value of option
   * \param[in] description is description of option
   */
  Option(std::string name,
         std::string type,
         std::string value,
         std::string description)
    : description_(description),
      name_(name),
      type_(type),
      value_string_(value){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~Option(){};
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print option information
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void print(const Reporter* reporter) const;
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get lower bound as a double
   * \return lower bound for option as double value
   */
  inline double lowerBoundAsDouble() const { return lower_bound_double_; };
  /**
   * Get lower bound as an integer
   * \return lower bound for option as integer value
   */
  inline int lowerBoundAsInteger() const { return lower_bound_int_; };
  /**
   * Get name
   * \return name of option as string
   */
  inline std::string name() const { return name_; };
  /**
   * Get type
   * \return type of option as string
   */
  inline std::string type() const { return type_; };
  /**
   * Get upper bound as a double
   * \return upper bound for option as double value
   */
  inline double upperBoundAsDouble() const { return upper_bound_double_; };
  /**
   * Get upper bound as an integer
   * \return upper bound for option as integer value
   */
  inline int upperBoundAsInteger() const { return upper_bound_int_; };
  /**
   * Get value as a bool
   * \return current value of option as bool
   */
  inline bool valueAsBool() const { return value_bool_; };
  /**
   * Get value as a double
   * \return current value of option as double value
   */
  inline double valueAsDouble() const { return value_double_; };
  /**
   * Get value as an integer
   * \return current value of option as integer value
   */
  inline int valueAsInteger() const { return value_int_; };
  /**
   * Get value as a string
   * \return current value of option as string
   */
  inline std::string valueAsString() const { return value_string_; };
  //@}

  /** @name Modify methods */
  //@{
  /**
   * Modify bool value
   * \param[in] value is new bool value of option to be set
   */
  inline void modifyBoolValue(bool value) { value_bool_ = value; };
  /**
   * Modify double value
   * \param[in] value is new double value of option to be set
   */
  inline void modifyDoubleValue(double value) { value_double_ = value; };
  /**
   * Modify integer value
   * \param[in] value is new integer value of option to be set
   */
  inline void modifyIntegerValue(int value) { value_int_ = value; };
  /**
   * Modify string value
   * \param[in] value is new string value of option to be set
   */
  inline void modifyStringValue(std::string value) { value_string_ = value; };
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Constructor (no arguments)
   */
  Option();
  /**
   * Copy constructor
   */
  Option(const Option&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Option&);
  //@}

  /**
   * Private members
   */
  //@{
  std::string description_;   /**< Description of option */
  double lower_bound_double_; /**< Lower bound for value as double */
  int lower_bound_int_;       /**< Lower bound for value as int */
  std::string name_;          /**< Name of option */
  std::string type_;          /**< Type of option */
  double upper_bound_double_; /**< Upper bound for value as double */
  int upper_bound_int_;       /**< Upper bound for value as int */
  bool value_bool_;           /**< Value as bool */
  double value_double_;       /**< Value as double */
  int value_int_;             /**< Value as int */
  std::string value_string_;  /**< Value as string */
  //@}

}; // end Option

/**
 * Options class
 */
class Options
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Options(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~Options(){};
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print options information
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void print(const Reporter* reporter) const;
  //@}

  /** @name Add methods */
  //@{
  /**
   * Add bool option
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is default value of option
   * \param[in] description is description of option
   */
  bool addBoolOption(const Reporter* reporter,
                     std::string name,
                     bool value,
                     std::string description);
  /**
   * Add double option
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is default value of option
   * \param[in] lower_bound is lower bound for option values
   * \param[in] upper_bound is upper bound for option values
   * \param[in] description is description of option
   */
  bool addDoubleOption(const Reporter* reporter,
                       std::string name,
                       double value,
                       double lower_bound,
                       double upper_bound,
                       std::string description);
  /**
   * Add integer option
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is default value of option
   * \param[in] lower_bound is lower bound for option values
   * \param[in] upper_bound is upper bound for option values
   * \param[in] description is description of option
   */
  bool addIntegerOption(const Reporter* reporter,
                        std::string name,
                        int value,
                        int lower_bound,
                        int upper_bound,
                        std::string description);
  /**
   * Add string option
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is default value of option
   * \param[in] description is description of option
   */
  bool addStringOption(const Reporter* reporter,
                       std::string name,
                       std::string value,
                       std::string description);
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get lower bound as a double
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is lower bound value of option as double
   */
  bool lowerBoundAsDouble(const Reporter* reporter,
                          std::string name,
                          double& value) const;
  /**
   * Get lower bound as an integer
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is lower bound value of option as integer
   */
  bool lowerBoundAsInteger(const Reporter* reporter,
                           std::string name,
                           int& value) const;
  /**
   * Get upper bound as a double
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is upper bound value of option as double
   */
  bool upperBoundAsDouble(const Reporter* reporter,
                          std::string name,
                          double& value) const;
  /**
   * Get upper bound as an integer
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is upper bound value of option as integer
   */
  bool upperBoundAsInteger(const Reporter* reporter,
                           std::string name,
                           int& value) const;
  /**
   * Get value as a bool
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is value of option as bool
   */
  bool valueAsBool(const Reporter* reporter,
                   std::string name,
                   bool& value) const;
  /**
   * Get value as a double
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is value of option as double
   */
  bool valueAsDouble(const Reporter* reporter,
                     std::string name,
                     double& value) const;
  /**
   * Get value as an integer
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is value of option as integer
   */
  bool valueAsInteger(const Reporter* reporter,
                      std::string name,
                      int& value) const;
  /**
   * Get value as a string
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[out] value is value of option as string
   */
  bool valueAsString(const Reporter* reporter,
                     std::string name,
                     std::string& value) const;
  //@}

  /** Modify methods */
  //@{
  /**
   * Modify from file
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] file_name is default file name as string
   */
  void modifyOptionsFromFile(const Reporter* reporter,
                             std::string file_name = "nonopt.opt");
  /**
   * Modify bool value
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is value of option to be set as bool
   */
  bool modifyBoolValue(const Reporter* reporter,
                       std::string name,
                       bool value);
  /**
   * Modify double value
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is value of option to be set as double
   */
  bool modifyDoubleValue(const Reporter* reporter,
                         std::string name,
                         double value);
  /**
   * Modify integer value
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is value of option to be set as integer
   */
  bool modifyIntegerValue(const Reporter* reporter,
                          std::string name,
                          int value);
  /**
   * Modify string value
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] name is name of option
   * \param[in] value is value of option to be set as string
   */
  bool modifyStringValue(const Reporter* reporter,
                         std::string name,
                         std::string value);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Options(const Options&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Options&);
  //@}

  /** @name Private members */
  //@{
  std::vector<std::shared_ptr<Option>> option_list_; /**< Vector of options */
  //@}

}; // end Options

} // namespace NonOpt

#endif /* __NONOPTOPTION_HPP__ */
