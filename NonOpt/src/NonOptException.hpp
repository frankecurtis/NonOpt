// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTEXCEPTION_HPP__
#define __NONOPTEXCEPTION_HPP__

#include <string>

#include "NonOptReporter.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Reporter;

/**
 * Exception class
 */
class Exception
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   * \param[in] message is exception identifier message
   * \param[in] file_name is file name in which exception occurred
   * \param[in] line_number is line number at which exception occurred
   * \param[in] type is type of exception
   */
  Exception(std::string message,
            std::string file_name,
            int line_number,
            std::string type = "Exception")
      : message_(message),
        file_name_(file_name),
        line_number_(line_number),
        type_(type) {}
  /**
   * Copy constructor
   * \param[in] copy is Exception object to copy
   */
  Exception(const Exception &copy)
      : message_(copy.message_),
        file_name_(copy.file_name_),
        line_number_(copy.line_number_),
        type_(copy.type_) {}
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~Exception() {}
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] type is type for which to report exception
   * \param[in] level is level at which to report exception
   */
  void print(const Reporter *reporter,
             ReportType type,
             ReportLevel level = R_BASIC) const
  {
    reporter->printf(type, level,
                     "\n"
                     "Exception of type \"%s\" in file \"%s\" at line %d\n"
                     "Exception message: %s\n",
                     type_.c_str(), file_name_.c_str(), line_number_,
                     message_.c_str());
  }
  //@}

  /** @name Get methods */
  //@{
  /**
   * Message
   * \return exception message as string
   */
  const std::string &message() const { return message_; };
  //@}

  /**
   * Private
   */
 private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Constructor w/ no arguments
   */
  Exception();
  /**
   * Overloaded equals operator
   */
  void operator=(const Exception &);
  //@}

  /** @name Private members */
  //@{
  std::string message_;   /**< Message of exception */
  std::string file_name_; /**< File name in which exception occurred */
  int line_number_;       /**< Line number at which exception occurred */
  std::string type_;      /**< Type of exception */
  //@}

};  // end Exception

#define THROW_EXCEPTION(__except_type, __message) \
  throw __except_type((__message), (__FILE__), (__LINE__));

#define ASSERT_EXCEPTION(__condition, __except_type, __message) \
  if (!(__condition)) {                                         \
    std::string new_message = #__condition;                     \
    new_message += " evaluated false: ";                        \
    new_message += __message;                                   \
    throw __except_type((new_message), (__FILE__), (__LINE__)); \
  }

#define DECLARE_EXCEPTION(__except_type)                                       \
  class __except_type : public Exception                                       \
  {                                                                            \
   public:                                                                     \
    __except_type(std::string message, std::string file_name, int line_number) \
        : Exception(message, file_name, line_number, #__except_type) {}        \
    __except_type(const __except_type &copy) : Exception(copy) {}              \
                                                                               \
   private:                                                                    \
    __except_type();                                                           \
    void operator=(const __except_type &);                                     \
  }

}  // namespace NonOpt

#endif /* __NONOPTEXCEPTION_HPP__ */
