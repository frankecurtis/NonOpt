// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSTRATEGY_HPP__
#define __NONOPTSTRATEGY_HPP__

#include <string>

#include "NonOptOptions.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;

/**
 * Strategy class
 */
class Strategy
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Strategy(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~Strategy(){};

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void addOptions(Options* options,
                          const Reporter* reporter) = 0;
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void setOptions(const Options* options,
                          const Reporter* reporter) = 0;
  //@}

  /** @name Initialization method */
  //@{
  /**
   * Initialize strategy
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void initialize(const Options* options,
                          Quantities* quantities,
                          const Reporter* reporter) = 0;
  //@}

  /** @name Get method */
  //@{
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Strategy(const Strategy&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Strategy&);
  //@}

}; // end Strategy

} // namespace NonOpt

#endif /* __NONOPTSTRATEGY_HPP__ */
