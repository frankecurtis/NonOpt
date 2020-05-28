// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTPOINTSETUPDATE_HPP__
#define __NONOPTPOINTSETUPDATE_HPP__

#include <memory>
#include <string>

#include "NonOptEnumerations.hpp"
#include "NonOptOptions.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptStrategies.hpp"
#include "NonOptStrategy.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;
class Strategies;
class Strategy;

/**
 * PointSetUpdate class
 */
class PointSetUpdate : public Strategy
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  PointSetUpdate(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~PointSetUpdate(){};

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
   * Get iteration header string
   * \return string of header values
   */
  virtual std::string iterationHeader() = 0;
  /**
   * Get iteration null values string
   * \return string of null values
   */
  virtual std::string iterationNullValues() = 0;
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  /**
   * Get status
   * \return current status of strategy
   */
  inline PS_Status status() { return status_; };
  //@}

  /** @name Set method */
  //@{
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(PS_Status status) { status_ = status; };
  //@}

  /** @name Point set update method */
  //@{
  /**
   * Run point set update
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  virtual void updatePointSet(const Options* options,
                              Quantities* quantities,
                              const Reporter* reporter,
                              Strategies* strategies) = 0;
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  PointSetUpdate(const PointSetUpdate&);
  /**
   * Overloaded equals operator
   */
  void operator=(const PointSetUpdate&);
  //@}

  /** @name Private members */
  //@{
  PS_Status status_; /**< Termination status */
  //@}

}; // end PointSetUpdate

} // namespace NonOpt

#endif /* __NONOPTPOINTSETUPDATE_HPP__ */
