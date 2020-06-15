// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTAPPROXIMATEHESSIANUPDATE_HPP__
#define __NONOPTAPPROXIMATEHESSIANUPDATE_HPP__

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
 * ApproximateHessianUpdate class
 */
class ApproximateHessianUpdate : public Strategy
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  ApproximateHessianUpdate(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~ApproximateHessianUpdate(){};

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
  inline AH_Status status() { return status_; };
  //@}

  /** @name Set method */
  //@{
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(AH_Status status) { status_ = status; };
  //@}

  /** @name Update method */
  //@{
  /**
   * Run approximate Hessian update
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  virtual void updateApproximateHessian(const Options* options,
                                        Quantities* quantities,
                                        const Reporter* reporter,
                                        Strategies* strategies) = 0;
  //@}

protected:
  /** @name Protected members */
  //@{
  bool initial_update_performed_; /**< Indicates initial update has been performed */
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  ApproximateHessianUpdate(const ApproximateHessianUpdate&);
  /**
   * Overloaded equals operator
   */
  void operator=(const ApproximateHessianUpdate&);
  //@}

  /** @name Private members */
  //@{
  AH_Status status_; /**< Termination status */
  //@}

}; // set ApproximateHessianUpdate

} // namespace NonOpt

#endif /* __NONOPTAPPROXIMATEHESSIANUPDATE_HPP__ */
