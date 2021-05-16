// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTDERIVATIVECHECKER_HPP__
#define __NONOPTDERIVATIVECHECKER_HPP__

#include "NonOptEnumerations.hpp"
#include "NonOptOptions.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptStrategy.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;
class Strategy;

/**
 * DerivativeChecker class
 */
class DerivativeChecker : public Strategy
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  DerivativeChecker(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~DerivativeChecker(){};
  //@}

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

  /** @name Get methods */
  //@{
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  /**
   * Get status
   * \return current status of derivative checker
   */
  inline DE_Status status() { return status_; };
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(DE_Status status) { status_ = status; };
  //@}

  /** @name Check methods */
  //@{
  /**
   * Check derivatives at a point
   * \param[in] options is pointer to Options
   * \param[in] quantities is pointer to Quantities
   * \param[in] reporter is pointer to Reporter
   */
  virtual void checkDerivatives(const Options* options,
                                Quantities* quantities,
                                const Reporter* reporter) = 0;
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  DerivativeChecker(const DerivativeChecker&);
  /**
   * Overloaded equals operator
   */
  void operator=(const DerivativeChecker&);
  //@}

  /** @name Private members */
  //@{
  DE_Status status_; /**< DerivativeChecker status */
  //@}

}; // end DerivativeChecker

} // namespace NonOpt

#endif /* __NONOPTDERIVATIVECHECKER_HPP__ */
