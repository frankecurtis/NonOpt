// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTDERIVATIVECHECKER_HPP__
#define __NONOPTDERIVATIVECHECKER_HPP__

#include "NonOptPoint.hpp"
#include "NonOptReporter.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Point;
class Reporter;

/**
 * DerivativeChecker class
 */
class DerivativeChecker
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  DerivativeChecker()
      : increment_(1e-6){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~DerivativeChecker(){};
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set increment
   */
  inline void setIncrement(double increment) { if (increment > 0.0) { increment_ = increment; } }
  //@}

  /** @name Derivative checker methods */
  //@{
  /**
   * Check derivatives at a point
   * \param[in] reporter is pointer to Reporter
   * \param[in] point is pointer to Point
   */
  void checkDerivatives(const Reporter *reporter,
                        const std::shared_ptr<Point> point) const;
  //@}

 private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  DerivativeChecker(const DerivativeChecker &);
  /**
   * Overloaded equals operator
   */
  void operator=(const DerivativeChecker &);
  //@}

  /** @name Private members */
  //@{
  double increment_;
  //@}

};  // end DerivativeChecker

}  // namespace NonOpt

#endif /* __NONOPTDERIVATIVECHECKER_HPP__ */
