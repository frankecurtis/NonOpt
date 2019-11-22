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
      : increment_(1e-8),
        tolerance_(1e-4){};
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
  /**
   * Set tolerance
   */
  inline void setTolerance(double tolerance) { if (tolerance > 0.0) { tolerance_ = tolerance; } }
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
  double tolerance_;
  //@}

};  // end DerivativeChecker

}  // namespace NonOpt

#endif /* __NONOPTDERIVATIVECHECKER_HPP__ */
