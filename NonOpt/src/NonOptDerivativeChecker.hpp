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
  DerivativeChecker(double increment)
      : increment_(increment){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~DerivativeChecker();
  //@}

  /** @name Derivative checker methods */
  //@{
  /**
   * Check derivatives at a point
   * \param[in] point is pointer to Point
   */
  void checkDerivativesAtPoint(const Reporter *reporter,
                               const Point *point) const;
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
