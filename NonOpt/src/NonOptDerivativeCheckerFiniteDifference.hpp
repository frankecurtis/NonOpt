// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTDERIVATIVECHECKERFINITEDIFFERENCE_HPP__
#define __NONOPTDERIVATIVECHECKERFINITEDIFFERENCE_HPP__

#include "NonOptDerivativeChecker.hpp"

namespace NonOpt
{

/**
 * DerivativeCheckerFiniteDifference class
 */
class DerivativeCheckerFiniteDifference : public DerivativeChecker
{

public:
  /** @name Constructors */
  //@{
  /**
   * Construct DerivativeCheckerFiniteDifference
   */
  DerivativeCheckerFiniteDifference(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destruct
   */
  ~DerivativeCheckerFiniteDifference(){};
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   */
  void addOptions(Options* options);
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   */
  void setOptions(Options* options);
  //@}

  /** @name Initialization method */
  //@{
  /**
   * Initialize strategy
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void initialize(const Options* options,
                  Quantities* quantities,
                  const Reporter* reporter);
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "FiniteDifference"; };
  //@}

  /** @name Check methods */
  //@{
  /**
   * Check derivatives at a point
   * \param[in] options is pointer to Options
   * \param[in] quantities is pointer to Quantities
   * \param[in] reporter is pointer to Reporter
   */
  void checkDerivatives(const Options* options,
                        Quantities* quantities,
                        const Reporter* reporter);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  DerivativeCheckerFiniteDifference(const DerivativeCheckerFiniteDifference&);
  /**
   * Overloaded equals operator
   */
  void operator=(const DerivativeCheckerFiniteDifference&);
  //@}

  /** @name Private members */
  //@{
  bool check_derivatives_;
  double increment_;
  double tolerance_;
  //@}

}; // end DerivativeCheckerFiniteDifference

} // namespace NonOpt

#endif /* __NONOPTDERIVATIVECHECKERFINITEDIFFERENCE_HPP__ */
