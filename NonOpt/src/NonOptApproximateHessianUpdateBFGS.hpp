// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTAPPROXIMATEHESSIANUPDATEBFGS_HPP__
#define __NONOPTAPPROXIMATEHESSIANUPDATEBFGS_HPP__

#include "NonOptApproximateHessianUpdate.hpp"

namespace NonOpt
{

/**
 * ApproximateHessianUpdateBFGS class
 */
class ApproximateHessianUpdateBFGS : public ApproximateHessianUpdate
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  ApproximateHessianUpdateBFGS(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~ApproximateHessianUpdateBFGS(){};

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
   * Get iteration header values
   * \return string of header values
   */
  std::string iterationHeader() { return "Up. Fact. U?"; };
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues() { return "--------- --"; };
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "BFGS"; };
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
  void updateApproximateHessian(const Options* options,
                                Quantities* quantities,
                                const Reporter* reporter,
                                Strategies* strategies);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  ApproximateHessianUpdateBFGS(const ApproximateHessianUpdateBFGS&);
  /**
   * Overloaded equals operator
   */
  void operator=(const ApproximateHessianUpdateBFGS&);
  //@}

  /** @name Private members */
  //@{
  bool fail_on_tolerance_violation_;
  double correction_threshold_1_;
  double correction_threshold_2_;
  double norm_tolerance_;
  double product_tolerance_;
  //@}

  /** @name Private methods */
  //@{
  /**
   * Evaluate scalar for self-correcting BFGS update
   * \param[in,out] s is Vector representing iterate displacement
   * \param[in,out] y is Vector representing gradient displacement
   * \param[out] scalar is resulting scalar for update
   */
  void evaluateSelfCorrectingScalar(Vector& s,
                                    Vector& y,
                                    double& scalar);
  //@}

}; // end ApproximateHessianUpdateBFGS

} // namespace NonOpt

#endif /* __NONOPTAPPROXIMATEHESSIANUPDATEBFGS_HPP__ */
