// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTINVERSEHESSIANUPDATEBFGS_HPP__
#define __NONOPTINVERSEHESSIANUPDATEBFGS_HPP__

#include "NonOptInverseHessianUpdate.hpp"

namespace NonOpt
{

/**
 * InverseHessianUpdateBFGS class
 */
class InverseHessianUpdateBFGS : public InverseHessianUpdate
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  InverseHessianUpdateBFGS(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~InverseHessianUpdateBFGS(){};

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void addOptions(Options* options,
                  const Reporter* reporter);
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void setOptions(const Options* options,
                  const Reporter* reporter);
  //@}

  /** @name Initialize method */
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
  std::string iterationHeader() { return " Correction  U?"; };
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues() { return "-----------  --"; };
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "BFGS"; };
  //@}

  /** @name Inverse Hessian update method */
  //@{
  /**
   * Run inverse Hessian update
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  void updateInverseHessian(const Options* options,
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
  InverseHessianUpdateBFGS(const InverseHessianUpdateBFGS&);
  /**
   * Overloaded equals operator
   */
  void operator=(const InverseHessianUpdateBFGS&);
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

};  // end InverseHessianUpdateBFGS

}  // namespace NonOpt

#endif /* __NONOPTINVERSEHESSIANUPDATEBFGS_HPP__ */
