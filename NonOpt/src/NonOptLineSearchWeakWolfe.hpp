// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTLINESEARCHWEAKWOLFE_HPP__
#define __NONOPTLINESEARCHWEAKWOLFE_HPP__

#include "NonOptLineSearch.hpp"

namespace NonOpt
{

/**
 * LineSearchWeakWolfe class
 */
class LineSearchWeakWolfe : public LineSearch
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  LineSearchWeakWolfe(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~LineSearchWeakWolfe(){};

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
  std::string iterationHeader() { return "  Stepsize "; };
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues() { return "-----------"; };
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "WeakWolfe"; };
  //@}

  /** @name Line search method */
  //@{
  /**
   * Run line search
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  void runLineSearch(const Options* options,
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
  LineSearchWeakWolfe(const LineSearchWeakWolfe&);
  /**
   * Overloaded equals operator
   */
  void operator=(const LineSearchWeakWolfe&);
  //@}

  /** @name Private members */
  //@{
  bool fail_on_small_interval_;
  double stepsize_initial_;
  double stepsize_minimum_;
  double stepsize_maximum_;
  double stepsize_sufficient_decrease_threshold_;
  double stepsize_curvature_threshold_;
  double stepsize_decrease_factor_;
  double stepsize_increase_factor_;
  double stepsize_bound_tolerance_;
  //@}

}; // end LineSearchWeakWolfe

} // namespace NonOpt

#endif /* __NONOPTLINESEARCHWEAKWOLFE_HPP__ */
