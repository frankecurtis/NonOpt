// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSYMMETRICTERMINATIONBASIC_HPP__
#define __NONOPTSYMMETRICTERMINATIONBASIC_HPP__

#include "NonOptTermination.hpp"

namespace NonOpt
{

/**
 * TerminationBasic class
 */
class TerminationBasic : public Termination
{

public:
  /** @name Constructors */
  //@{
  /**
   * Construct TerminationBasic
   */
  TerminationBasic(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destruct
   */
  ~TerminationBasic(){};
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
   * Get iteration header values
   * \return string of header values
   */
  std::string iterationHeader() { return " |Grad.|  |G. Cmb.|"; };
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues() { return "--------- ---------"; };
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "Basic"; };
  //@}

  /** @name Check condition methods */
  //@{
  /**
   * Check all conditions
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  void checkConditions(const Options* options,
                       Quantities* quantities,
                       const Reporter* reporter,
                       Strategies* strategies);
  /**
   * Check conditions for use in direction computation
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  void checkConditionsDirectionComputation(const Options* options,
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
  TerminationBasic(const TerminationBasic&);
  /**
   * Overloaded equals operator
   */
  void operator=(const TerminationBasic&);
  //@}

  /** @name Private members */
  //@{
  int objective_similarity_counter_;
  int objective_similarity_limit_;
  double objective_reference_;
  double objective_similarity_tolerance_;
  double objective_tolerance_;
  double stationarity_reference_;
  double stationarity_tolerance_factor_;
  //@}

}; // end TerminationBasic

} // namespace NonOpt

#endif /* __NONOPTSYMMETRICTERMINATIONBASIC_HPP__ */
