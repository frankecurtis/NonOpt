// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSYMMETRICTERMINATIONSECONDQP_HPP__
#define __NONOPTSYMMETRICTERMINATIONSECONDQP_HPP__

#include "NonOptTermination.hpp"

namespace NonOpt
{

/**
  * TerminationSecondQP class
  */
class TerminationSecondQP : public Termination
{

public:
  /** @name Constructors */
  //@{
  /**
    * Construct TerminationSecondQP
    */
  TerminationSecondQP(){};
  //@}

  /** @name Destructor */
  //@{
  /**
    * Destruct
    */
  ~TerminationSecondQP(){};
  //@}

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
  std::string iterationHeader() { return " |Grad.|  |G. C. H| |G. C. I|"; };
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues() { return "--------- --------- ---------"; };
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "SecondQP"; };
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
  TerminationSecondQP(const TerminationSecondQP&);
  /**
    * Overloaded equals operator
    */
  void operator=(const TerminationSecondQP&);
  //@}

  /** @name Private members */
  //@{
  int objective_similarity_counter_;
  int objective_similarity_limit_;
  double objective_reference_;
  double objective_similarity_tolerance_;
  double stationarity_reference_;
  double stationarity_tolerance_factor_;
  //@}

  /** @name Private methods */
  //@{
  void solveQP(const Options* options,
               Quantities* quantities,
               const Reporter* reporter,
               Strategies* strategies);
  //@}

}; // end TerminationSecondQP

} // namespace NonOpt

#endif /* __NONOPTSYMMETRICTERMINATIONSECONDQP_HPP__ */
