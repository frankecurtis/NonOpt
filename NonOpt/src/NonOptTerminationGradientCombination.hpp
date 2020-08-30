// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTSYMMETRICTERMINATIONGRADIENTCOMBINATION_HPP__
#define __NONOPTSYMMETRICTERMINATIONGRADIENTCOMBINATION_HPP__

#include "NonOptTermination.hpp"

namespace NonOpt
{

/**
  * TerminationGradientCombination class
  */
class TerminationGradientCombination : public Termination
{

public:
  /** @name Constructors */
  //@{
  /**
    * Construct TerminationGradientCombination
    */
  TerminationGradientCombination(){};
  //@}

  /** @name Destructor */
  //@{
  /**
    * Deconstruct
    */
  ~TerminationGradientCombination(){};
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
  std::string name() { return "GradientCombination"; };
  //@}

  /** @name Check condition methods */
  //@{
  /**
   * Check final termination conditions
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  bool checkFinal(const Options* options,
                  Quantities* quantities,
                  const Reporter* reporter,
                  Strategies* strategies) const;
  /**
   * Check radii not final
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  bool checkRadiiNotFinal(const Options* options,
                          Quantities* quantities,
                          const Reporter* reporter,
                          Strategies* strategies) const;
  /**
   * Check radii update
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  bool checkRadiiUpdate(const Options* options,
                        Quantities* quantities,
                        const Reporter* reporter,
                        Strategies* strategies) const;
  //@}

private:
  /** @name Default compiler generated methods
    * (Hidden to avoid implicit creation/calling.)
    */
  //@{
  /**
    * Copy constructor
    */
  TerminationGradientCombination(const TerminationGradientCombination&);
  /**
    * Overloaded equals operator
    */
  void operator=(const TerminationGradientCombination&);
  //@}

  /** @name Private members */
  //@{
  double stationarity_reference_;
  double stationarity_tolerance_;
  double stationarity_tolerance_factor_;
  //@}

}; // end TerminationGradientCombination

} // namespace NonOpt

#endif /* __NONOPTSYMMETRICMATRIXDENSE_HPP__ */
