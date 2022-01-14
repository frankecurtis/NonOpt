// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTDIRECTIONCOMPUTATIONGRADIENT_HPP__
#define __NONOPTDIRECTIONCOMPUTATIONGRADIENT_HPP__

#include "NonOptDirectionComputation.hpp"

namespace NonOpt
{

/**
 * DirectionComputationGradient class
 */
class DirectionComputationGradient : public DirectionComputation
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  DirectionComputationGradient(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~DirectionComputationGradient(){};

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
   * Get iteration header string
   * \return string of header values
   */
  std::string iterationHeader();
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues();
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "Gradient"; };
  //@}

  /** @name Direction computation method */
  //@{
  /**
   * Run direction computation
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  void computeDirection(const Options* options,
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
  DirectionComputationGradient(const DirectionComputationGradient&);
  /**
   * Overloaded equals operator
   */
  void operator=(const DirectionComputationGradient&);
  //@}

  /** @name Private members */
  //@{
  bool fail_on_QP_failure_;
  //@}

  /** @name Private methods */
  //@{
  /**
   * Converts QP solution to Quantity's direction
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  void convertQPSolutionToStep(Quantities* quantities,
                               Strategies* strategies);
  //@}

}; // end DirectionComputationGradient

} // namespace NonOpt

#endif /* __NONOPTDIRECTIONCOMPUTATIONGRADIENT_HPP__ */
