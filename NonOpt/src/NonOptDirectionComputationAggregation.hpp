// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTDirectionComputationAggregation_HPP__
#define __NONOPTDirectionComputationAggregation_HPP__

#include "NonOptDirectionComputation.hpp"
#include "NonOptRandomNumberGenerator.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class RandomNumberGenerator;

/**
 * DirectionComputationAggregation class
 */
class DirectionComputationAggregation : public DirectionComputation
{

 public:
  /** @name Constructor */
  //@{
  /**
   * Constructor
   */
  DirectionComputationAggregation(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~DirectionComputationAggregation(){};

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
  std::string name() { return "GradientCombinationAgg"; };
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
  DirectionComputationAggregation(const DirectionComputationAggregation&);
  /**
   * Overloaded equals operator
   */
  void operator=(const DirectionComputationAggregation&);
  //@}

  /** @name Private members */
  //@{
  bool fail_on_iteration_limit_;
  bool fail_on_QP_failure_;
  bool try_shortened_step_;
  bool used_full_last_;
  double downshift_constant_;
  double random_sample_fraction_;
  double shortened_stepsize_;
  double step_acceptance_tolerance_;
  double full_size_factor_;
  double * yk;
  double * gamm;
  bool do_agg_next_;
  int inner_iteration_limit_;
  RandomNumberGenerator random_number_generator_;
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

};  // end DirectionComputationAggregation

}  // namespace NonOpt

#endif /* __NONOPTDirectionComputationAggregation_HPP__ */
