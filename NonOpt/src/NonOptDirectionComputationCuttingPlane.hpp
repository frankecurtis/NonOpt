// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTDIRECTIONCOMPUTATIONCUTTINGPLANE_HPP__
#define __NONOPTDIRECTIONCOMPUTATIONCUTTINGPLANE_HPP__

#include "NonOptDirectionComputation.hpp"

namespace NonOpt
{

/**
 * DirectionComputationCuttingPlane class
 */
class DirectionComputationCuttingPlane : public DirectionComputation
{

public:
  /** @name Constructor */
  //@{
  /**
   * Constructor
   */
  DirectionComputationCuttingPlane(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~DirectionComputationCuttingPlane(){};

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
  std::string name() { return "CuttingPlane"; };
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
  DirectionComputationCuttingPlane(const DirectionComputationCuttingPlane&);
  /**
   * Overloaded equals operator
   */
  void operator=(const DirectionComputationCuttingPlane&);
  //@}

  /** @name Private members */
  //@{
  bool add_far_points_;
  bool fail_on_iteration_limit_;
  bool fail_on_QP_failure_;
  bool try_aggregation_;
  bool try_gradient_step_;
  bool try_shortened_step_;
  double aggregation_size_threshold_;
  double downshift_constant_;
  double gradient_stepsize_;
  double shortened_stepsize_;
  double step_acceptance_tolerance_;
  int inner_iteration_limit_;
  int qp_small_limit_;
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

}; // end DirectionComputationCuttingPlane

} // namespace NonOpt

#endif /* __NONOPTDIRECTIONCOMPUTATIONCUTTINGPLANE_HPP__ */
