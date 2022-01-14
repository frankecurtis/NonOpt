// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTTERMINATION_HPP__
#define __NONOPTTERMINATION_HPP__

#include "NonOptOptions.hpp"
#include "NonOptQuantities.hpp"
#include "NonOptReporter.hpp"
#include "NonOptStrategies.hpp"
#include "NonOptStrategy.hpp"

namespace NonOpt
{

/**
 * Forward declarations
 */
class Options;
class Quantities;
class Reporter;
class Strategy;
class Strategies;

/**
 * Termination class
 */
class Termination : public Strategy
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Termination(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~Termination(){};
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   */
  virtual void addOptions(Options* options) = 0;
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   */
  virtual void setOptions(Options* options) = 0;
  //@}

  /** @name Initialization method */
  //@{
  /**
   * Initialize strategy
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void initialize(const Options* options,
                          Quantities* quantities,
                          const Reporter* reporter) = 0;
  //@}

  /** @name Get methods */
  //@{
  /**
   * Iteration header string
   * \return string of header values
   */
  virtual std::string iterationHeader() = 0;
  /**
   * Iteration null values string
   * \return string of null values
   */
  virtual std::string iterationNullValues() = 0;
  /**
   * Name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  /**
   * Status
   * \return current status of termination
   */
  inline TE_Status status() { return status_; };
  /**
   * Termination indicator based on objective tolerance
   * \return bool with termination indicator
   */
  inline bool const terminateObjective() const { return terminate_objective_; };
  /**
   * Termination indicator based on objective similarity
   * \return bool with termination indicator
   */
  inline bool const terminateObjectiveSimilarity() const { return terminate_objective_similarity_; };
  /**
   * Termination indicator based on stationarity
   * \return bool with termination indicator
   */
  inline bool const terminateStationary() const { return terminate_stationary_; };
  /**
   * Update radii indicator
   * \return bool with update indicator
   */
  inline bool const updateRadii() const { return update_radii_; };
  /**
   * Update radii indicator for use in direction computation
   * \return bool with update indicator
   */
  inline bool const updateRadiiDirectionComputation() const { return update_radii_direction_computation_; };
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set status
   * \param[in] status is new status to be set
   */
  inline void setStatus(TE_Status status) { status_ = status; };
  //@}

  /** @name Check methods */
  //@{
  /**
   * Check all conditions
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  virtual void checkConditions(const Options* options,
                               Quantities* quantities,
                               const Reporter* reporter,
                               Strategies* strategies) = 0;
  /**
   * Check conditions for use in direction computation
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  virtual void checkConditionsDirectionComputation(const Options* options,
                                                   Quantities* quantities,
                                                   const Reporter* reporter,
                                                   Strategies* strategies) = 0;
  //@}

protected:
  /** @name Protected members */
  //@{
  bool terminate_objective_;                /**< Indicator for termination based on objective tolerance */
  bool terminate_objective_similarity_;     /**< Indicator for termination based on objective similarity */
  bool terminate_stationary_;               /**< Indicator for termination based on stationarity */
  bool update_radii_;                       /**< Indicator for radii update */
  bool update_radii_direction_computation_; /**< Indicator for radii update for use in direction computation */
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Termination(const Termination&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Termination&);
  //@}

  /** @name Private members */
  //@{
  TE_Status status_; /**< Termination status */
  //@}

}; // end Termination

} // namespace NonOpt

#endif /* __NONOPTTERMINATION_HPP__ */
