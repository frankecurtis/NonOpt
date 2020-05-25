// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTPOINTSETUPDATEPROXIMITY_HPP__
#define __NONOPTPOINTSETUPDATEPROXIMITY_HPP__

#include "NonOptPointSetUpdate.hpp"

namespace NonOpt
{

/**
 * PointSetUpdateProximity class
 */
class PointSetUpdateProximity : public PointSetUpdate
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  PointSetUpdateProximity(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~PointSetUpdateProximity(){};

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
  std::string iterationHeader() { return ""; };
  /**
   * Get iteration null values string
   * \return string of null values
   */
  std::string iterationNullValues() { return ""; };
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "Proximity"; };
  //@}

  /** @name Point set update method */
  //@{
  /**
   * Run point set update
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in,out] strategies is pointer to Strategies object from NonOpt
   */
  void updatePointSet(const Options* options,
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
  PointSetUpdateProximity(const PointSetUpdateProximity&);
  /**
   * Overloaded equals operator
   */
  void operator=(const PointSetUpdateProximity&);
  //@}

  /** @name Private members */
  //@{
  double envelope_factor_;
  double size_factor_;
  //@}

}; // end PointSetUpdateProximity

} // namespace NonOpt

#endif /* __NONOPTPOINTSETUPDATEPROXIMITY_HPP__ */
