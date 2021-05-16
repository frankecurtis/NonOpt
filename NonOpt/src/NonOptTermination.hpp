// Copyright (C) 2019 Frank E. Curtis
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
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void addOptions(Options* options,
                          const Reporter* reporter) = 0;
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  virtual void setOptions(const Options* options,
                          const Reporter* reporter) = 0;
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

  /** @name Get method */
  //@{
  /**
   * Get iteration header string
   * \return string of header values
   */
  virtual std::string iterationHeader() = 0;
  /**
   * Get iteration null values string
   * \return string of null values
   */
  virtual std::string iterationNullValues() = 0;
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  virtual std::string name() = 0;
  //@}

  /** @name Check condition methods */
  //@{
  /**
   * Check objective similarity
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  virtual bool checkObjectiveSimilarity(const Options* options,
                                        Quantities* quantities,
                                        const Reporter* reporter,
                                        Strategies* strategies) = 0;
  /**
   * Check radii final
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  virtual bool checkRadiiFinal(const Options* options,
                               Quantities* quantities,
                               const Reporter* reporter,
                               Strategies* strategies) const = 0;
  /**
   * Check radii update
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  virtual bool checkRadiiUpdate(const Options* options,
                                Quantities* quantities,
                                const Reporter* reporter,
                                Strategies* strategies) const = 0;
  /**
   * Check stationarity final conditions
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   * \param[in] strategies is pointer to Strategies object from NonOpt
   * \return bool to indicate conditions satisfied or not
   */
  virtual bool checkStationarityFinal(const Options* options,
                                      Quantities* quantities,
                                      const Reporter* reporter,
                                      Strategies* strategies) const = 0;

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

}; // end Termination

} // namespace NonOpt

#endif /* __NONOPTTERMINATION_HPP__ */
