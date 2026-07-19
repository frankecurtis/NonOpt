// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef NONOPTRANDOMNUMBERGENERATOR_HPP_
#define NONOPTRANDOMNUMBERGENERATOR_HPP_

#include <random>

namespace NonOpt
{

/**
 * RandomNumberGenerator class
 */
class RandomNumberGenerator
{

public:
  /** @name Constructors */
  //@{
  /**
   * Declare RandomNumberGenerator
   */
  RandomNumberGenerator() { generator.seed(0); };
  //@}

  /** @name Destructor */
  //@{
  /**
   * Delete
   */
  ~RandomNumberGenerator(){};
  //@}

  /** @name Set methods */
  //@{
  inline void resetSeed() { generator.seed(0); };
  //@}

  /** @name Generate methods */
  //@{
  /**
   * Generate standard normally distributed value
   * \return double value from standard normal distribution
   */
  double generateStandardNormal()
  {
    std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
  };
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  RandomNumberGenerator(const RandomNumberGenerator&) = delete;
  /**
   * Overloaded equals operator
   */
  void operator=(const RandomNumberGenerator&) = delete;
  //@}

  /** @name * Private members */
  //@{
  std::default_random_engine generator; /**< Random number generator */
  //@}

}; // end RandomNumberGenerator

} // namespace NonOpt

#endif /* NONOPTRANDOMNUMBERGENERATOR_HPP_ */
