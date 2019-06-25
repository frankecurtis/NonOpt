// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

// Description : Implementation for NonOpt of the objective
//                 f(x) = sum_{i=1..n-1}
//                          max{x_i^4+x_{i+1}^2,
//                              (2-x_i)^2+(2-x_{i+1})^2,
//                              2*exp(-x_i+x_{i+1})}
//               with initial point
//                 x_i = 2.0 for all i = 1..n
// Notes       : THIS PROBLEM IS CONVEX
//               Optimal value: 2*(n-1)

#ifndef __CHAINEDCB3_1_HPP__
#define __CHAINEDCB3_1_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * ChainedCB3_1 class
 */
class ChainedCB3_1 : public Problem
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  ChainedCB3_1();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~ChainedCB3_1();
  //@}

  /** @name Get methods */
  //@{
  /**
   * Number of variables
   * \param[out] n is the number of variables, an integer (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool numberOfVariables(int& n);
  /**
   * Initial point
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[out] x is the initial point/iterate, a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool initialPoint(int n,
                    double* x);
  //@}

  /** @name Evaluate methods */
  //@{
  /**
   * Evaluates objective
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateObjective(int n,
                         const double* x,
                         double& f);
  /**
   * Evaluates gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] g is the gradient value at "x", a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateGradient(int n,
                        const double* x,
                        double* g);
  //@}

  /** @name Finalize methods */
  //@{
  /**
   * Finalizes solution
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is the final point/iterate, a constant double array
   * \param[in] f is the objective value at "x", a constant double
   * \param[in] g is the gradient value at "x", a constant double array
   * \return indicator of success (true) or failure (false)
   */
  bool finalizeSolution(int n,
                        const double* x,
                        double f,
                        const double* g);
  //@}

 private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  ChainedCB3_1(const ChainedCB3_1&);
  /**
   * Overloaded equals operator
   */
  void operator=(const ChainedCB3_1&);
  //@}

};  // end ChainedCB3_1

#endif /* __CHAINEDCB3_1_HPP__ */
