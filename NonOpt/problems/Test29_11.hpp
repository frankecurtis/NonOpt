// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

// Description : Implementation for NonOpt of the objective
//                 f(x) = sum_{k=1..2*(n-1)} |f_k(x)|
//               where
//                 i = (k+1)/2,
//                 f_k(x) = x_i + x_{i+1}*((1 + x_{i+1})*x_{i+1} - 14) - 29 for mod(k,2) = 0, and
//                 f_k(x) = x_i + x_{i+1}*((5 - x_{i+1})*x_{i+1} - 2) - 13 for mod(k,1)
//               with initial point
//                 x_i = 0.5 for all i = 1..n-1 and
//                 x_n = -2.0
// Notes       : THIS PROBLEM IS NONCONVEX

#ifndef __TEST29_11_HPP__
#define __TEST29_11_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * Test29_11 class
 */
class Test29_11 : public Problem
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Test29_11(int n);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~Test29_11();
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
  Test29_11(const Test29_11&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Test29_11&);
  //@}

  /** @name Private members */
  //@{
  int number_of_variables_; /**< Number of variables */
  //@}

};  // end Test29_11

#endif /* __Test29_11_HPP__ */
