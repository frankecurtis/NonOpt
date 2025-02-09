// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

// Description : Implementation for NonOpt of the objective
//                 f(x) = sum_{k=1..2*(n-2)}
//                          |y_l + sum_{h=1..3} h^2/l prod_{j=1..4} (x_{i+j}/|x_{i+j}|)*|x_{i+j}|^{j/(hl)}
//               where
//                 i = 2*div(k+3,4) - 2
//                 l = mod(k-1,4) + 1,
//                 y_1 = -14.4, y_2 = -6.8, y_3 = -4.2, and y_4 = -3.2
//               with initial point
//                 x_i = 0.8 for mod(i,4) == 0,
//                 x_i = -0.8 for mod(i,4) == 1,
//                 x_i = 1.2 for mod(i,4) == 2, and
//                 x_i = -1.2 for mod(i,4) == 3
// Notes       : THIS PROBLEM IS NONCONVEX

#ifndef __TEST29_13_HPP__
#define __TEST29_13_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * Test29_13 class
 */
class Test29_13 : public Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Test29_13(int n);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~Test29_13();
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
   * Evaluates objective and gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   * \param[out] g is the gradient value at "x", a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateObjectiveAndGradient(int n,
                                    const double* x,
                                    double& f,
                                    double* g);
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
  Test29_13(const Test29_13&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Test29_13&);
  //@}

  /** @name Private members */
  //@{
  int number_of_variables_; /**< Number of variables */
  //@}

}; // end Test29_13

#endif /* __Test29_13_HPP__ */
