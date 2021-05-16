// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

// Description : Implementation for NonOpt of the objective
//                 f(x) = sum_{i=1..n-1}
//                          (|x_i|^{x_{i+1}^2+1} + |x_{i+1}|^{x_i^2+1})
//               with initial point
//                 x_i = -1.0 for all i = 1..n when mod(i,2) == 1
//                 x_i =  1.0 for all i = 1..n when mod(i,2) == 0
// Notes       : THIS PROBLEM IS NONCONVEX
//               Optimal value: 0.0

#ifndef __BROWNFUNCTION_2_HPP__
#define __BROWNFUNCTION_2_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * BrownFunction_2 class
 */
class BrownFunction_2 : public Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  BrownFunction_2(int n);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~BrownFunction_2();
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
  BrownFunction_2(const BrownFunction_2&);
  /**
   * Overloaded equals operator
   */
  void operator=(const BrownFunction_2&);
  //@}

  /** @name Private members */
  //@{
  int number_of_variables_; /**< Number of variables */
  //@}

}; // end BrownFunction_2

#endif /* __BROWNFUNCTION_2_HPP__ */
