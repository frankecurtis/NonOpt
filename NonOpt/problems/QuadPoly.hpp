// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

// Description : Implementation for NonOpt of the objective
//                 f(x) = c'*x + 0.5*x'*Q*x + max(b + A*x)
//               with initial point x = 0.0
// Notes       : THIS PROBLEM IS CONVEX
//               Optimal value: 0.0

#ifndef __QUADPOLY_HPP__
#define __QUADPOLY_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * QuadPoly class
 */
class QuadPoly : public Problem
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   * \param[in] n is the number of variables
   * \param[in] m is the number of affine functions
   * \param[in] a is the number of active affine functions at the solution
   * \param[in] s is the scaling factor for the quadratic term
   */
  QuadPoly(int n,
           int m,
           int a,
           double s);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~QuadPoly();
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
   * Constructor (no arguments)
   */
  QuadPoly();
  /**
   * Copy constructor
   */
  QuadPoly(const QuadPoly&);
  /**
   * Overloaded equals operator
   */
  void operator=(const QuadPoly&);
  //@}

  /** @name Private members */
  //@{
  int number_of_active_affine_; /**< Number of active affine functions at solution */
  int number_of_affine_;        /**< Number of affine functions in max term */
  int number_of_variables_;     /**< Number of variables */
  double* constant_;            /**< Constant term in max, i.e., "b" */
  double* linear_;              /**< Linear term, i.e., "c" */
  double* matrix_;              /**< Linear term in max, i.e., "A" */
  double* symmetric_matrix_;    /**< Quadratic term, i.e., "Q" */
  //@}

};  // end QuadPoly

#endif /* __QUADPOLY_HPP__ */
