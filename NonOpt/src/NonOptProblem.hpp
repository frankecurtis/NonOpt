// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTPROBLEM_HPP__
#define __NONOPTPROBLEM_HPP__

#include <iostream>

namespace NonOpt
{

/**
 * Problem class
 */
class Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  Problem(){};
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  virtual ~Problem(){};

  /** @name Get methods */
  //@{
  /**
   * Returns number of variables
   * \param[out] n is the number of variables, an integer (return value)
   */
  virtual bool numberOfVariables(int& n) = 0;
  /**
   * Returns initial point
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[out] x is the initial point/iterate, a double array (return value)
   */
  virtual bool initialPoint(int n,
                            double* x) = 0;
  //@}

  /** @name Evaluate methods */
  //@{
  /**
   * Evaluates objective
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   */
  virtual bool evaluateObjective(int n,
                                 const double* x,
                                 double& f) = 0;
  /**
   * Evaluates objective and gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   * \param[out] g is the gradient value at "x", a double array (return value)
   */
  virtual bool evaluateObjectiveAndGradient(int n,
                                            const double* x,
                                            double& f,
                                            double* g)
  {
    ///////////////////////////////////////////////
    // Default method if not overwritten by user //
    ///////////////////////////////////////////////

    // Evaluate function
    bool objective_evaluation_success = evaluateObjective(n, x, f);

    // Evaluate gradient
    bool gradient_evaluation_success = evaluateGradient(n, x, g);

    // Return
    return (objective_evaluation_success && gradient_evaluation_success);
  }
  /**
   * Evaluates gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] g is the gradient value at "x", a double array (return value)
   */
  virtual bool evaluateGradient(int n,
                                const double* x,
                                double* g) = 0;
  //@}

  /** @name Finalize methods */
  //@{
  /**
   * Finalizes solution
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is the final point/iterate, a constant double array
   * \param[in] f is the objective value at "x", a constant double
   * \param[in] g is the gradient value at "x", a constant double array
   */
  virtual bool finalizeSolution(int n,
                                const double* x,
                                double f,
                                const double* g) = 0;
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  Problem(const Problem&);
  /**
   * Overloaded equals operator
   */
  void operator=(const Problem&);
  //@}

}; // end Problem

} // namespace NonOpt

#endif /* __NONOPTPROBLEM_HPP__ */
