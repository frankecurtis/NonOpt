// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __AMPLPROBLEM_HPP__
#define __AMPLPROBLEM_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * AMPLProblem class
 */
class AMPLProblem : public Problem
{

 public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   * \param[in] stub is name of AMPL stub file
   */
  AMPLProblem(char* stub);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~AMPLProblem();
  //@}

  /** @name Get methods */
  //@{
  /**
   * Stub, i.e., name of problem
   * \return stub as character array
   */
  char* stub() { return stub_; };
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
  AMPLProblem();
  /**
   * Copy constructor
   */
  AMPLProblem(const AMPLProblem&);
  /**
   * Overloaded equals operator
   */
  void operator=(const AMPLProblem&);
  //@}

  /** @name Private members */
  //@{
  char* stub_; /**< Stub, i.e., name of problem */
  //@}

};  // end AMPLProblem

#endif /* __AMPLPROBLEM_HPP__ */
