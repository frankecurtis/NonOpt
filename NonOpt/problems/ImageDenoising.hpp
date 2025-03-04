// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __IMAGEDENOISING_HPP__
#define __IMAGEDENOISING_HPP__

#include "NonOptProblem.hpp"

using namespace NonOpt;

/**
 * ImageDenoising class
 */
class ImageDenoising : public Problem
{

public:
  /** @name Constructors */
  //@{
  /**
   * Constructor
   */
  ImageDenoising(int rows,
                 int cols,
                 double* image,
                 int regularizer,
                 double regularization,
                 double smoothing);
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~ImageDenoising();
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
  ImageDenoising(const ImageDenoising&);
  /**
   * Overloaded equals operator
   */
  void operator=(const ImageDenoising&);
  //@}

  /** @name Private members */
  //@{
  int number_of_variables_; /**< Number of variables */
  int rows_;                /**< Number of rows in image */
  int cols_;                /**< Number of cols in image */
  int regularizer_;         /**< Regularizer option */
  double* image_;           /**< Image array in grayscale */
  double regularization_;   /**< Regularization parameter */
  double smoothing_;        /**< Smoothing parameter */
  //@}

  /** @name Private methods */
  //@{
  double termFunction(const double t);
  double termDerivative(const double t);
  //@}

}; // end ImageDenoising

#endif /* __IMAGEDENOISING_HPP__ */
