// Copyright (C) 2024 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Lara Zebiane

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <iomanip>

#include "NonOptBLASLAPACK.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptQPSolverInteriorPoint.hpp"

namespace NonOpt
{

// Constructor
QPSolverInteriorPoint::QPSolverInteriorPoint()
  : inexact_solution_tolerance_(0.0),
    scalar_(0.0) {}

// Destructor
QPSolverInteriorPoint::~QPSolverInteriorPoint()
{

  // Delete arrays

} // end destructor

// Add options
void QPSolverInteriorPoint::addOptions(Options* options)
{

  // Add bool options
  options->addBoolOption("QPIPM_allow_inexact_termination",
                         false,
                         "Indicator for whether to allow early termination.\n"
                         "Default     : false.");
  // Add double options
  options->addDoubleOption("QPIPM_barrier_parameter_factor",
                           1.00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for setting barrier parameter in predictor-corrector method.\n"
                           "Default     : 1.00.");
  options->addDoubleOption("QPIPM_barrier_parameter_initial",
                           5e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Initial barrier parameter.\n"
                           "Default     : 5e-01.");
  options->addDoubleOption("QPIPM_barrier_parameter_maximum",
                           1e+04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Maximum barrier parameter.\n"
                           "Default     : 1e+02.");
  options->addDoubleOption("QPIPM_barrier_parameter_minimum",
                           1e-12,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Minimum barrier parameter.\n"
                           "Default     : 1e-10.");
  options->addDoubleOption("QPIPM_fraction_to_boundary_tolerance",
                           1e-03,
                           0.0,
                           1.0,
                           "Fraction-to-boundary parameter for line searches.\n"
                           "Default     : 1e-03.");
  options->addDoubleOption("QPIPM_kkt_tolerance",
                           1e-08,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining optimality.  If the KKT error\n"
                           "              computed by the algorithm falls below this tolerance, then\n"
                           "              the algorithm terminates with a message of success.\n"
                           "Default     : 1e-08.");
  options->addDoubleOption("QPIPM_inexact_termination_descent_tolerance",
                           1e-04,
                           0.0,
                           1.0,
                           "Descent direction tolerance for inexactness conditions.\n"
                           "Default     : 1e-04.");
  options->addDoubleOption("QPIPM_inexact_termination_initialization_factor",
                           0.0,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Determines number of iterations (fraction of gradient list size)\n"
                           "              to perform before inexact termination conditions are checked.\n"
                           "Default     : 2.5e-01.");
  options->addDoubleOption("QPIPM_inexact_termination_ratio_minimum",
                           1e-02,
                           0.0,
                           1.0,
                           "Minimum value for ratio used in inexact termination condition.\n"
                           "Default     : 1e-02.");
  options->addDoubleOption("QPIPM_solution_initialization_factor",
                           1e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Factor for initializing interior solution vectors.\n"
                           "Default     : 1e-01.");
  // Add integer options
  options->addIntegerOption("QPIPM_inexact_termination_check_interval",
                            4,
                            1,
                            NONOPT_INT_INFINITY,
                            "Number of iterations to perform between checks of inexact\n"
                            "              termination conditions.\n"
                            "Default     : 4.");
  options->addIntegerOption("QPIPM_iteration_limit",
                            2e+02,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of iterations.\n"
                            "Default     : 1e+03.");

} // end addOptions

// Set options
void QPSolverInteriorPoint::setOptions(Options* options)
{

  // Read bool options
  options->valueAsBool("QPIPM_allow_inexact_termination", allow_inexact_termination_);
  // Read double options
  options->valueAsDouble("QPIPM_barrier_parameter_factor", barrier_parameter_factor_);
  options->valueAsDouble("QPIPM_barrier_parameter_initial", barrier_parameter_initial_);
  options->valueAsDouble("QPIPM_barrier_parameter_maximum", barrier_parameter_maximum_);
  options->valueAsDouble("QPIPM_barrier_parameter_minimum", barrier_parameter_minimum_);
  options->valueAsDouble("QPIPM_fraction_to_boundary_tolerance", fraction_to_boundary_tolerance_);
  options->valueAsDouble("QPIPM_kkt_tolerance", kkt_tolerance_);
  options->valueAsDouble("QPIPM_inexact_termination_descent_tolerance", inexact_termination_descent_tolerance_);
  options->valueAsDouble("QPIPM_inexact_termination_initialization_factor", inexact_termination_initialization_factor_);
  options->valueAsDouble("QPIPM_inexact_termination_ratio_minimum", inexact_termination_ratio_minimum_);
  options->valueAsDouble("QPIPM_solution_initialization_factor", solution_initialization_factor_);
  // Read integer options
  options->valueAsInteger("QPIPM_inexact_termination_check_interval", inexact_termination_check_interval_);
  options->valueAsInteger("QPIPM_iteration_limit", iteration_limit_);

} // end setOptions

// Initialize
void QPSolverInteriorPoint::initialize(const Options* options,
                                       Quantities* quantities,
                                       const Reporter* reporter)
{
  initializeData(quantities->numberOfVariables());
}

// Initialize data
void QPSolverInteriorPoint::initializeData(int gamma_length)
{

  // Set length parameters
  gamma_length_ = gamma_length;

  // Initialize problem data
  matrix_ = nullptr;
  vector_list_.clear();
  vector_.clear();
  scalar_ = -1.0;
  
  // Set status
  setStatus(QP_UNSET);

  // Set null solution
  setNullSolution();

} // end initializeData

// Get combination norm
double QPSolverInteriorPoint::combinationNormInf()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, combination_.values(), &increment);

  // Set norm
  double combination_norm_inf = fabs(combination_.values()[i - 1]);

  // Set maximum
  if (std::isnan(combination_norm_inf) || combination_norm_inf > NONOPT_DOUBLE_INFINITY) {
    combination_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return combination_norm_inf;

} // end combinationNormInf

// Get translated combination norm
double QPSolverInteriorPoint::combinationTranslatedNormInf()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, combination_translated_.values(), &increment);

  // Set norm
  double combination_translated_norm_inf = fabs(combination_translated_.values()[i - 1]);

  // Set maximum
  if (std::isnan(combination_translated_norm_inf) || combination_translated_norm_inf > NONOPT_DOUBLE_INFINITY) {
    combination_translated_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return combination_translated_norm_inf;

} // end combinationTranslatedNormInf

// Get translated combination 2-norm square
double QPSolverInteriorPoint::combinationTranslatedNorm2Squared()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Set norm
  double combination_translated_norm_2 = dnrm2_(&length, combination_translated_.values(), &increment);

  // Set maximum
  if (std::isnan(combination_translated_norm_2) || combination_translated_norm_2 > NONOPT_DOUBLE_INFINITY) {
    combination_translated_norm_2 = NONOPT_DOUBLE_INFINITY;
  }

  // Return 2-norm square
  return combination_translated_norm_2 * combination_translated_norm_2;

} // end combinationTranslatedNorm2Squared

// Get dual objective quadratic value
double QPSolverInteriorPoint::dualObjectiveQuadraticValue()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Set quadratic value
  double dual_objective_quadratic_value = -ddot_(&length, primal_solution_.values(), &increment, combination_translated_.values(), &increment);

  // Set maximum
  if (std::isnan(dual_objective_quadratic_value) || dual_objective_quadratic_value > NONOPT_DOUBLE_INFINITY) {
    dual_objective_quadratic_value = NONOPT_DOUBLE_INFINITY;
  }

  // Return
  return dual_objective_quadratic_value;

} // end dualObjectiveQuadraticValue

// Get objective quadratic value for feasible dual step
double QPSolverInteriorPoint::dualObjectiveQuadraticValueScaled()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Set quadratic value
  double dual_objective_quadratic_value_scaled = -ddot_(&length, primal_solution_.values(), &increment, combination_translated_.values(), &increment);

  // Scale so value corresponds to feasible dual step
  dual_objective_quadratic_value_scaled *= pow(primal_solution_projection_scalar_, 2.0);

  // Set maximum
  if (std::isnan(dual_objective_quadratic_value_scaled) || dual_objective_quadratic_value_scaled > NONOPT_DOUBLE_INFINITY) {
    dual_objective_quadratic_value_scaled = NONOPT_DOUBLE_INFINITY;
  }

  // Return
  return dual_objective_quadratic_value_scaled;

} // end dualObjectiveQuadraticValueScaled

// Get dual solution
void QPSolverInteriorPoint::dualSolution(double omega[], double gamma[])
{

  // Set inputs for BLASLAPACK
  int length = (int)vector_.size();
  int increment = 1;

  // Copy values
  dcopy_(&length, omega_.values(), &increment, omega, &increment);

  // Set inputs for BLASLAPACK
  length = gamma_length_;

  // Copy values
  dcopy_(&length, gamma_.values(), &increment, gamma, &increment);

} // end dualSolution

// Get dual solution, omega part
void QPSolverInteriorPoint::dualSolutionOmega(double omega[])
{

  // Set inputs for BLASLAPACK
  int length = (int)vector_.size();
  int increment = 1;

  // Copy values
  dcopy_(&length, omega_.values(), &increment, omega, &increment);

} // end dualSolutionOmega

// Get primal solution
void QPSolverInteriorPoint::primalSolution(double d[])
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, primal_solution_.values(), &increment, d, &increment);

} // end primalSolution

// Get feasible primal solution
void QPSolverInteriorPoint::primalSolutionFeasible(double d_feasible[])
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, primal_solution_feasible_.values(), &increment, d_feasible, &increment);

} // end primalSolutionFeasible

// Get primal solution inf-norm
double QPSolverInteriorPoint::primalSolutionNormInf()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, primal_solution_.values(), &increment);

  // Set norm
  double primal_solution_norm_inf = fabs(primal_solution_.values()[i - 1]);

  // Set maximum
  if (std::isnan(primal_solution_norm_inf) || primal_solution_norm_inf > NONOPT_DOUBLE_INFINITY) {
    primal_solution_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return primal_solution_norm_inf;

} // end primalSolutionNormInf

// Get primal solution 2-norm square
double QPSolverInteriorPoint::primalSolutionNorm2Squared()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Set norm
  double primal_solution_norm_2 = dnrm2_(&length, primal_solution_.values(), &increment);

  // Set maximum
  if (std::isnan(primal_solution_norm_2) || primal_solution_norm_2 > NONOPT_DOUBLE_INFINITY) {
    primal_solution_norm_2 = NONOPT_DOUBLE_INFINITY;
  }

  // Return 2-norm square
  return primal_solution_norm_2 * primal_solution_norm_2;

} // end primalSolutionNorm2Squared

// Get feasible primal solution norm
double QPSolverInteriorPoint::primalSolutionFeasibleNormInf()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, primal_solution_feasible_.values(), &increment);

  // Set norm
  double primal_solution_feasible_norm_inf = fabs(primal_solution_feasible_.values()[i - 1]);

  // Set maximum
  if (std::isnan(primal_solution_feasible_norm_inf) || primal_solution_feasible_norm_inf > NONOPT_DOUBLE_INFINITY) {
    primal_solution_feasible_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return primal_solution_feasible_norm_inf;

} // end primalSolutionFeasibleNormInf

// Initialize data
void QPSolverInteriorPoint::setNullSolution()
{

  // Algorithm parameters
  iteration_count_ = 0;
  kkt_error_ = -NONOPT_DOUBLE_INFINITY;
  dual_objective_reference_ = -NONOPT_DOUBLE_INFINITY;
  primal_directional_derivative_feasible_best_ = NONOPT_DOUBLE_INFINITY;
  primal_objective_feasible_best_ = NONOPT_DOUBLE_INFINITY;
  primal_objective_reference_ = NONOPT_DOUBLE_INFINITY;
  primal_objective_simple_ = NONOPT_DOUBLE_INFINITY;
  primal_quadratic_feasible_best_ = -NONOPT_DOUBLE_INFINITY;

  // Initialize dual values
  combination_.setLength(gamma_length_);
  combination_.scale(0.0);
  combination_translated_.setLength(gamma_length_);
  combination_translated_.scale(0.0);
  primal_solution_.setLength(gamma_length_);
  primal_solution_.scale(0.0);
  primal_solution_feasible_.setLength(gamma_length_);
  primal_solution_feasible_.scale(0.0);
  primal_solution_feasible_best_.setLength(gamma_length_);
  primal_solution_feasible_best_.scale(0.0);
  primal_solution_simple_.setLength(gamma_length_);
  primal_solution_simple_.scale(0.0);
  primal_solution_projection_scalar_ = 0.0;

} // end setNullSolution

// Add vectors
void QPSolverInteriorPoint::addData(const std::vector<std::shared_ptr<Vector>> vector_list,
                                    const std::vector<double> vector)
{

  // Loop through new elements
  for (int i = 0; i < (int)vector.size(); i++) {
    vector_list_.push_back(vector_list[i]);
    vector_.push_back(vector[i]);
  }

} // end addData

// Inexact termination condition
bool QPSolverInteriorPoint::inexactTerminationCondition(const Quantities* quantities,
                                                        const Reporter* reporter)
{

  /////////////////////////////////////////////////////////////////////////////
  // EVALUATE DUAL OBJECTIVE CORRESPONDING TO CURRENT DUAL SOLUTION ESTIMATE //
  /////////////////////////////////////////////////////////////////////////////

  // Initialize dual objective
  double dual_objective = -0.5 * dualObjectiveQuadraticValue();

  // Add linear term
  for (int i = 0; i < (int)vector_.size(); i++) {
    dual_objective += vector_[i] * omega_.values()[i];
  }

  // Add norm term
  if (scalar_ < NONOPT_DOUBLE_INFINITY) {
    for (int i = 0; i < gamma_length_; i++) {
      dual_objective -= scalar_ * gamma_.values()[i];
    }
  } // end if

  //////////////////////////////////////////////////////////////////////////////////
  // EVALUATE PRIMAL OBJECTIVE CORRESPONDING TO PRIMAL FEASIBLE SOLUTION ESTIMATE //
  //////////////////////////////////////////////////////////////////////////////////

  // Initialize primal objective
  double primal_objective_feasible = 0.5 * dualObjectiveQuadraticValueScaled();

  // Add max of linear terms
  double bPlusGd_max = -NONOPT_DOUBLE_INFINITY;
  for (int i = 0; i < (int)vector_.size(); i++) {
    double bgd = vector_[i] + vector_list_[i]->innerProduct(primal_solution_feasible_);
    if (bgd > bPlusGd_max) {
      bPlusGd_max = bgd;
    }
  } // end for
  primal_objective_feasible += bPlusGd_max;

  ////////////////////////////////////////////////////////////////////////////////
  // EVALUATE PRIMAL OBJECTIVE CORRESPONDING TO SIMPLE PRIMAL FEASIBLE SOLUTION //
  ////////////////////////////////////////////////////////////////////////////////

  // Check iteration count
  if (iteration_count_ == 1) {

    // Set simple primal feasible solution based on average of gradients
    for (int i = 0; i < (int)vector_list_.size(); i++) {
      primal_solution_simple_.addScaledVector(1.0 / (double)vector_list_.size(), *vector_list_[i].get());
    }
    double projection_scalar_ = fmin(1.0, scalar_ / primal_solution_simple_.normInf());
    primal_solution_simple_.scale(projection_scalar_);

    // Initialize primal objective
    primal_objective_simple_ = 0.5 * matrix_->innerProduct(primal_solution_simple_);

    // Add max of linear terms
    primal_objective_reference_ = -NONOPT_DOUBLE_INFINITY;
    double bPlusGd_max = -NONOPT_DOUBLE_INFINITY;
    for (int i = 0; i < (int)vector_.size(); i++) {
      if (vector_[i] > primal_objective_reference_) {
        primal_objective_reference_ = vector_[i];
      }
      double bgd = vector_[i] + vector_list_[i]->innerProduct(primal_solution_simple_);
      if (bgd > bPlusGd_max) {
        bPlusGd_max = bgd;
      }
    } // end for
    primal_objective_simple_ += bPlusGd_max;

    // Check for reversion to zero solution
    if (primal_objective_simple_ > primal_objective_reference_) {
      primal_solution_simple_.scale(0.0);
      primal_objective_simple_ = primal_objective_reference_;
    } // end if

  } // end if

  ///////////////////////////////////////////////////
  // UPDATE BEST PRIMAL FEASIBLE SOLUTION ESTIMATE //
  ///////////////////////////////////////////////////

  // Check iteration count
  if (iteration_count_ == 1) {
    dual_objective_reference_ = dual_objective;
    if (primal_objective_feasible <= primal_objective_simple_) {
      primal_objective_feasible_best_ = primal_objective_feasible;
      primal_solution_feasible_best_.copy(primal_solution_feasible_);
      primal_directional_derivative_feasible_best_ = vector_list_[0]->innerProduct(primal_solution_feasible_best_);
      primal_quadratic_feasible_best_ = matrix_->innerProduct(primal_solution_feasible_best_);
    } // end if
    else {
      primal_objective_feasible_best_ = primal_objective_simple_;
      primal_solution_feasible_best_.copy(primal_solution_simple_);
      primal_directional_derivative_feasible_best_ = vector_list_[0]->innerProduct(primal_solution_feasible_best_);
      primal_quadratic_feasible_best_ = matrix_->innerProduct(primal_solution_feasible_best_);
    } // end else
  }   // end if
  else {
    if (primal_objective_feasible < primal_objective_feasible_best_) {
      primal_objective_feasible_best_ = primal_objective_feasible;
      primal_solution_feasible_best_.copy(primal_solution_feasible_);
      primal_directional_derivative_feasible_best_ = vector_list_[0]->innerProduct(primal_solution_feasible_best_);
      primal_quadratic_feasible_best_ = matrix_->innerProduct(primal_solution_feasible_best_);
    } // end if
  }   // end else

  //////////////////////
  // CHECK CONDITIONS //
  //////////////////////

  // Set objective ratio
  double objective_ratio = NONOPT_DOUBLE_INFINITY;
  if (primal_objective_feasible_best_ < primal_objective_reference_) {
    objective_ratio = (primal_objective_reference_ - dual_objective_reference_) / (primal_objective_reference_ - primal_objective_feasible_best_);
  }

  // Set inexact termination ratio value
  double inexact_termination_ratio = 1.0 - (pow(quantities->inexactTerminationFactor(), 2.0) + 2 * quantities->inexactTerminationFactor()) / (objective_ratio - 1.0);

  // Initialize return value
  bool condition_bool = false;

  // Check "zero d" condition
  if (primalSolutionNormInf() <= inexact_solution_tolerance_ &&
      combinationNormInf() <= inexact_solution_tolerance_ &&
      combinationTranslatedNormInf() <= inexact_solution_tolerance_) {
    condition_bool = true;
  }
  // Check first "nonzero d" condition
  else if (vector_list_[0]->innerProduct(primal_solution_feasible_) <= -inexact_termination_descent_tolerance_ * dualObjectiveQuadraticValue()) {
    if ((pow(quantities->inexactTerminationFactor(), 2.0) + 2 * quantities->inexactTerminationFactor()) * (primal_objective_reference_ - primal_objective_feasible_best_) >=
        primal_objective_feasible_best_ - dual_objective) {
      primal_solution_.copy(primal_solution_feasible_best_);
      condition_bool = true;
    } // end if
    // Check second "nonzero d" condition
    else if (dual_objective - dual_objective_reference_ >=
             fmax(inexact_termination_ratio, inexact_termination_ratio_minimum_) * (primal_objective_feasible_best_ - dual_objective_reference_)) {
      primal_solution_.copy(primal_solution_feasible_best_);
      condition_bool = true;
    } // end if
  }

  // Return
  return condition_bool;

} // end inexactTerminationCondition

// Solve
void QPSolverInteriorPoint::solveQP(const Options* options,
                                    const Reporter* reporter,
                                    Quantities* quantities)
{

  // Initialize values
  setStatus(QP_UNSET);
  setNullSolution();
  iteration_count_ = 0;
  kkt_error_ = -NONOPT_DOUBLE_INFINITY;

  // try to solve QP, terminate on any exception
  try {

    // Check quantity compatibility
    if (!checkQuantityCompatibility()) {
      THROW_EXCEPTION(QP_INPUT_ERROR_EXCEPTION, "QP solve unsuccessful. Input error.");
    }

    // Print message
    reporter->printf(R_QP, R_PER_ITERATION, "\n");
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "Entering main iteration loop\n");

    // Set sizes
    int omega_length = (int)vector_.size();
    int total_length;
    if (scalar_ < NONOPT_DOUBLE_INFINITY) {
      total_length = 2 * gamma_length_ + omega_length;
    }
    else {
      total_length = omega_length;
    }

    // Print sizes
    reporter->printf(R_QP, R_PER_ITERATION, "========================\n"
                                            " Problem size quantities\n"
                                            "========================\n");
    reporter->printf(R_QP, R_PER_ITERATION, "n (gamma length) : %8d\n", gamma_length_);
    reporter->printf(R_QP, R_PER_ITERATION, "m (omega length) : %8d\n", omega_length);
    reporter->printf(R_QP, R_PER_ITERATION, "KKT matrix size  : %8d\n", total_length);

    // Construct A transpose
    At_.setLength(total_length);
    for (int i = 0; i < omega_length; i++) {
      At_.set(i,1.0);
    }
    for (int i = omega_length; i < total_length; i++) {
      At_.set(i,0.0);
    }

    // Set b
    b_ = 1.0;
 
    // Construct c
    c_.setLength(total_length);
    for (int i = 0; i < omega_length; i++) {
      c_.set(i,-vector_[i]);
    }
    for (int i = omega_length ; i < total_length; i++) {
      c_.set(i,scalar_);
    }

    // Construct W*G
    std::vector<std::shared_ptr<Vector>> WG;
    for (int i = 0; i < omega_length; i++) {
      std::shared_ptr<Vector> product(new Vector(gamma_length_));
      matrix_->matrixVectorProductOfInverse(*vector_list_[i], *product);
      WG.push_back(product);
    }

    // Initialize Q
    Q_.setAsDiagonal(total_length, 1.0);

    // Set (1,1) block of Q
    for (int i = 0; i < omega_length; i++) {
      for (int j = i; j < omega_length; j++) {
        Q_.setElement(i, j, vector_list_[i]->innerProduct(*WG[j]));
      }
    }

    // Set remaining blocks of Q if scalar_ < infinity
    if (scalar_ != NONOPT_DOUBLE_INFINITY) {
      for (int i = omega_length; i < omega_length + gamma_length_; i++) {
        for (int j = 0; j < omega_length; j++) {
          double value = WG[j]->values()[i - omega_length];
          Q_.setElement(i, j, value);
          Q_.setElement(i + gamma_length_, j, -value);
        } // end for
      } // end for
      for (int i = omega_length; i < omega_length + gamma_length_; i++) {
        for (int j = i; j < omega_length  + gamma_length_; j++) {
          double value = matrix_->elementOfInverse(i - omega_length, j - omega_length);
          Q_.setElement(i, j, value);
          Q_.setElement(i + gamma_length_, j + gamma_length_, value);
          Q_.setElement(i + gamma_length_, j, -value);
          if (i != j) {
            Q_.setElement(j + gamma_length_, i, -value);
          } // end if
        } // end for
      } // end for
    } // end for

    // Initialize barrier parameter for initial iterate
    mu_ = barrier_parameter_initial_;

    // Initialize primal-dual solution
    theta_.setLength(total_length);
    theta_.scale(0.0);
    u_ = 0.0;
    v_.setLength(total_length);
    v_.scale(0.0);

    // Initialize omega uniformly
    for (int i = 0; i < omega_length; i++){
      theta_.set(i, 1.0 / omega_length);
      v_.set(i, mu_ / theta_.values()[i]);
    }

    // Initialize sigma and rho if scalar_ < infinity
    if (scalar_ != NONOPT_DOUBLE_INFINITY) {
    
      // Compute G*omega
      Vector Gomega;
      Gomega.setLength(gamma_length_);
      Gomega.scale(0.0);
      for (int i = 0; i < gamma_length_; i++) {
        for (int j = 0; j < omega_length; j++) {
          Gomega.set(i, Gomega.values()[i] += theta_.values()[j] * vector_list_[j]->values()[i]);
        }
      }

      // Initialize sigma and rho
      for (int i = 0; i < gamma_length_; i++) {
        if (Gomega.values()[i] <= -scalar_) {
          theta_.set(i + omega_length, fmax(solution_initialization_factor_, scalar_ - Gomega.values()[i]));
          v_.set(i + omega_length, mu_ / theta_.values()[i + omega_length]);
          theta_.set(i + omega_length + gamma_length_, solution_initialization_factor_);
          v_.set(i + omega_length + gamma_length_, mu_ / theta_.values()[i + omega_length + gamma_length_]);
        } // end if
        else if (Gomega.values()[i] >= scalar_) {
          theta_.set(i + omega_length, solution_initialization_factor_);
          v_.set(i + omega_length, mu_ / theta_.values()[i + omega_length]);
          theta_.set(i + omega_length + gamma_length_, fmax(solution_initialization_factor_, -scalar_ - Gomega.values()[i]));
          v_.set(i + omega_length + gamma_length_, mu_ / theta_.values()[i + omega_length + gamma_length_]);
        } // end else if
        else {
          theta_.set(i + omega_length, solution_initialization_factor_);
          theta_.set(i + omega_length + gamma_length_, solution_initialization_factor_);
          v_.set(i + omega_length, mu_ / theta_.values()[i + omega_length]);
          v_.set(i + omega_length + gamma_length_, mu_ / theta_.values()[i + omega_length + gamma_length_]);
        } // end else
      } // end for

    } // end if

    // Set barrier parameter
    mu_= fmax(barrier_parameter_minimum_, fmin(theta_.innerProduct(v_) / total_length, barrier_parameter_maximum_));

    // Initialize residual vectors
    r_dual_.setLength(total_length);
    r_comp_.setLength(total_length);

    // Set linear system base matrix
    SymmetricMatrixDense LS_matrix_base;
    LS_matrix_base.setAsDiagonal(total_length + 1, 0.0);
    for (int i = 0; i < total_length ; i++) {
      for (int j = i; j < total_length ; j++) {
        LS_matrix_base.setElement(i, j, -Q_.element(i, j));
      } // end for
    } // end for
    for (int i = 0; i < omega_length; i++) {
      LS_matrix_base.setElement(i, total_length, At_.values()[i]);
    }

    // Initialize linear system matrix
    SymmetricMatrixDense LS_matrix;
    LS_matrix.setAsDiagonal(total_length + 1, 0.0);

    // Initialize rewritable quantities
    Vector rhs(total_length + 1, 0.0);
    Vector d(total_length + 1, 0.0);
    Vector dtheta(total_length, 0.0);
    double du = 0.0;
    Vector dv(total_length, 0.0);
    double atheta, au, av;
    Vector ptheta(total_length, 0.0);
    Vector pv(total_length, 0.0);
    double mu_aff;
    double mu_factor;

    // Iteration loop
    while (true) {

      // Print message
      if (iteration_count_ % 20 == 0) {
        reporter->printf(R_QP, R_PER_ITERATION, "========================================\n"
                                                "  Iter.   ||Dual||   ||Prim||   ||Comp||\n"
                                                "========================================\n");
      }
      reporter->printf(R_QP, R_PER_ITERATION, " %6d", iteration_count_);

      // Evaluate primal vectors
      evaluatePrimalVectors();

      // Increment iteration counter
      iteration_count_++;

      // Construct dual residual
      Q_.matrixVectorProduct(theta_, r_dual_);
      r_dual_.addScaledVector( 1.0, c_);
      r_dual_.addScaledVector(- u_,At_);
      r_dual_.addScaledVector(-1.0, v_);

      // Construct primal residual
      r_prim_ = b_ - At_.innerProduct(theta_);

      // Construct complementarity residual
      for (int i = 0; i < total_length; i++) {
        r_comp_.set(i, -theta_.values()[i]*v_.values()[i]);
      }

      // Print message
      reporter->printf(R_QP, R_PER_ITERATION, "  %+.2e  %+.2e  %+.2e", r_dual_.normInf(), fabs(r_prim_), r_comp_.normInf());

      // Evaluate KKT error
      kkt_error_ = fmax(r_dual_.normInf(), fmax(fabs(r_prim_), r_comp_.normInf()));

      // Check for successful solve
      if (kkt_error_ <= kkt_tolerance_) {
        THROW_EXCEPTION(QP_SUCCESS_EXCEPTION, "QP solve successful.");
      }

      // Check for inexact termination
      if (allow_inexact_termination_ &&
          iteration_count_ >= (int)ceil(inexact_termination_initialization_factor_ * (double)vector_.size()) &&
          (iteration_count_ - (int)ceil(inexact_termination_initialization_factor_ * (double)vector_.size())) % inexact_termination_check_interval_ == 0 &&
          inexactTerminationCondition(quantities, reporter)) {
        THROW_EXCEPTION(QP_SUCCESS_EXCEPTION, "QP solve successful.");
      }

      // Check for iteration limit
      if (iteration_count_ >= iteration_limit_) { 
        THROW_EXCEPTION(QP_ITERATION_LIMIT_EXCEPTION, "QP solve unsuccessful. Iteration limit reached.");
      }

      // Check for CPU time limit
      if ((clock() - quantities->startTime()) / (double)CLOCKS_PER_SEC >= quantities->cpuTimeLimit()) {
        THROW_EXCEPTION(QP_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
      }

      ////////////////////
      // PREDICTOR STEP //
      ////////////////////

      // Set linear system matrix
      for (int i = 0; i < total_length + 1; i++) {
        for (int j = i; j < total_length + 1; j++) {
          if (j == i && i < total_length) {
            LS_matrix.setElement(i, i, -Q_.element(i, i) - v_.values()[i] / theta_.values()[i]);
          }
          else {
            LS_matrix.setElement(i, j, LS_matrix_base.element(i,j));
          }
        }
      }

      // Set right-hand side vector
      for (int i = 0; i < total_length; i++) {
        rhs.set(i, r_dual_.values()[i] + v_.values()[i]);
      }
      rhs.set(total_length, r_prim_);

      // Store factorization values
      int* ipiv = new int[total_length + 1];

      // Solve linear system
      solveLinearSystem(total_length + 1, LS_matrix.valuesModifiable(), rhs.valuesModifiable(), ipiv, d.valuesModifiable());

      // Parse solution
      for (int i = 0; i < total_length; i++) {
        dtheta.set(i, d.values()[i]);
        dv.set(i, -v_.values()[i] - v_.values()[i]*dtheta.values()[i]/theta_.values()[i]);
      }
      du = d.values()[total_length];

      // Calculate step sizes
      computeStepSizes(dtheta, du, dv, atheta, au, av);

      // Set predictor values
      ptheta.linearCombination(1.0, theta_, atheta, dtheta);
      pv.linearCombination(1.0, v_, av, dv);

      // Set mu_factor
      mu_aff = barrier_parameter_factor_ * ptheta.innerProduct(pv) / (double)total_length;
      mu_factor = fmax(barrier_parameter_minimum_ / mu_, fmin(pow(mu_aff / mu_, 3), barrier_parameter_maximum_ / mu_));
      
      ////////////////////
      // CORRECTOR STEP //
      ////////////////////

      // Set right-hand side vector
      for (int i = 0; i < total_length; i++) {
        rhs.set(i, r_dual_.values()[i] + v_.values()[i] - mu_factor * mu_ / theta_.values()[i]);
      }
      rhs.set(total_length, r_prim_);

      // Solve linear system
      solveLinearSystemReuseFactorization(total_length + 1, LS_matrix.valuesModifiable(), rhs.valuesModifiable(), ipiv, d.valuesModifiable());

      // Delete factorization info
      delete[] ipiv;

      // Parse solution
      for (int i = 0; i < total_length; i++) {
        dtheta.set(i, d.values()[i]);
        dv.set(i, -v_.values()[i] + mu_factor * mu_ / theta_.values()[i] - v_.values()[i] * dtheta.values()[i] / theta_.values()[i]);
      }
      du = d.values()[total_length];

      // Calculate step sizes
      computeStepSizes(dtheta, du, dv, atheta, au, av);

      // Update the elements
      theta_.addScaledVector(atheta, dtheta);
      u_ += au * du;
      v_.addScaledVector(av, dv);
              
      // Update barrier parameter
      mu_ = fmax(barrier_parameter_minimum_, fmin(theta_.innerProduct(v_) / total_length, barrier_parameter_maximum_));

    } //end while

  } // end try

  // catch exceptions
  catch (QP_SUCCESS_EXCEPTION& exec) {
    setStatus(QP_SUCCESS);
  } catch (QP_CPU_TIME_LIMIT_EXCEPTION& exec) {
    setStatus(QP_CPU_TIME_LIMIT);
    THROW_EXCEPTION(NONOPT_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
  } catch (QP_FACTORIZATION_ERROR_EXCEPTION& exec) {
    setStatus(QP_FACTORIZATION_ERROR);
  } catch (QP_INPUT_ERROR_EXCEPTION& exec) {
    setStatus(QP_INPUT_ERROR);
  } catch (QP_ITERATION_LIMIT_EXCEPTION& exec) {
    setStatus(QP_ITERATION_LIMIT);
  } catch (QP_NAN_ERROR_EXCEPTION& exec) {
    setStatus(QP_NAN_ERROR);
  }

  // Print new line
  reporter->printf(R_QP, R_PER_ITERATION, "\n");

  // Print finalizing
  reporter->printf(R_QP, R_PER_INNER_ITERATION, "Finalizing solution\n");

  // Finalize solution
  finalizeSolution();

} // end solveQP

// Solve hot
void QPSolverInteriorPoint::solveQPHot(const Options* options,
                                       const Reporter* reporter,
                                       Quantities* quantities)
{

  // No hot start, just solve QP
  solveQP(options,reporter,quantities);

} // end solveQPHot

// Print method
void QPSolverInteriorPoint::printData(const Reporter* reporter)
{
  // Print matrix values
  reporter->printf(R_QP, R_BASIC, "\nMATRIX:\n");
  for (int i = 0; i < gamma_length_; i++) {
    for (int j = 0; j < gamma_length_; j++) {
      reporter->printf(R_QP, R_BASIC, " %+23.16e", matrix_->elementOfInverse(i, j));
    }
    reporter->printf(R_QP, R_BASIC, "\n");
  } // end for

  // Print vector list
  reporter->printf(R_QP, R_BASIC, "VECTOR LIST:\n");
  for (int j = 0; j < gamma_length_; j++) {
    for (int i = 0; i < (int)vector_list_.size(); i++) {
      reporter->printf(R_QP, R_BASIC, " %+23.16e", vector_list_[i]->values()[j]);
    }
    reporter->printf(R_QP, R_BASIC, "\n");
  } // end for

  // Print vector
  reporter->printf(R_QP, R_BASIC, "VECTOR:\n");
  for (int i = 0; i < (int)vector_.size(); i++) {
    reporter->printf(R_QP, R_BASIC, " %+23.16e\n", vector_[i]);
  }

  // Print scalar
  reporter->printf(R_QP, R_BASIC, "SCALAR:\n");
  reporter->printf(R_QP, R_BASIC, " %+23.16e\n", scalar_);

} // end printData

// Check quantities
bool QPSolverInteriorPoint::checkQuantityCompatibility()
{

  // Check for null pointer to matrix
  if (matrix_ == nullptr) {
    return false;
  }

  // Check matrix size
  if (gamma_length_ != matrix_->size()) {
    return false;
  }

  // Loop through vector list
  for (int i = 0; i < (int)vector_list_.size(); i++) {

    // Check vector size
    if (gamma_length_ != vector_list_[i]->length()) {
      return false;
    }

  } // end for

  // Check length of vector list and vector
  if ((int)vector_list_.size() != (int)vector_.size()) {
    return false;
  }

  // Check scalar (should be positive)
  if (scalar_ <= 0.0) {
    return false;
  }

  // Return success
  return true;

} // end checkQuantityCompatibility

// Evaluate dual vectors
void QPSolverInteriorPoint::evaluatePrimalVectors()
{

  // Set omega length
  int omega_length = (int)vector_.size();
  
  // Zero-out vectors
  omega_.setLength(omega_length);
  omega_.scale(0.0);
  gamma_.setLength(gamma_length_);
  gamma_.scale(0.0);

  // Set omega values
  for (int i = 0; i < omega_length; i++) {
    omega_.set(i, theta_.values()[i]);
  }

  // Set gamma values
  for (int i = 0; i < gamma_length_; i++) {
    gamma_.set(i, theta_.values()[i + omega_length] - theta_.values()[i + omega_length + gamma_length_]);
  }

  // Zero-out vectors
  combination_.scale(0.0);
  combination_translated_.scale(0.0);
  primal_solution_.scale(0.0);
  primal_solution_feasible_.scale(0.0);

  // Loop to compute gradient combination
  for (int i = 0; i < omega_length; i++) {
    combination_.addScaledVector(omega_.values()[i], *vector_list_[i]);
  }

  // Initialize gradient combination shifted
  combination_translated_.linearCombination(1.0, combination_, 1.0, gamma_);

  // Compute matrix-vector product
  matrix_->matrixVectorProductOfInverse(combination_translated_, primal_solution_);
  primal_solution_.scale(-1.0);

  // Compute feasible primal step by projection
  primal_solution_feasible_.copy(primal_solution_);
  primal_solution_projection_scalar_ = fmin(1.0, scalar_ / primal_solution_.normInf());
  primal_solution_feasible_.scale(primal_solution_projection_scalar_);

} // end evaluatePrimalVectors

// Finalize solution
void QPSolverInteriorPoint::finalizeSolution()
{

  // Evaluate primal vectors
  evaluatePrimalVectors();

} // end finalizeSolution

// Solve linear system
void QPSolverInteriorPoint::solveLinearSystem(int size,
                                              double matrix[],
                                              double right_hand_side[],
                                              int ipiv[],
                                              double solution[])
{

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  int n = size;
  int nrhs = 1;
  int info;
  int* ipiv_local = new int[size];
    
  // Query for optimal workspace
  double work_query;
  int lwork = -1;
  dsysv_(&upper_lower, &n, &nrhs, matrix, &n, ipiv_local, solution, &n, &work_query, &lwork, &info);
    
  // Allocate optimal workspace
  lwork = static_cast<int>(work_query);
  double* work = new double[lwork];

  // Copy right_hand_side to solution
  int increment1 = 1;
  dcopy_(&n, right_hand_side, &increment1, solution, &increment1);

  // Solve system
  dsysv_(&upper_lower, &n, &nrhs, matrix, &n, ipiv_local, solution, &n, work, &lwork, &info);

  // Copy ipiv_local to ipiv
  for (int i = 0; i < n; i++) {
    ipiv[i] = ipiv_local[i];
  }

  // Delete arrays
  delete[] ipiv_local;
  delete[] work;

} // end solveLinearSystem

// Solve linear system
void QPSolverInteriorPoint::solveLinearSystemReuseFactorization(int size,
                                                                double matrix[],
                                                                double right_hand_side[],
                                                                int ipiv[],
                                                                double solution[])
{

  // Set inputs for BLASLAPACK
  char upper_lower = 'L';
  int n = size;
  int nrhs = 1;
  int info;

  // Copy right_hand_side to solution
  int increment1 = 1;
  dcopy_(&n, right_hand_side, &increment1, solution, &increment1);

  // Solve system
  dsytrs_(&upper_lower, &n, &nrhs, matrix, &n, ipiv, solution, &n, &info);

} // end solveLinearSystem

void QPSolverInteriorPoint::computeStepSizes(const Vector& dtheta, 
                                             const double& du, 
                                             const Vector& dv, 
                                             double& atheta, 
                                             double& au, 
                                             double& av)
{

  // Set total length
  int total_length = theta_.length();

  // Compute r
  Vector r(total_length + 1, 0.0);
  r.set(0, r_prim_);
  for (int i = 1; i < total_length + 1; i++) {
    r.set(i, r_dual_.values()[i - 1]);
  }

  // Compute s
  Vector s(total_length + 1, 0.0);
  s.set(0, -At_.innerProduct(dtheta));
  Vector Qdtheta(total_length, 0.0);
  Q_.matrixVectorProduct(dtheta, Qdtheta);
  for (int i = 1; i < total_length + 1; i++) {
    s.set(i, Qdtheta.values()[i - 1]);
  }

  // Compute o
  Vector o(total_length + 1, 0.0);
  for (int i = 1; i < total_length + 1; i++) {
    o.set(i, -At_.values()[i-1]);
  }
  o.scale(du);

  // Compute p
  Vector p(total_length + 1, 0.0);
  for (int i = 1; i < total_length + 1; i++) {
    p.set(i, -dv.values()[i - 1]);
  }

   // Impose fraction-to-the-boundary rule
  double bar_atheta = 1.0;
  double bar_av = 1.0;
  for (int i = 0; i < total_length; i++){
    if (dtheta.values()[i] < 0) {
      bar_atheta = fmin(bar_atheta, (fraction_to_boundary_tolerance_ - 1.0) * theta_.values()[i] / dtheta.values()[i]);
    }
    if (dv.values()[i] < 0) {
      bar_av = fmin(bar_av, (fraction_to_boundary_tolerance_ - 1.0) * v_.values()[i] / dv.values()[i]);
    }
  }
  double bar_alpha = fmin(bar_atheta, bar_av);

  // Solve preliminary 2d QP
  Vector s_plus_p(total_length + 1, 0.0);
  s_plus_p.linearCombination(1.0, s, 1.0, p);
  double Q1_11 = s_plus_p.innerProduct(s_plus_p) + dtheta.innerProduct(dv);
  double Q1_12 = s_plus_p.innerProduct(o);
  double Q1_22 = o.innerProduct(o);
  double c1_1 = r.innerProduct(s_plus_p) + 0.5 * (dtheta.innerProduct(v_) + theta_.innerProduct(dv));
  double c1_2 = r.innerProduct(o);
  double atheta1 = fmax(0.0, fmin((c1_2 * Q1_12 / Q1_22 - c1_1) / (Q1_11 - Q1_12 * Q1_12 / Q1_22), bar_alpha));
  double au1 = (-c1_2 - atheta1 * Q1_12) / Q1_22;
  double av1 = atheta1;

  // Solve secondary 2d QP
  double atheta2, au2, av2;
  if (bar_atheta <= bar_av) {
    double Q2_11 = o.innerProduct(o);
    double Q2_12 = o.innerProduct(p);
    double Q2_22 = p.innerProduct(p);
    double c2_1 = r.innerProduct(o) + bar_atheta * s.innerProduct(o);
    double c2_2 = r.innerProduct(p) + 0.5 * theta_.innerProduct(dv) + bar_atheta * (s.innerProduct(p) + 0.5 * dtheta.innerProduct(dv));
    atheta2 = bar_atheta;
    av2 = fmax(bar_alpha, fmin((c2_1 * Q2_12 / Q2_11 - c2_2) / (Q2_22 - Q2_12 * Q2_12 / Q2_11), bar_av));
    au2 = (-c2_1 - av2 * Q2_12) / Q2_11;
  }
  else {
    double Q2_11 = s.innerProduct(s);
    double Q2_12 = s.innerProduct(o);
    double Q2_22 = o.innerProduct(o);
    double c2_1 = r.innerProduct(s) + 0.5 * dtheta.innerProduct(v_) + bar_av * (s.innerProduct(p) + 0.5 * dtheta.innerProduct(dv));
    double c2_2 = r.innerProduct(o) + bar_av * o.innerProduct(p);
    atheta2 = fmax(bar_alpha, fmin((c2_2 * Q2_12 / Q2_22 - c2_1) / (Q2_11 - Q2_12 * Q2_12 / Q2_22), bar_atheta));
    au2 = (-c2_2 - atheta2 * Q2_12) / Q2_22;
    av2 = bar_av;
  }

  // Compute trial points
  Vector theta1(total_length, 0.0);
  theta1.linearCombination(1.0, theta_, atheta1, dtheta);
  Vector theta2(total_length, 0.0);
  theta2.linearCombination(1.0, theta_, atheta2, dtheta);
  double u1 = u_ + au1 * du;
  double u2 = u_ + au2 * du;
  Vector v1(total_length, 0.0);
  v1.linearCombination(1.0, v_, av1, dv);
  Vector v2(total_length, 0.0);
  v2.linearCombination(1.0, v_, av2, dv);

  // Evaluate merit function values
  Vector r_dual1(total_length, 0.0);
  Q_.matrixVectorProduct(theta1, r_dual1);
  r_dual1.addScaledVector(1.0, c_);
  r_dual1.addScaledVector(-u1, At_);
  r_dual1.addScaledVector(-1.0, v1);
  double r_prim1 = At_.innerProduct(theta1) - b_;
  Vector r_dual2(total_length, 0.0);
  Q_.matrixVectorProduct(theta2, r_dual2);
  r_dual2.addScaledVector(1.0, c_);
  r_dual2.addScaledVector(-u2, At_);
  r_dual2.addScaledVector(-1.0, v2);
  double r_prim2 = At_.innerProduct(theta2) - b_;

  // Choose best step sizes
  if (pow(r_prim1, 2.0) + pow(r_dual1.norm2(), 2.0) + theta1.innerProduct(v1)
    < pow(r_prim2, 2.0) + pow(r_dual2.norm2(), 2.0) + theta2.innerProduct(v2)) {
    atheta = atheta1;
    au = au1;
    av = av1;
  }
  else {
    atheta = atheta2;
    au = au2;
    av = av2;
  }

} // end computeStepSizes

} // namespace NonOpt