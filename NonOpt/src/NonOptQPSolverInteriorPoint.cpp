// Copyright (C) 2024 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Lara Zebiane

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <iomanip> // Added by Lara (?)

#include "NonOptBLASLAPACK.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptQPSolverInteriorPoint.hpp"

namespace NonOpt
{

// Constructor
QPSolverInteriorPoint::QPSolverInteriorPoint()
  : inexact_solution_tolerance_(0.0),
    scalar_(0.0),
    factor_(nullptr),
    inner_solution_1_(nullptr),
    inner_solution_2_(nullptr),
    inner_solution_3_(nullptr),
    inner_solution_ls_(nullptr),
    inner_solution_trial_(nullptr),
    new_system_vector_(nullptr),
    right_hand_side_(nullptr),
    system_solution_(nullptr),
    system_solution_best_(nullptr) {}

// Destructor
QPSolverInteriorPoint::~QPSolverInteriorPoint()
{

  // Delete arrays
  if (inner_solution_1_ != nullptr) {
    delete[] inner_solution_1_;
    inner_solution_1_ = nullptr;
  } // end if
  if (inner_solution_2_ != nullptr) {
    delete[] inner_solution_2_;
    inner_solution_2_ = nullptr;
  } // end if
  if (inner_solution_3_ != nullptr) {
    delete[] inner_solution_3_;
    inner_solution_3_ = nullptr;
  } // end if
  if (inner_solution_ls_ != nullptr) {
    delete[] inner_solution_ls_;
    inner_solution_ls_ = nullptr;
  } // end if
  if (inner_solution_trial_ != nullptr) {
    delete[] inner_solution_trial_;
    inner_solution_trial_ = nullptr;
  } // end if
  if (new_system_vector_ != nullptr) {
    delete[] new_system_vector_;
    new_system_vector_ = nullptr;
  } // end if
  if (right_hand_side_ != nullptr) {
    delete[] right_hand_side_;
    right_hand_side_ = nullptr;
  } // end if
  if (system_solution_ != nullptr) {
    delete[] system_solution_;
    system_solution_ = nullptr;
  } // end if
  if (system_solution_best_ != nullptr) {
    delete[] system_solution_best_;
    system_solution_best_ = nullptr;
  } // end if
  if (factor_ != nullptr) {
    delete[] factor_;
    factor_ = nullptr;
  } // end if

} // end destructor

// Add options
void QPSolverInteriorPoint::addOptions(Options* options)
{

  // Add bool options
  options->addBoolOption("QPDAS_fail_on_factorization_error",
                         false,
                         "Indicator for whether to indicate failure on factorization error.\n"
                         "Default     : false.");
  options->addBoolOption("QPDAS_allow_inexact_termination",
                         false,
                         "Indicator for whether to allow early termination.\n"
                         "Default     : false.");
  // Add double options
  options->addDoubleOption("QPDAS_kkt_tolerance",
                           1e-08,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining optimality.  If the KKT error\n"
                           "              computed by the algorithm falls below this tolerance, then\n"
                           "              the algorithm terminates with a message of success.\n"
                           "Default     : 1e-08.");
  options->addDoubleOption("QPDAS_cholesky_tolerance",
                           1e-12,
                           0.0,
                           1.0,
                           "Tolerance for small (or negative values) when updating the\n"
                           "              Cholesky factorization of an augmented matrix.  If a diagonal\n"
                           "              value in the factorization falls below this tolerance, then\n"
                           "              the factorization may be re-computed from scratch and/or the\n"
                           "              diagonal value may be replaced by this tolerance value.\n"
                           "Default     : 1e-12.");
  options->addDoubleOption("QPDAS_inexact_termination_descent_tolerance",
                           1e-04,
                           0.0,
                           1.0,
                           "Descent direction tolerance for inexactness conditions.\n"
                           "Default     : 1e-04.");
  options->addDoubleOption("QPDAS_inexact_termination_initialization_factor",
                           2.5e-01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Determines number of iterations (fraction of gradient list size)\n"
                           "              to perform before inexact termination conditions are checked.\n"
                           "Default     : 2.5e-01.");
  options->addDoubleOption("QPDAS_inexact_termination_ratio_minimum",
                           1e-02,
                           0.0,
                           1.0,
                           "Minimum value for ratio used in inexact termination condition.\n"
                           "Default     : 1e-02.");
  options->addDoubleOption("QPDAS_linear_independence_tolerance",
                           1e-12,
                           0.0,
                           1.0,
                           "Tolerance for a linear independence check when adding a new\n"
                           "              vector to an augmented matrix.  If a residual value falls\n"
                           "              below this tolerance, then indices are exchanged to maintain\n"
                           "              linear independence of the augmented matrix.\n"
                           "Default     : 1e-12.");
  options->addDoubleOption("QPIPM_beta",
                           1e-2,
                           0.0,
                           1.0,
                           "beta has to be 0 and 1. TBC...\n"
                           "Default     : 0.5.");
  options->addDoubleOption("QPIPM_sigma",
                           0.002,
                           0.0,
                           1.0,
                           "sigma has to be 0 and 1. TBC...\n"
                           "Default     : 0.005.");   
  options->addDoubleOption("QPIPM_eps",
                           1e-16,
                           0.0,
                           1.0,
                           "eps is a very small number. TBC...\n"
                           "Default     : 1e-16.");                 
  options->addDoubleOption("QPIPM_tol_in",
                           1e-5,
                           0.0,
                           1.0,
                           "Tolerance of inner loop. TBC...\n"
                           "Default     : 1e-3.");    
  options->addDoubleOption("QPIPM_tol_out",
                           1e-6,
                           0.0,
                           1.0,
                           "Tolerance of outer loop. TBC...\n"
                           "Default     : 1e-4.");            
  options->addDoubleOption("QPIPM_mu_factor",
                           0.2,
                           0.0,
                           1.0,
                           "The mu factor has to be 0 and 1. TBC...\n"
                           "Default     : 0.5.");                                                                                                          
  // Add integer options
  options->addIntegerOption("QPDAS_inexact_termination_check_interval",
                            4,
                            1,
                            NONOPT_INT_INFINITY,
                            "Number of iterations to perform between checks of inexact\n"
                            "              termination conditions.\n"
                            "Default     : 4.");
  options->addIntegerOption("QPDAS_iteration_limit_minimum",
                            1e+00,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Minimum limit on the number of iterations.\n"
                            "Default     : 1e+00.");
  options->addIntegerOption("QPDAS_iteration_limit_maximum",
                            1e+06,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Maximum limit on the number of iterations.\n"
                            "Default     : 1e+06.");
  options->addIntegerOption("QPIPM_max_iter_in",
                            500,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Maximum number of iterations for the inner loop.\n"
                            "Default     : 500.");                   
  options->addIntegerOption("QPIPM_max_iter_out",
                            750,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Maximum number of iterations for the outer loop.\n"
                            "Default     : 150.");                                    

} // end addOptions

// Set options
void QPSolverInteriorPoint::setOptions(Options* options)
{

  // Read bool options
  options->valueAsBool("QPDAS_fail_on_factorization_error", fail_on_factorization_error_);
  options->valueAsBool("QPDAS_allow_inexact_termination", allow_inexact_termination_);

  // Read double options
  options->valueAsDouble("QPDAS_kkt_tolerance", kkt_tolerance_);
  options->valueAsDouble("QPDAS_cholesky_tolerance", cholesky_tolerance_);
  options->valueAsDouble("QPDAS_inexact_termination_descent_tolerance", inexact_termination_descent_tolerance_);
  options->valueAsDouble("QPDAS_inexact_termination_initialization_factor", inexact_termination_initialization_factor_);
  options->valueAsDouble("QPDAS_inexact_termination_ratio_minimum", inexact_termination_ratio_minimum_);
  options->valueAsDouble("QPDAS_linear_independence_tolerance", linear_independence_tolerance_);
  options->valueAsDouble("QPIPM_beta", beta_);
  options->valueAsDouble("QPIPM_sigma", sigma_);
  options->valueAsDouble("QPIPM_eps", eps_);
  options->valueAsDouble("QPIPM_tol_in", tol_in_);
  options->valueAsDouble("QPIPM_tol_out", tol_out_);
  options->valueAsDouble("QPIPM_mu_factor", mu_factor_);  
  // Read integer options
  options->valueAsInteger("QPDAS_inexact_termination_check_interval", inexact_termination_check_interval_);
  options->valueAsInteger("QPDAS_iteration_limit_minimum", iteration_limit_minimum_);
  options->valueAsInteger("QPDAS_iteration_limit_maximum", iteration_limit_maximum_);
  options->valueAsInteger("QPIPM_max_iter_in", max_iter_in_);
  options->valueAsInteger("QPIPM_max_iter_out", max_iter_out_);
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
  system_solution_length_ = 3 * gamma_length_ + 1;
  factor_length_ = 2 * system_solution_length_; // system_solution_length_ * system_solution_length_;

  // Initialize problem data
  matrix_ = nullptr;
  vector_list_.clear();
  vector_.clear();
  scalar_ = -1.0;
  
  // std::cout << "Hi1!!" << std::endl;

  // Delete arrays (in case they exist)
  if (inner_solution_1_ != nullptr) {
    delete[] inner_solution_1_;
    inner_solution_1_ = nullptr;
  } // end if
  if (inner_solution_2_ != nullptr) {
    delete[] inner_solution_2_;
    inner_solution_2_ = nullptr;
  } // end if
  if (inner_solution_3_ != nullptr) {
    delete[] inner_solution_3_;
    inner_solution_3_ = nullptr;
  } // end if
  if (inner_solution_ls_ != nullptr) {
    delete[] inner_solution_ls_;
    inner_solution_ls_ = nullptr;
  } // end if
  if (inner_solution_trial_ != nullptr) {
    delete[] inner_solution_trial_;
    inner_solution_trial_ = nullptr;
  } // end if
  if (new_system_vector_ != nullptr) {
    delete[] new_system_vector_;
    new_system_vector_ = nullptr;
  } // end if
  if (right_hand_side_ != nullptr) {
    delete[] right_hand_side_;
    right_hand_side_ = nullptr;
  } // end if
  if (system_solution_ != nullptr) {
    delete[] system_solution_;
    system_solution_ = nullptr;
  } // end if
  if (system_solution_best_ != nullptr) {
    delete[] system_solution_best_;
    system_solution_best_ = nullptr;
  } // end if
  if (factor_ != nullptr) {
    delete[] factor_;
    factor_ = nullptr;
  } // end if

  // Allocate arrays
  inner_solution_1_ = new double[system_solution_length_];
  inner_solution_2_ = new double[system_solution_length_];
  inner_solution_3_ = new double[system_solution_length_];
  inner_solution_ls_ = new double[system_solution_length_];
  inner_solution_trial_ = new double[system_solution_length_];
  new_system_vector_ = new double[system_solution_length_];
  right_hand_side_ = new double[system_solution_length_];
  system_solution_ = new double[system_solution_length_];
  system_solution_best_ = new double[system_solution_length_];
  factor_ = new double[factor_length_];

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

// Get KKT error from dual solution
// double QPSolverInteriorPoint::KKTErrorDual() //remove
// {

//   // Evaluate gradient combination
//   Vector Gomega(gamma_length_, 0.0);
//   for (int i = 0; i < (int)vector_list_.size(); i++) {
//     Gomega.addScaledVector(omega_.values()[i], *vector_list_[i].get());
//   }

//   // Declare dual vector
//   Vector d0(gamma_length_, 0.0);
//   d0.addScaledVector(-1.0, Gomega);
//   d0.addScaledVector(-1.0, gamma_);
//   Vector d(gamma_length_);
//   matrix_->matrixVectorProductOfInverse(d0, d);

//   // Evaluate gradient inner products
//   Vector Gtd;
//   Gtd.setLength((int)vector_.size());
//   for (int i = 0; i < (int)vector_list_.size(); i++) {
//     Gtd.values()[i] = vector_list_[i]->innerProduct(d);
//   }

//   // Evaluate dual scalar
//   double z = omega_.innerProduct(Gtd);
//   for (int i = 0; i < (int)vector_list_.size(); i++) {
//     z = z + omega_.values()[i] * vector_[i];
//   }

//   // Initialize KKT error
//   double kkt_error = 0.0;

//   // Loop through KKT conditions, update error
//   for (int j = 0; j < (int)vector_.size(); j++) {
//     double term = -z + vector_[j] + Gtd.values()[j];
//     if (term > 0.0 && term > kkt_error) {
//       kkt_error = term;
//     }
//   } // end for
//   for (int i = 0; i < gamma_length_; i++) {
//     double term = -scalar_ + d.values()[i];
//     if (term > 0.0 && term > kkt_error) {
//       kkt_error = term;
//     }
//     term = -scalar_ - d.values()[i];
//     if (term > 0.0 && term > kkt_error) {
//       kkt_error = term;
//     }
//   } // end for
//   for (int j = 0; j < (int)vector_.size(); j++) {
//     double term = -omega_.values()[j];
//     if (term > 0.0 && term > kkt_error) {
//       kkt_error = term;
//     }
//   } // end for
//   double sum = 0.0;
//   for (int j = 0; j < (int)vector_.size(); j++) {
//     sum = sum + omega_.values()[j];
//   }
//   if (fabs(sum - 1.0) > kkt_error) {
//     kkt_error = fabs(sum - 1.0);
//   }
//   for (int i = 0; i < gamma_length_; i++) {
//     double term1 = gamma_.values()[i];
//     double term2 = scalar_ - d.values()[i];
//     if (term1 > 0.0 && term1 * term2 > kkt_error) {
//       kkt_error = term1 * term2;
//     }
//     term1 = -gamma_.values()[i];
//     term2 = scalar_ + d.values()[i];
//     if (term1 > 0.0 && term1 * term2 > kkt_error) {
//       kkt_error = term1 * term2;
//     }
//   } // end for

//   // Set maximum
//   if (std::isnan(kkt_error) || kkt_error > NONOPT_DOUBLE_INFINITY) {
//     kkt_error = NONOPT_DOUBLE_INFINITY;
//   }

//   // Return
//   return kkt_error;

// } // end KKTErrorDual

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

  // Clear solution arrays
  omega_positive_.clear();
  omega_positive_best_.clear();
  gamma_positive_.clear();
  gamma_positive_best_.clear();
  gamma_negative_.clear();
  gamma_negative_best_.clear();

  // Initialize dual values
  multiplier_ = 0.0;
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

  // Finalize solution to set omega and gamma
  finalizeSolution(reporter);

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

// Lara's header

void printHeader(const Reporter* reporter) {
    reporter->printf(R_QP, R_PER_INNER_ITERATION,
        "%-20s%-20s%-25s%-15s%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-20s\n",
        "Out_Ite",  "Min_Val_x", "Min_Val_z", "mu", 
        "Inf||r_p||", "Inf||r_d||", "Inf||r_c||",
        "Inf||d_x||", "Inf||d_y||", "Inf||d_z||", 
         "alph_x", "alph_y", "alph_z"
    );
}



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

    // ********* Lara's START! *********

    int omega_length = (int)vector_.size();

    int total_length;

    if (scalar_ == NONOPT_DOUBLE_INFINITY) {
        total_length = omega_length;
    } else {
        total_length = 2 * gamma_length_ + omega_length;
    }

    
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma length:  %d\n", gamma_length_);
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "vector list size length:  %d\n", vector_list_.size());
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "vector size length:  %d\n", vector_.size());
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "omega length:  %d\n", omega_length);
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "total length:  %d\n", total_length);

    printHeader(reporter);

    matrix_At_.setLength(total_length);
    matrix_At_.scale(0.0);
    for (int i = 0; i < omega_length; i++){
        matrix_At_.set(i,1.0);
    }

    scalar_b_ *= 0;
    scalar_b_ = 1.0;
 
    vector_c_.setLength(total_length);
    vector_c_.scale(0.0); 
    for (int i = 0; i < omega_length; i++){
        vector_c_.set(i,-vector_[i]);
    }
    for (int i = omega_length ; i < total_length; i++){
        vector_c_.set(i,scalar_);
    }   


    // Construct matrix Q
    std::vector<std::shared_ptr<Vector>> WG;
    for (auto& elem : WG) {
        elem = nullptr;
    }
    for (int i = 0; i < omega_length; i++) {
      std::shared_ptr<Vector> product(new Vector(gamma_length_));
      matrix_->matrixVectorProductOfInverse(*vector_list_[i], *product);
      WG.push_back(product);
    }

    matrix_Q_.setAsDiagonal(total_length, 1.0);

    for (int i = 0; i < omega_length; i++){
          for (int j = i; j < omega_length; j++){
            double value = vector_list_[i]->innerProduct(*WG[j]);
            matrix_Q_.setElement(i, j, value);
          }
    }

    if (scalar_ != NONOPT_DOUBLE_INFINITY){
      for (int i = omega_length; i < omega_length + gamma_length_; i++){
            for (int j = 0; j < omega_length; j++){
              double value = WG[j]->values()[i-omega_length];
              matrix_Q_.setElement(i, j, value);
              matrix_Q_.setElement(i+gamma_length_, j, -value);
            }
      }  
      for (int i = omega_length; i < omega_length + gamma_length_; i++){
            for (int j = i; j < omega_length  + gamma_length_; j++){
              double value = matrix_->elementOfInverse(i-omega_length, j-omega_length);
              matrix_Q_.setElement(i, j, value);
              matrix_Q_.setElement(i + gamma_length_, j + gamma_length_, value);
              matrix_Q_.setElement(i + gamma_length_, j, -value);
              if (i != j){
                matrix_Q_.setElement(j + gamma_length_, i, -value);
              }
            }
      }
    }

    mu_ = 1e-1;
    vector_x_.setLength(total_length);
    vector_x_.scale(0.0);
    vector_z_.setLength(total_length);
    vector_z_.scale(0.0);
    scalar_y_ *= 0.0;

    for (int i = 0; i < omega_length; i++){
      vector_x_.set(i, 1.0 / omega_length);
      vector_z_.set(i, 1.0);
    }

    Vector Gomega_new;
    Gomega_new.setLength(gamma_length_);
    Gomega_new.scale(0.0);
    for (int i = 0; i < gamma_length_; i++) {
      for (int j = 0; j < omega_length; j++){
      Gomega_new.set(i, Gomega_new.values()[i] += vector_x_.values()[j] * vector_list_[j]->values()[i] );}
    }

    for (int i = 0; i < gamma_length_; i++){
      if (Gomega_new.values()[i] <= - scalar_){
        vector_x_.set(i + omega_length, fmax(init_param_, scalar_ - Gomega_new.values()[i]) );
        vector_z_.set(i + omega_length, mu_ / vector_x_.values()[i + omega_length]);
        vector_x_.set(i + omega_length + gamma_length_, init_param_);
        vector_z_.set(i + omega_length + gamma_length_, mu_ / vector_x_.values()[i + omega_length + gamma_length_]);
      }
      else{
        if (Gomega_new.values()[i] >= scalar_)
        {
          vector_x_.set(i + omega_length, init_param_);
          vector_z_.set(i + omega_length, mu_ / vector_x_.values()[i + omega_length]);
          vector_x_.set(i + omega_length + gamma_length_, fmax(init_param_, -scalar_ - Gomega_new.values()[i]));
          vector_z_.set(i + omega_length + gamma_length_, mu_ / vector_x_.values()[i + omega_length + gamma_length_]);
        }
        else{
          vector_x_.set(i + omega_length, init_param_);
          vector_x_.set(i + omega_length + gamma_length_, init_param_);
          vector_z_.set(i + omega_length, mu_ / vector_x_.values()[i + omega_length]);
          vector_z_.set(i + omega_length + gamma_length_, mu_ / vector_x_.values()[i + omega_length + gamma_length_]);
        }
        
      }
    }



    mu_= fmax( 1e-10, vector_x_.innerProduct(vector_z_) / total_length) ;

    iter_count_ = 0;

    while (iter_count_ < max_iter_out_)
    {
      iter_count_++;

      if (iter_count_ % 20 == 0) {
          printHeader(reporter);
      }


      // Constructing residuals here
      Vector QX(total_length, 0.0);
      matrix_Q_.matrixVectorProduct(vector_x_, QX);
      Vector oper1(total_length, 0.0);
      oper1.copy(QX);
      oper1.addScaledVector(-scalar_y_, matrix_At_);
      r_dual_.setLength(total_length);
      r_dual_.scale(0.0);
      r_dual_.copy(vector_c_);  
      r_dual_.addScaledVector(-1, vector_z_);
      r_dual_.addScaledVector(1, oper1);

      r_cent_.setLength(total_length);
      r_cent_.scale(0.0);
      for ( int i=0 ; i < total_length; i++){
        r_cent_.set(i, -vector_x_.values()[i]*vector_z_.values()[i]);
      }

      r_pri_ *= 0;
      r_pri_=-matrix_At_.innerProduct(vector_x_) + scalar_b_;

      Vector residual( 2 * total_length + 1, 0.0);
      for (int i = 0; i < total_length; i++){
        residual.set(i, r_dual_.values()[i]);
        residual.set(i + total_length + 1, r_cent_.values()[i] );
      }
      residual.set(total_length, r_pri_);

      if (residual.normInf() <= tol_out_){ //kkterror
        kkt_error_ = 0;
        kkt_error_ = residual.normInf();
        break;
      }
      

      // ******* Predictor Step *******
      SymmetricMatrixDense LHS_matrix_P;
      LHS_matrix_P.setAsDiagonal(total_length + 1, 0);
      for (int i = 0; i < total_length ; i++){
        for (int j = i; j < total_length ; j++){
          LHS_matrix_P.setElement(i, j, -matrix_Q_.element(i, j));
        }
        LHS_matrix_P.setElement(i, i, LHS_matrix_P.element(i,i) - vector_z_.values()[i] / vector_x_.values()[i] );
      }
      for (int i = 0; i < omega_length ; i++){
        LHS_matrix_P.setElement(i, total_length, 1);
      }


      Vector r_mu(total_length, 0.0);
      r_mu.copy( r_cent_ ); 
      Vector RHS_vector_P(total_length + 1, 0.0);
      for (int i = 0; i < total_length; i++){
        RHS_vector_P.set( i, r_dual_.values()[i] - r_mu.values()[i] / vector_x_.values()[i]);
      }
      RHS_vector_P.set(total_length, r_pri_);


      Vector delta_P( total_length + 1, 0.0 );
      solveLinearSystem(total_length + 1, LHS_matrix_P.valuesModifiable(),
                        RHS_vector_P.valuesModifiable(), delta_P.valuesModifiable());


      Vector delta_x_P(total_length, 0.0);
      Vector delta_z_P(total_length, 0.0);
      for (int i = 0; i < total_length; i++){
        delta_x_P.set(i, delta_P.values()[i]);
        delta_z_P.set(i, (r_mu.values()[i] - vector_z_.values()[i]*delta_x_P.values()[i])/vector_x_.values()[i] );
      }
      double delta_y_P = delta_P.values()[total_length];


      double alpha_x_P, alpha_y_P, alpha_z_P;
      calculate_step_sizes(delta_x_P, delta_y_P, delta_z_P, r_pri_, r_dual_, vector_x_, vector_z_, vector_c_, matrix_Q_, matrix_At_, scalar_y_, scalar_b_, beta_, alpha_x_P, alpha_y_P, alpha_z_P);


      Vector vector_x_new(total_length, 0.0);
      vector_x_new.copy(vector_x_);
      vector_x_new.addScaledVector(alpha_x_P, delta_x_P);

      Vector vector_z_new(total_length, 0.0);
      vector_z_new.copy(vector_z_);
      vector_z_new.addScaledVector(alpha_z_P, delta_z_P);

      double mu_aff = mu_aff_factor_ * vector_x_new.innerProduct(vector_z_new) / total_length;
      sigma_ *= 0;
      sigma_ = fmax( 1e-10 / mu_ , pow( (mu_aff / mu_) , 3) );
      


      // ******* CORRECTOR STEP *******
      SymmetricMatrixDense LHS_matrix;
      LHS_matrix.setAsDiagonal(total_length + 1, 0);
      for (int i = 0; i < total_length ; i++){
        for (int j = i; j < total_length ; j++){
          LHS_matrix.setElement(i, j, -matrix_Q_.element(i, j));
        }
        LHS_matrix.setElement(i, i, LHS_matrix.element(i,i) - vector_z_.values()[i] / vector_x_.values()[i] );
      }
      for (int i = 0; i < omega_length ; i++){
        LHS_matrix.setElement(i, total_length, 1);
      }


      Vector r_mu_C(total_length, 0.0);
      Vector vector_ones_TL(total_length, 1.0);
      r_mu_C.copy( r_cent_ );  
      r_mu_C.addScaledVector(sigma_ * mu_, vector_ones_TL);

      Vector RHS_vector(total_length + 1, 0.0);
      for (int i = 0; i < total_length; i++){
        RHS_vector.set( i, r_dual_.values()[i] - (r_mu_C.values()[i])/ vector_x_.values()[i]);
      }
      RHS_vector.set(total_length, r_pri_);

      Vector delta( total_length + 1, 0.0);
      solveLinearSystem(total_length + 1, LHS_matrix.valuesModifiable(),
                        RHS_vector.valuesModifiable(), delta.valuesModifiable());


      Vector delta_x(total_length, 0.0);
      Vector delta_z(total_length, 0.0);
      for (int i = 0; i < total_length; i++){
        delta_x.set(i, delta.values()[i]);
        delta_z.set(i, (r_mu_C.values()[i] - vector_z_.values()[i]*delta_x.values()[i])/vector_x_.values()[i] );
      }
      double delta_y = delta.values()[total_length];



      double alpha_x, alpha_y, alpha_z;
      calculate_step_sizes(delta_x, delta_y, delta_z, r_pri_, r_dual_, vector_x_, vector_z_, vector_c_, matrix_Q_, matrix_At_, scalar_y_, scalar_b_, beta_, alpha_x, alpha_y, alpha_z);


      
      // // vector_x_.print(reporter, "vector_x_before");

      // std::cout << "Beg of vector x.\n" << std::endl;

      // // // Print each primal_solution value in the loop
      // for (int i = 0; i < total_length; i++) {
      //     std::cout << "x[" << i << ", " << "] = " << vector_x_.values()[i] << std::endl;
      // }


      // Update the elements
      vector_x_.addScaledVector(alpha_x, delta_x);
      scalar_y_ += alpha_y*delta_y;
      vector_z_.addScaledVector(alpha_z, delta_z);

      
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
              "%-20d%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e%-20.6e\n",
              iter_count_,
              vector_x_.min(),
              vector_z_.min(),
              mu_,
              abs(r_pri_),
              r_dual_.normInf(),
              r_cent_.normInf(),
              delta_x.normInf(),
              abs(delta_y),
              delta_z.normInf(),
              alpha_x,
              alpha_y,
              alpha_z
          );
        

      mu_=fmax( vector_x_.innerProduct(vector_z_) / total_length, 1e-10 );

    }//end while

    // std::cout << "Accumulated time for Part 1: " << part1_time << " seconds" << std::endl;
    // std::cout << "Accumulated time for Part 2: " << part2_time << " seconds" << std::endl;
    // std::cout << "Accumulated time for Part 3: " << part3_time << " seconds" << std::endl;
    // std::cout << "Accumulated time for Part 4: " << part4_time << " seconds" << std::endl;

    if (kkt_error_ >= -tol_out_) {  // why -tol?
      THROW_EXCEPTION(QP_SUCCESS_EXCEPTION, "QP solve successful.");
    }


    if (iter_count_ >= max_iter_out_) { 
      THROW_EXCEPTION(QP_ITERATION_LIMIT_EXCEPTION, "QP solve unsuccessful. Iteration limit reached.");
    }
      
    
    if ((clock() - quantities->startTime()) / (double)CLOCKS_PER_SEC >= quantities->cpuTimeLimit()) {
      THROW_EXCEPTION(QP_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
    }


    // ****** END for Lara ******
    // functions added/modified: finalizeSolution, calculate_step_sizes, and solveLinearSystem.


    // ********** TO BE ERASED FROM HERE ************

    // Initialize minimum index and value
    int index = -1;
    double value = NONOPT_DOUBLE_INFINITY;

    // Loop through vector list
    for (int i = 0; i < (int)vector_list_.size(); i++) {

      // Compute objective value 0.5*g_i^T*W*g_i-b_i
      double t = 0.5 * matrix_->innerProductOfInverse(*vector_list_[i]) - vector_[i];

      // Check for minimum
      if (i == 0 || t < value) {

        // Update minimum's index and value
        index = i;
        value = t;

      } // end if

    } // end for

    // Check for failed index choice
    if (index == -1) {
      THROW_EXCEPTION(QP_INPUT_ERROR_EXCEPTION, "QP solve unsuccessful. Input error.");
    }

    // Update positive set
    omega_positive_.push_back(index);

    // Set factor
    factor_[0] = sqrt(1.0 + matrix_->innerProductOfInverse(*vector_list_[index]));

    // Set dual multiplier b_i-g_i^T*W*g_i
    multiplier_ = 1.0 + vector_[index] - pow(factor_[0], 2);

    // Set system solution
    system_solution_[0] = 1.0;
    system_solution_best_[0] = 1.0;
    inner_solution_1_[0] = 1.0 / factor_[0];
    inner_solution_2_[0] = vector_[index] / factor_[0];

  
    // Print message
    reporter->printf(R_QP, R_PER_ITERATION, "\n");
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "Entering main iteration loop\n");

    // Set iteration limit
    int iteration_limit = fmax(iteration_limit_minimum_, fmin(2 * ((int)vector_list_.size() + gamma_length_), iteration_limit_maximum_));

    // Iteration loop
    while (true) {

      // Print header information
      if (iteration_count_ == 0) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "======================================================================================================\n");
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "Starting iteration %6d and Inner iteration %6d\n", quantities->iterationCounter(), quantities->innerIterationCounter());
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "======================================================================================================\n");
      } // end if

      // Print message
      if (iteration_count_ % 20 == 0) {
        reporter->printf(R_QP, R_PER_ITERATION, "=======================================================\n"
                                                "  Iter.    |S|     |P|     |N|    min(KKT)  Set changes\n"
                                                "=======================================================\n");
      }
      reporter->printf(R_QP, R_PER_ITERATION, " %6d  %6d  %6d  %6d", iteration_count_, (int)omega_positive_.size(), (int)gamma_positive_.size(), (int)gamma_negative_.size());
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "=======================================\n"
                                                    "Starting iteration %8d of %8d\n"
                                                    "=======================================\n",
                       iteration_count_,
                       iteration_limit);
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "omega_positive (%6d elements):", (int)omega_positive_.size());
      for (int i = 0; i < (int)omega_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_positive (%6d elements):", (int)gamma_positive_.size());
      for (int i = 0; i < (int)gamma_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_negative (%6d elements):", (int)gamma_negative_.size());
      for (int i = 0; i < (int)gamma_negative_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

      // Update best solution
      bool real_solution = updateBestSolution();

      // Evaluate primal vectors
      evaluatePrimalVectors();

      // Check nan error
      if (!real_solution) {
        THROW_EXCEPTION(QP_NAN_ERROR_EXCEPTION, "QP solve unsuccessful.  NaN error.");
      }

      // Initialize solution and KKT vectors
      std::vector<double> kkt_residual_omega((int)vector_.size(), 0.0);
      std::vector<double> kkt_residual_gamma_positive(gamma_length_, 0.0);
      std::vector<double> kkt_residual_gamma_negative(gamma_length_, 0.0);

      // Declare KKT minimum element set
      int kkt_residual_minimum_set;
      int kkt_residual_minimum_index;

      // Evaluate omega's KKT error components  -g_i^T*W*g_i-g_i^Td
      for (int i = 0; i < (int)vector_.size(); i++) {
        kkt_residual_omega[i] = multiplier_ - vector_[i] - vector_list_[i]->innerProduct(primal_solution_);
      }

      // Zero-out omega's KKT error components for positive set
      for (int i = 0; i < (int)omega_positive_.size(); i++) {
        kkt_residual_omega[omega_positive_[i]] = 0.0;
      }

      // Evaluate gamma's KKT error components for positive side
      for (int i = 0; i < gamma_length_; i++) {
        kkt_residual_gamma_positive[i] = scalar_ - primal_solution_.values()[i];
      }

      // Zero-out gamma's KKT error components for positive side and positive set
      for (int i = 0; i < (int)gamma_positive_.size(); i++) {
        kkt_residual_gamma_positive[gamma_positive_[i]] = 0.0;
      }

      // Evaluate gamma's KKT error components for negative side
      for (int i = 0; i < gamma_length_; i++) {
        kkt_residual_gamma_negative[i] = scalar_ + primal_solution_.values()[i];
      }

      // Zero-out gamma's KKT error components for negative side and negative set
      for (int i = 0; i < (int)gamma_negative_.size(); i++) {
        kkt_residual_gamma_negative[gamma_negative_[i]] = 0.0;
      }

      // Determine minimum element indices
      int kkt_residual_omega_minimum_index = distance(kkt_residual_omega.begin(), min_element(kkt_residual_omega.begin(), kkt_residual_omega.end()));
      int kkt_residual_gamma_positive_minimum_index = distance(kkt_residual_gamma_positive.begin(), min_element(kkt_residual_gamma_positive.begin(), kkt_residual_gamma_positive.end()));
      int kkt_residual_gamma_negative_minimum_index = distance(kkt_residual_gamma_negative.begin(), min_element(kkt_residual_gamma_negative.begin(), kkt_residual_gamma_negative.end()));

      // Set minimum elements
      double kkt_residual_omega_minimum = kkt_residual_omega[kkt_residual_omega_minimum_index];
      double kkt_residual_gamma_positive_minimum = kkt_residual_gamma_positive[kkt_residual_gamma_positive_minimum_index];
      double kkt_residual_gamma_negative_minimum = kkt_residual_gamma_negative[kkt_residual_gamma_negative_minimum_index];

      // Determine value and type of most negative KKT component
      if (kkt_residual_omega_minimum < -kkt_tolerance_) {
        kkt_error_ = kkt_residual_omega_minimum;
        kkt_residual_minimum_set = 1;
        kkt_residual_minimum_index = kkt_residual_omega_minimum_index;
      } // end if
      else {
        if (kkt_residual_omega_minimum <= kkt_residual_gamma_positive_minimum) {
          if (kkt_residual_omega_minimum <= kkt_residual_gamma_negative_minimum) {
            kkt_error_ = kkt_residual_omega_minimum;
            kkt_residual_minimum_set = 1;
            kkt_residual_minimum_index = kkt_residual_omega_minimum_index;
          } // end if
          else {
            kkt_error_ = kkt_residual_gamma_negative_minimum;
            kkt_residual_minimum_set = 3;
            kkt_residual_minimum_index = kkt_residual_gamma_negative_minimum_index;
          } // end else
        }   // end if
        else {
          if (kkt_residual_gamma_positive_minimum <= kkt_residual_gamma_negative_minimum) {
            kkt_error_ = kkt_residual_gamma_positive_minimum;
            kkt_residual_minimum_set = 2;
            kkt_residual_minimum_index = kkt_residual_gamma_positive_minimum_index;
          } // end if
          else {
            kkt_error_ = kkt_residual_gamma_negative_minimum;
            kkt_residual_minimum_set = 3;
            kkt_residual_minimum_index = kkt_residual_gamma_negative_minimum_index;
          } // end else
        }   // end else
      }     // end else

      // Print message
      reporter->printf(R_QP, R_PER_ITERATION, "  %+.2e", kkt_error_);
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Set of minimum KKT element is %d\n"
                                                    "Index of minimum KKT element is %d with value %+.16e\n",
                       kkt_residual_minimum_set,
                       kkt_residual_minimum_index,
                       kkt_error_);

      // Increment iteration counter
      iteration_count_++;

      // Check for successful solve
      if (kkt_error_ >= -kkt_tolerance_) {
        THROW_EXCEPTION(QP_SUCCESS_EXCEPTION, "QP solve successful.");
      }

      // Check for inexact termination
      if (allow_inexact_termination_ &&
          iteration_count_ >= (int)ceil(inexact_termination_initialization_factor_ * (double)vector_.size()) &&
          (iteration_count_ - (int)ceil(inexact_termination_initialization_factor_ * (double)vector_.size())) % inexact_termination_check_interval_ == 0 &&
          inexactTerminationCondition(quantities, reporter)) {
        THROW_EXCEPTION(QP_SUCCESS_EXCEPTION, "QP solve successful.");
      }

      // Check for iteration limit ** use for outer **
      if (iteration_count_ >= iteration_limit) { 
        THROW_EXCEPTION(QP_ITERATION_LIMIT_EXCEPTION, "QP solve unsuccessful. Iteration limit reached.");
      }

      // Check for CPU time limit
      if ((clock() - quantities->startTime()) / (double)CLOCKS_PER_SEC >= quantities->cpuTimeLimit()) {
        THROW_EXCEPTION(QP_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
      }

      // Print message
      reporter->printf(R_QP, R_PER_ITERATION, "  %d", kkt_residual_minimum_set);
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be added to index set %d\n", kkt_residual_minimum_index, kkt_residual_minimum_set);

      // Evaluate new system vector
      evaluateSystemVector(kkt_residual_minimum_set, kkt_residual_minimum_index, new_system_vector_);

      // Print message
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Solving intermediate system for least squares\n");

      // Solve intermediate system for least squares
      solveSystemTranspose(new_system_vector_, inner_solution_3_);

      // Declare new diagonal value
      double new_diagonal_squared;

      // Check which set is being updated
      if (kkt_residual_minimum_set == 1) {
        new_diagonal_squared = 1.0 + matrix_->innerProductOfInverse(*vector_list_[kkt_residual_minimum_index]);
      }
      else {
        new_diagonal_squared = matrix_->elementOfInverse(kkt_residual_minimum_index, kkt_residual_minimum_index);
      }

      // Set inputs for BLASLAPACK
      int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
      int increment1 = 1;

      // Compute new diagonal for factor (squared)
      double rho2 = fmax(0.0, new_diagonal_squared - ddot_(&length, inner_solution_3_, &increment1, inner_solution_3_, &increment1));

      // Compute comparison value of new diagonal for factor (squared)
      double rhoT = cholesky_tolerance_ * new_diagonal_squared;

      // Initialize linear independence check boolean
      bool linear_independence_flag = false;

      // Initialize value to add in augmentation
      double augmentation_value = 0.0;

      // Check for sufficiently large new diagonal for factor (squared)
      if (rho2 > rhoT) {
        linear_independence_flag = true;
      }
      else {

        // Print message
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "Solving least squares system\n");

        // Solve intermediate system for least squares
        solveSystem(inner_solution_3_, inner_solution_ls_);

        // Set inputs for BLASLAPACK
        int length = (int)omega_positive_.size();
        int increment0 = 0;
        int increment1 = 1;
        double value = 1.0;

        // Declare residual value
        double residual = ddot_(&length, &value, &increment0, inner_solution_ls_, &increment1);

        // Check for sufficiently negative value
        if (residual - 1.0 < -linear_independence_tolerance_) {
          linear_independence_flag = true;
        }

      } // end else

      // Print message
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Linear independence check yields %d\n", linear_independence_flag);

      // Check for column exchange
      if (linear_independence_flag == false) {

        // Print message
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "Exchange! Looking for index to delete\n");

        // Initialize minimum index and values
        int delete_index = -1;
        double delete_value = NONOPT_DOUBLE_INFINITY;
        int delete_set = -1;

        // Compute index with minimum value
        if ((int)omega_positive_.size() > 1 || (int)gamma_positive_.size() + (int)gamma_negative_.size() == 0) {
          for (int i = 0; i < (int)omega_positive_.size(); i++) {
            if (inner_solution_ls_[i] > 0.0) {
              double temporary_scalar = system_solution_[i] / inner_solution_ls_[i];
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 1;
              } // end if
            }   // end if
          }     // end for
        }       // end if
        else {
          for (int i = 0; i < (int)gamma_positive_.size(); i++) {
            if (inner_solution_ls_[(int)omega_positive_.size() + i] > 0.0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + i] / inner_solution_ls_[(int)omega_positive_.size() + i];
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 2;
              } // end if
            }   // end if
          }     // end for
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            if (inner_solution_ls_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] > 0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] / inner_solution_ls_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i];
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 3;
              } // end if
            }   // end if
          }     // end for
        }       // end else

        // Update augmentation value
        augmentation_value = delete_value;

        // Print message
        reporter->printf(R_QP, R_PER_ITERATION, "  -%d", delete_set);
        if (delete_set == 1) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be deleted from omega's positive set\n", omega_positive_[delete_index]);
        }
        else if (delete_set == 2) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be deleted from gamma's positive set\n", gamma_positive_[delete_index]);
        }
        else if (delete_set == 3) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be deleted from gamma's negative set\n", gamma_negative_[delete_index]);
        }
        else {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "Uh oh!  Search for element to delete failed!\n");
        }

        // Check if search for element to delete was successful
        if (delete_set > 0) {

          // Set inputs for BLASLAPACK
          int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
          double value = -delete_value;
          int increment = 1;

          // Update solution
          daxpy_(&length, &value, inner_solution_ls_, &increment, system_solution_, &increment);

          // Perform set deletion
          setDelete(reporter, delete_set, delete_index, inner_solution_1_, inner_solution_2_);

        } // end if

        // Print sets
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "omega_positive (%6d elements):", (int)omega_positive_.size());
        for (int i = 0; i < (int)omega_positive_.size(); i++) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
        }
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_positive (%6d elements):", (int)gamma_positive_.size());
        for (int i = 0; i < (int)gamma_positive_.size(); i++) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
        }
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_negative (%6d elements):", (int)gamma_negative_.size());
        for (int i = 0; i < (int)gamma_negative_.size(); i++) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
        }
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

      } // end if (linear_independence_flag == false)

      // Print message
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Performing set augmentation\n");

      // Perform set augmentation
      bool augment_success = setAugment(reporter, kkt_residual_minimum_set, kkt_residual_minimum_index, new_system_vector_, inner_solution_1_, inner_solution_2_, augmentation_value);

      // Print sets
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "omega_positive (%6d elements):", (int)omega_positive_.size());
      for (int i = 0; i < (int)omega_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_positive (%6d elements):", (int)gamma_positive_.size());
      for (int i = 0; i < (int)gamma_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_negative (%6d elements):", (int)gamma_negative_.size());
      for (int i = 0; i < (int)gamma_negative_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

      // Check for Cholesky update error
      if (!augment_success) {

        // Check for failure on factorization error
        if (fail_on_factorization_error_) {
          THROW_EXCEPTION(QP_FACTORIZATION_ERROR_EXCEPTION, "QP solve unsuccessful. Factorization error.");
        }

        // Print message
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "Re-doing Cholesky factorization from scratch!\n");

        // Try to compute from scratch
        choleskyFromScratch(reporter);

      } // end if

      // Compute multiplier
      evaluatePrimalMultiplier(inner_solution_1_, inner_solution_2_);

      // Set inner iteration limit
      int inner_iteration_limit = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();

      // Subproblem solution loop
      for (int inner_iteration_count = 0; inner_iteration_count < inner_iteration_limit; inner_iteration_count++) {

        // Print message
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "Solving subproblem\n");

        // Set inputs for BLASLAPACK
        int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
        int increment = 1;
        double value = 1.0 - multiplier_;

        // Set right-hand side vector
        dcopy_(&length, inner_solution_2_, &increment, right_hand_side_, &increment);
        daxpy_(&length, &value, inner_solution_1_, &increment, right_hand_side_, &increment);

        // Solve subproblem
        solveSystem(right_hand_side_, inner_solution_trial_);

        // Declare boolean for correct signs
        bool correct_signs = true;

        // Check signs of subproblem solution elements
        if (correct_signs) {
          for (int i = 0; i < (int)omega_positive_.size() + (int)gamma_positive_.size(); i++) {
            if (inner_solution_trial_[i] <= 0.0) {
              correct_signs = false;
              break;
            } // end if
          }   // end for
        }     // end if
        if (correct_signs) {
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            if (inner_solution_trial_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] >= 0.0) {
              correct_signs = false;
              break;
            } // end if
          }   // end for
        }     // end if

        // Check feasibility of subproblem solution
        if (correct_signs) {

          // Print message
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "Feasible subproblem solution! Breaking loop\n");

          // Replace current solution
          dcopy_(&length, inner_solution_trial_, &increment, system_solution_, &increment);

          // Break subproblem solution loop
          break;

        } // end if
        else {

          // Print message
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "Infeasible subproblem solution! Looking for index to delete\n");

          // Initialize minimum index and values
          int delete_index = -1;
          double delete_value = NONOPT_DOUBLE_INFINITY;
          int delete_set = -1;

          // Compute index with minimum value
          for (int i = 0; i < (int)omega_positive_.size(); i++) {
            if (inner_solution_trial_[i] < 0.0) {
              double temporary_scalar = system_solution_[i] / (system_solution_[i] - inner_solution_trial_[i]);
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 1;
              } // end if
            }   // end if
          }     // end for
          for (int i = 0; i < (int)gamma_positive_.size(); i++) {
            if (inner_solution_trial_[(int)omega_positive_.size() + i] < 0.0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + i] / (system_solution_[(int)omega_positive_.size() + i] - inner_solution_trial_[(int)omega_positive_.size() + i]);
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 2;
              } // end if
            }   // end if
          }     // end for
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            if (inner_solution_trial_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] > 0.0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] / (system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] - inner_solution_trial_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i]);
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 3;
              } // end if
            }   // end if
          }     // end for

          // Print message
          reporter->printf(R_QP, R_PER_ITERATION, "  %+d", -delete_set);
          if (delete_set == 1) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be deleted from omega's positive set\n", omega_positive_[delete_index]);
          }
          else if (delete_set == 2) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be deleted from gamma's positive set\n", gamma_positive_[delete_index]);
          }
          else if (delete_set == 3) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, "Index %d will be deleted from gamma's negative set\n", gamma_negative_[delete_index]);
          }
          else {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, "Uh oh!  Search for element to delete failed!\n");
          }

          // Check if search for element to delete was successful
          if (delete_set > 0) {

            // Update value
            delete_value = fmin(1.0, delete_value);

            // Set inputs for BLASLAPACK
            int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
            double value = 1.0 - delete_value;
            int increment = 1;

            // Scale solution
            dscal_(&length, &value, system_solution_, &increment);

            // Update solution
            daxpy_(&length, &delete_value, inner_solution_trial_, &increment, system_solution_, &increment);

            // Perform set deletion
            setDelete(reporter, delete_set, delete_index, inner_solution_1_, inner_solution_2_);

            // Evaluate multiplier
            evaluatePrimalMultiplier(inner_solution_1_, inner_solution_2_);

          } // end if

          // Print sets
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "omega_positive (%6d elements):", (int)omega_positive_.size());
          for (int i = 0; i < (int)omega_positive_.size(); i++) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
          }
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_positive (%6d elements):", (int)gamma_positive_.size());
          for (int i = 0; i < (int)gamma_positive_.size(); i++) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
          }
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "gamma_negative (%6d elements):", (int)gamma_negative_.size());
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
          }
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

        } // end else (correct_signs)

      } // end for (subproblem solution loop)

      // Print new line
      reporter->printf(R_QP, R_PER_ITERATION, "\n");

    } // end of main iteration loop

  } // end try

  // catch exceptions
  // catch (QP_SUCCESS_EXCEPTION& exec) { //this
  //   setStatus(QP_SUCCESS);
  // } catch (QP_CPU_TIME_LIMIT_EXCEPTION& exec) { //this
  //   setStatus(QP_CPU_TIME_LIMIT);
  //   THROW_EXCEPTION(NONOPT_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
  // }  catch (QP_ITERATION_LIMIT_EXCEPTION& exec) {  //this
  //   setStatus(QP_ITERATION_LIMIT);
  // }

  catch (QP_SUCCESS_EXCEPTION& exec) { //this
    setStatus(QP_SUCCESS);
  } catch (QP_CPU_TIME_LIMIT_EXCEPTION& exec) { //this
    setStatus(QP_CPU_TIME_LIMIT);
    THROW_EXCEPTION(NONOPT_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
  } catch (QP_FACTORIZATION_ERROR_EXCEPTION& exec) {
    setStatus(QP_FACTORIZATION_ERROR);
  } catch (QP_INPUT_ERROR_EXCEPTION& exec) {
    setStatus(QP_INPUT_ERROR);
  } catch (QP_ITERATION_LIMIT_EXCEPTION& exec) {  //this
    setStatus(QP_ITERATION_LIMIT);
  } catch (QP_NAN_ERROR_EXCEPTION& exec) {
    setStatus(QP_NAN_ERROR);
  }





  // Print new line
  reporter->printf(R_QP, R_PER_ITERATION, "\n");

  // Print finalizing
  reporter->printf(R_QP, R_PER_INNER_ITERATION, "Finalizing solution\n");

  // Finalize solution
  finalizeSolution(reporter); // updated

  

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

// Update best solution
bool QPSolverInteriorPoint::updateBestSolution()
{

  // Initialize boolean
  bool real_solution = true;

  // Loop through solution
  for (int i = 0; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size(); i++) {
    if (std::isnan(system_solution_[i])) {
      real_solution = false;
    }
  }

  // Confirm real solution
  if (real_solution) {

    // Update best sets
    omega_positive_best_ = omega_positive_;
    gamma_positive_best_ = gamma_positive_;
    gamma_negative_best_ = gamma_negative_;

    // Update solution
    for (int i = 0; i < (int)omega_positive_best_.size() + (int)gamma_positive_best_.size() + (int)gamma_negative_best_.size(); i++) {
      system_solution_best_[i] = system_solution_[i];
    }

  } // end if

  // Return
  return real_solution;

} // end updateBestSolution

// Cholesky augmentation
bool QPSolverInteriorPoint::choleskyAugment(double system_vector[],
                                            int index,
                                            double solution1[],
                                            double value1,
                                            double solution2[],
                                            double value2)
{

  // Initialize return value
  bool success_without_factorization_error = true;

  // Set new length
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();

  // Resize factor_?
  while ((length - 1) * (system_solution_length_ + 1) >= factor_length_ + 1) {

    // Declare temporary vector
    double* factor_temp = new double[factor_length_];

    // Set inputs for BLASLAPACK
    int increment = 1;

    // Copy factor_ to temporary vector
    dcopy_(&factor_length_, factor_, &increment, factor_temp, &increment);

    // Delete factor_ and resize
    delete[] factor_;
    factor_ = nullptr;
    int factor_length_new = 2 * factor_length_;
    factor_ = new double[factor_length_new];

    // Copy temporary vector back to factor_
    dcopy_(&factor_length_, factor_temp, &increment, factor_, &increment);

    // Set new factor length
    factor_length_ = factor_length_new;

  } // end while

  // "Add" zero values to R by shifting values
  for (int i = 0; i < length; i++) {
    for (int j = length - 1; j > index; j--) {
      factor_[i * system_solution_length_ + j] = factor_[i * system_solution_length_ + j - 1];
    }
    factor_[i * system_solution_length_ + index] = 0.0;
  } // end for
  for (int j = 0; j < length; j++) {
    for (int i = length - 1; i > index; i--) {
      factor_[i * system_solution_length_ + j] = factor_[(i - 1) * system_solution_length_ + j];
    }
    factor_[index * system_solution_length_ + j] = 0.0;
  } // end for

  // "Add" zero values to solutions by shifting values
  for (int i = length - 1; i > index; i--) {
    solution1[i] = solution1[i - 1];
    solution2[i] = solution2[i - 1];
  } // end for
  solution1[index] = 0.0;
  solution2[index] = 0.0;

  // Update Cholesky factor (first part)
  for (int i = 0; i < index; i++) {
    factor_[i * system_solution_length_ + index] = system_vector[i] / factor_[i * system_solution_length_ + i];
    for (int j = i + 1; j < index; j++) {
      system_vector[j] = system_vector[j] - factor_[i * system_solution_length_ + j] * factor_[i * system_solution_length_ + index];
    }
  } // end for

  // Compute temporary scalar
  double temporary_scalar = 0.0;
  for (int i = 0; i < index; i++) {
    temporary_scalar = temporary_scalar + factor_[i * system_solution_length_ + index] * factor_[i * system_solution_length_ + index];
  }

  // Check diagonal element tolerances
  double diagonal_squared = system_vector[index] - temporary_scalar;
  if (diagonal_squared >= cholesky_tolerance_) {
    factor_[index * system_solution_length_ + index] = sqrt(diagonal_squared);
  }
  else {
    factor_[index * system_solution_length_ + index] = sqrt(cholesky_tolerance_);
    success_without_factorization_error = false;
  } // end else

  // Update Cholesky factor (third part)
  for (int i = index + 1; i < length; i++) {
    temporary_scalar = 0.0;
    for (int j = 0; j < index; j++) {
      temporary_scalar = temporary_scalar + factor_[j * system_solution_length_ + i] * factor_[j * system_solution_length_ + index];
    }
    factor_[index * system_solution_length_ + i] = (system_vector[i] - temporary_scalar) / factor_[index * system_solution_length_ + index];
  } // end for

  // Update solution1
  temporary_scalar = 0.0;
  for (int i = 0; i < index; i++) {
    temporary_scalar = temporary_scalar + solution1[i] * factor_[i * system_solution_length_ + index];
  }
  solution1[index] = (value1 - temporary_scalar) / factor_[index * system_solution_length_ + index];

  // Update solution2
  temporary_scalar = 0.0;
  for (int i = 0; i < index; i++) {
    temporary_scalar = temporary_scalar + solution2[i] * factor_[i * system_solution_length_ + index];
  }
  solution2[index] = (value2 - temporary_scalar) / factor_[index * system_solution_length_ + index];

  // Declare and initialize temporary vector
  double* temporary_vector = new double[length - index + 1];
  for (int i = 0; i < length - index - 1; i++) {
    temporary_vector[i] = factor_[index * system_solution_length_ + index + i + 1];
  }
  temporary_vector[length - index - 1] = solution1[index];
  temporary_vector[length - index] = solution2[index];

  // Update Cholesky factor (last part)
  for (int i = index + 1; i < length; i++) {

    // Declare and initialize scalar
    double t_squared = factor_[i * system_solution_length_ + i] * factor_[i * system_solution_length_ + i] - temporary_vector[i - index - 1] * temporary_vector[i - index - 1];
    double t = 0.0;

    // Check diagonal element tolerances
    if (t_squared >= cholesky_tolerance_) {
      t = sqrt(t_squared);
    }
    else {
      t = sqrt(cholesky_tolerance_);
      success_without_factorization_error = false;
    }

    // Declare and set rotation values
    double c = t / factor_[i * system_solution_length_ + i];
    double s = temporary_vector[i - index - 1] / factor_[i * system_solution_length_ + i];

    // Update factor values
    factor_[i * system_solution_length_ + i] = t;
    for (int j = i + 1; j < length; j++) {
      factor_[i * system_solution_length_ + j] = (factor_[i * system_solution_length_ + j] - s * temporary_vector[j - index - 1]) / c;
    }

    // Update solution values
    solution1[i] = (solution1[i] - s * temporary_vector[length - index - 1]) / c;
    solution2[i] = (solution2[i] - s * temporary_vector[length - index]) / c;

    // Update temporary vector values
    for (int j = i + 1; j < length; j++) {
      temporary_vector[j - index - 1] = c * temporary_vector[j - index - 1] - s * factor_[i * system_solution_length_ + j];
    }
    temporary_vector[length - index - 1] = c * temporary_vector[length - index - 1] - s * solution1[i];
    temporary_vector[length - index] = c * temporary_vector[length - index] - s * solution2[i];

  } // end for

  // Delete temporary vector
  if (temporary_vector != nullptr) {
    delete[] temporary_vector;
    temporary_vector = nullptr;
  } // end if

  // Return
  return success_without_factorization_error;

} // end choleskyAugment

// Cholesky deletion
void QPSolverInteriorPoint::choleskyDelete(int index,
                                           double solution1[],
                                           double solution2[])
{

  // Set new length
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();

  // "Delete" column of R by shifting columns to the left
  for (int i = 0; i <= length; i++) {
    for (int j = index; j < length; j++) {
      factor_[i * system_solution_length_ + j] = factor_[i * system_solution_length_ + j + 1];
    }
  } // end for

  // Apply Givens rotations
  for (int i = index; i < length; i++) {

    // Compute scalars
    double p = factor_[i * system_solution_length_ + i];
    double q = factor_[(i + 1) * system_solution_length_ + i];
    double c = p / sqrt(p * p + q * q);
    double s = q / sqrt(p * p + q * q);

    // Loop over columns
    for (int j = i; j < length; j++) {

      // Perform rotation for factor values
      double y = c * factor_[i * system_solution_length_ + j] + s * factor_[(i + 1) * system_solution_length_ + j];
      double z = -s * factor_[i * system_solution_length_ + j] + c * factor_[(i + 1) * system_solution_length_ + j];

      // Set new factor values
      factor_[i * system_solution_length_ + j] = y;
      factor_[(i + 1) * system_solution_length_ + j] = z;

    } // end for

    // Set new zero value
    factor_[(i + 1) * system_solution_length_ + i] = 0.0;

    // Perform rotation for solution1 values
    double y = c * solution1[i] + s * solution1[i + 1];
    double z = -s * solution1[i] + c * solution1[i + 1];

    // Set new solution1 values
    solution1[i] = y;
    solution1[i + 1] = z;

    // Perform rotation for solution2 values
    y = c * solution2[i] + s * solution2[i + 1];
    z = -s * solution2[i] + c * solution2[i + 1];

    // Set new solution2 values
    solution2[i] = y;
    solution2[i + 1] = z;

  } // end for

} // end choleskyDelete

// Cholesky factorization from scratch
void QPSolverInteriorPoint::choleskyFromScratch(const Reporter* reporter)
{

  // Set size of matrix
  int size = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();

  // Print message
  reporter->printf(R_QP, R_PER_INNER_ITERATION, "Re-doing Cholesky of size %d\n", size);

  // Declare matrix and initialize
  double* matrix = new double[size * size];
  for (int i = 0; i < size * size; i++) {
    matrix[i] = 0.0;
  }

  // Declare WG vectors
  std::vector<std::shared_ptr<Vector>> WG;

  // Compute WG matrix
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    std::shared_ptr<Vector> product(new Vector(gamma_length_));
    matrix_->matrixVectorProductOfInverse(*vector_list_[omega_positive_[i]], *product);
    WG.push_back(product);
  } // end for

  // Compute (1,1)-block (upper triangle)
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    for (int j = i; j < (int)omega_positive_.size(); j++) {
      matrix[i * size + j] = 1.0 + vector_list_[omega_positive_[i]]->innerProduct(*WG[j]);
    }
  } // end for

  // Compute (1,2)-block
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    for (int j = 0; j < (int)gamma_positive_.size(); j++) {
      matrix[i * size + (int)omega_positive_.size() + j] = WG[i]->values()[gamma_positive_[j]];
    }
  } // end for

  // Compute (1,3)-block
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    for (int j = 0; j < (int)gamma_negative_.size(); j++) {
      matrix[i * size + (int)omega_positive_.size() + (int)gamma_positive_.size() + j] = WG[i]->values()[gamma_negative_[j]];
    }
  } // end for

  // Compute (2,2)-block (upper triangle)
  for (int i = 0; i < (int)gamma_positive_.size(); i++) {
    for (int j = i; j < (int)gamma_positive_.size(); j++) {
      matrix[((int)omega_positive_.size() + i) * size + ((int)omega_positive_.size() + j)] = matrix_->elementOfInverse(gamma_positive_[i], gamma_positive_[j]);
    }
  } // end for

  // Compute (2,3)-block (upper triangle)
  for (int i = 0; i < (int)gamma_positive_.size(); i++) {
    for (int j = 0; j < (int)gamma_negative_.size(); j++) {
      matrix[((int)omega_positive_.size() + i) * size + ((int)omega_positive_.size() + (int)gamma_positive_.size() + j)] = matrix_->elementOfInverse(gamma_positive_[i], gamma_negative_[j]);
    }
  } // end for

  // Compute (3,3)-block (upper triangle)
  for (int i = 0; i < (int)gamma_negative_.size(); i++) {
    for (int j = i; j < (int)gamma_negative_.size(); j++) {
      matrix[((int)omega_positive_.size() + (int)gamma_positive_.size() + i) * size + ((int)omega_positive_.size() + (int)gamma_positive_.size() + j)] = matrix_->elementOfInverse(gamma_negative_[i], gamma_negative_[j]);
    }
  } // end for

  // Set inputs for BLASLAPACK
  char upper_lower = 'L'; // Use "lower" since Fortran uses column-major ordering
  int flag = 0;

  // Compute factorization
  dpotf2_(&upper_lower, &size, matrix, &size, &flag);

  // Check tolerance on diagonal
  for (int i = 0; i < size; i++) {
    if (matrix[i * size + i] < cholesky_tolerance_) {
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "ARGH! Replacing Cholesky diagonal of %+23.16e with %+23.16e.\n", matrix[i * size + i], cholesky_tolerance_);
      matrix[i * size + i] = cholesky_tolerance_;
    } // end if
  }   // end for

  // Reset factorization
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      factor_[i * system_solution_length_ + j] = matrix[i * size + j];
    }
  } // end for

  // Set right-hand side 1
  double* right_hand_side = new double[size];
  for (int i = 0; i < size; i++) {
    right_hand_side[i] = 0.0;
  }
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    right_hand_side[i] = 1.0;
  }

  // Re-solve system 1
  solveSystemTranspose(right_hand_side, inner_solution_1_);

  // Set right-hand side 2
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    right_hand_side[i] = vector_[omega_positive_[i]];
  }
  for (int i = 0; i < (int)gamma_positive_.size(); i++) {
    right_hand_side[(int)omega_positive_.size() + i] = -scalar_;
  }
  for (int i = 0; i < (int)gamma_negative_.size(); i++) {
    right_hand_side[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] = scalar_;
  }

  // Re-solve system 2
  solveSystemTranspose(right_hand_side, inner_solution_2_);

  // Delete temporary matrices
  if (matrix != nullptr) {
    delete[] matrix;
    matrix = nullptr;
  } // end if
  if (right_hand_side != nullptr) {
    delete[] right_hand_side;
    right_hand_side = nullptr;
  } // end if

} // end choleskyFromScratch

// Evaluate dual vectors
void QPSolverInteriorPoint::evaluatePrimalVectors()
{

  // Zero-out vectors
  combination_.scale(0.0);
  combination_translated_.scale(0.0);
  primal_solution_.scale(0.0);
  primal_solution_feasible_.scale(0.0);

  // Loop to compute gradient combination
  for (int i = 0; i < (int)omega_positive_best_.size(); i++) {
    combination_.addScaledVector(system_solution_best_[i], *vector_list_[omega_positive_best_[i]]);
  }

  // Initialize gradient combination shifted
  combination_translated_.copy(combination_);

  // Loops to set gradient combination shifted values
  for (int i = 0; i < (int)gamma_positive_best_.size(); i++) {
    combination_translated_.set(gamma_positive_best_[i], combination_.values()[gamma_positive_best_[i]] + system_solution_best_[(int)omega_positive_best_.size() + i]);
  }
  for (int i = 0; i < (int)gamma_negative_best_.size(); i++) {
    combination_translated_.set(gamma_negative_best_[i], combination_.values()[gamma_negative_best_[i]] + system_solution_best_[(int)omega_positive_best_.size() + (int)gamma_positive_best_.size() + i]);
  }

  // Compute matrix-vector product
  matrix_->matrixVectorProductOfInverse(combination_translated_, primal_solution_);
  primal_solution_.scale(-1.0);

  // Compute feasible primal step by projection
  primal_solution_feasible_.copy(primal_solution_);
  primal_solution_projection_scalar_ = fmin(1.0, scalar_ / primal_solution_.normInf());
  primal_solution_feasible_.scale(primal_solution_projection_scalar_);

} // end evaluatePrimalVectors

// Evaluate primal multiplier
void QPSolverInteriorPoint::evaluatePrimalMultiplier(double solution1[],
                                                     double solution2[])
{

  // Set inputs for BLASLAPACK
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  int increment = 1;

  // Evaluate inner products
  double solution1_norm_squared = ddot_(&length, solution1, &increment, solution1, &increment);
  double solution1_solution2 = ddot_(&length, solution1, &increment, solution2, &increment);

  // Evaluate multiplier
  multiplier_ = (solution1_norm_squared + solution1_solution2 - 1.0) / solution1_norm_squared;

} // end evaluatePrimalMultiplier

// Evaluate system vector
void QPSolverInteriorPoint::evaluateSystemVector(int set,
                                                 int index,
                                                 double system_vector[])
{

  // Set inputs for BLASLAPACK
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  double value = 0.0;
  int increment0 = 0;
  int increment1 = 1;

  // Initialize values
  dcopy_(&length, &value, &increment0, system_vector, &increment1);

  // Check set
  if (set == 1) {

    // Declare temporary vector
    Vector temporary_vector(gamma_length_);

    // Evaluate temporary vector
    matrix_->matrixVectorProductOfInverse(*vector_list_[index], temporary_vector);

    // Set "omega" values, i.e.,
    for (int i = 0; i < (int)omega_positive_.size(); i++) {
      system_vector[i] = system_vector[i] + ddot_(&gamma_length_, vector_list_[omega_positive_[i]]->values(), &increment1, temporary_vector.values(), &increment1) + 1.0;
    }

    // Set "gamma positive" values, i.e.,
    for (int i = 0; i < (int)gamma_positive_.size(); i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(gamma_positive_[i], col);
      system_vector[(int)omega_positive_.size() + i] = vector_list_[index]->innerProduct(col);
    } // end for

    // Set "gamma negative" values, i.e.,
    for (int i = 0; i < (int)gamma_negative_.size(); i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(gamma_negative_[i], col);
      system_vector[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] = vector_list_[index]->innerProduct(col);
    } // end for

  } // end if

  else {

    // Set "omega" values, i.e., ith value = G(:,omega_positive_[i])'*W(:,kkt_residual_minimum_index)
    for (int i = 0; i < (int)omega_positive_.size(); i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(index, col);
      system_vector[i] = vector_list_[omega_positive_[i]]->innerProduct(col);
    } // end for

    // Set "gamma positive" values, i.e., W(gamma_positive_[i],kkt_residual_minimum_index)
    for (int i = 0; i < (int)gamma_positive_.size(); i++) {
      system_vector[(int)omega_positive_.size() + i] = matrix_->elementOfInverse(gamma_positive_[i], index);
    }

    // Set "gamma negative" values, i.e., W(gamma_positive_[i],kkt_residual_minimum_index)
    for (int i = 0; i < (int)gamma_negative_.size(); i++) {
      system_vector[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] = matrix_->elementOfInverse(gamma_negative_[i], index);
    }

  } // end else

} // end evaluateSystemVector

// Finalize solution
void QPSolverInteriorPoint::finalizeSolution(const Reporter* reporter)
{

  int omega_length = (int)vector_.size();

  gamma_.setLength(gamma_length_);
  gamma_.scale(0.0);
  for (int i = 0; i < gamma_length_; i++){
    gamma_.set(i, vector_x_.values()[i + omega_length + gamma_length_] - vector_x_.values()[i + omega_length]);
  }


  Vector Gomega_new;
  Gomega_new.setLength(gamma_length_);
  Gomega_new.scale(0.0);
  for (int i = 0; i < gamma_length_; i++) {
    for (int j = 0; j < omega_length; j++){
    Gomega_new.set(i, Gomega_new.values()[i] += vector_x_.values()[j] * vector_list_[j]->values()[i] );}
  }

  primal_solution_.setLength(gamma_length_);
  primal_solution_.scale(0.0);
  primal_solution_.addScaledVector(-1.0, Gomega_new);

  if (scalar_ != NONOPT_DOUBLE_INFINITY){
    primal_solution_.addScaledVector(1.0, gamma_);
  }

  // std::cout << "Beg of Primal Solution values ( -W ( Gomega+Gamma ) ).\n" << std::endl;
  // // Print each primal_solution value in the loop
  // for (int i = 0; i < gamma_length_; i++) {
  //     std::cout << "PM(" << i << ", " << ") = " << primal_solution_.values()[i] << std::endl;
  // }
  // std::cout << "End of Primal Solution values ( -W ( Gomega+Gamma ) ).\n" << std::endl;

} // end finalizeSolution

// Resize system solution
void QPSolverInteriorPoint::resizeSystemSolution()
{

  // Declare temp vectors
  double* inner_solution_1_temp = new double[system_solution_length_];
  double* inner_solution_2_temp = new double[system_solution_length_];
  double* inner_solution_3_temp = new double[system_solution_length_];
  double* inner_solution_ls_temp = new double[system_solution_length_];
  double* inner_solution_trial_temp = new double[system_solution_length_];
  double* new_system_vector_temp = new double[system_solution_length_];
  double* right_hand_side_temp = new double[system_solution_length_];
  double* system_solution_temp = new double[system_solution_length_];
  double* system_solution_best_temp = new double[system_solution_length_];

  // Set inputs for BLASLAPACK
  int increment = 1;

  // Initialize values
  dcopy_(&system_solution_length_, inner_solution_1_, &increment, inner_solution_1_temp, &increment);
  dcopy_(&system_solution_length_, inner_solution_2_, &increment, inner_solution_2_temp, &increment);
  dcopy_(&system_solution_length_, inner_solution_3_, &increment, inner_solution_3_temp, &increment);
  dcopy_(&system_solution_length_, inner_solution_ls_, &increment, inner_solution_ls_temp, &increment);
  dcopy_(&system_solution_length_, inner_solution_trial_, &increment, inner_solution_trial_temp, &increment);
  dcopy_(&system_solution_length_, new_system_vector_, &increment, new_system_vector_temp, &increment);
  dcopy_(&system_solution_length_, right_hand_side_, &increment, right_hand_side_temp, &increment);
  dcopy_(&system_solution_length_, system_solution_, &increment, system_solution_temp, &increment);
  dcopy_(&system_solution_length_, system_solution_best_, &increment, system_solution_best_temp, &increment);

  // Delete arrays (in case they exist)
  if (inner_solution_1_ != nullptr) {
    delete[] inner_solution_1_;
    inner_solution_1_ = nullptr;
  } // end if
  if (inner_solution_2_ != nullptr) {
    delete[] inner_solution_2_;
    inner_solution_2_ = nullptr;
  } // end if
  if (inner_solution_3_ != nullptr) {
    delete[] inner_solution_3_;
    inner_solution_3_ = nullptr;
  } // end if
  if (inner_solution_ls_ != nullptr) {
    delete[] inner_solution_ls_;
    inner_solution_ls_ = nullptr;
  } // end if
  if (inner_solution_trial_ != nullptr) {
    delete[] inner_solution_trial_;
    inner_solution_trial_ = nullptr;
  } // end if
  if (new_system_vector_ != nullptr) {
    delete[] new_system_vector_;
    new_system_vector_ = nullptr;
  } // end if
  if (right_hand_side_ != nullptr) {
    delete[] right_hand_side_;
    right_hand_side_ = nullptr;
  } // end if
  if (system_solution_ != nullptr) {
    delete[] system_solution_;
    system_solution_ = nullptr;
  } // end if
  if (system_solution_best_ != nullptr) {
    delete[] system_solution_best_;
    system_solution_best_ = nullptr;
  } // end if

  // Set new system solution length
  int system_solution_length_new = 2 * system_solution_length_;

  // Declare temp vectors
  inner_solution_1_ = new double[system_solution_length_new];
  inner_solution_2_ = new double[system_solution_length_new];
  inner_solution_3_ = new double[system_solution_length_new];
  inner_solution_ls_ = new double[system_solution_length_new];
  inner_solution_trial_ = new double[system_solution_length_new];
  new_system_vector_ = new double[system_solution_length_new];
  right_hand_side_ = new double[system_solution_length_new];
  system_solution_ = new double[system_solution_length_new];
  system_solution_best_ = new double[system_solution_length_new];

  // Initialize values
  dcopy_(&system_solution_length_, inner_solution_1_temp, &increment, inner_solution_1_, &increment);
  dcopy_(&system_solution_length_, inner_solution_2_temp, &increment, inner_solution_2_, &increment);
  dcopy_(&system_solution_length_, inner_solution_3_temp, &increment, inner_solution_3_, &increment);
  dcopy_(&system_solution_length_, inner_solution_ls_temp, &increment, inner_solution_ls_, &increment);
  dcopy_(&system_solution_length_, inner_solution_trial_temp, &increment, inner_solution_trial_, &increment);
  dcopy_(&system_solution_length_, new_system_vector_temp, &increment, new_system_vector_, &increment);
  dcopy_(&system_solution_length_, right_hand_side_temp, &increment, right_hand_side_, &increment);
  dcopy_(&system_solution_length_, system_solution_temp, &increment, system_solution_, &increment);
  dcopy_(&system_solution_length_, system_solution_best_temp, &increment, system_solution_best_, &increment);

  // Delete arrays (in case they exist)
  if (inner_solution_1_temp != nullptr) {
    delete[] inner_solution_1_temp;
    inner_solution_1_temp = nullptr;
  } // end if
  if (inner_solution_2_temp != nullptr) {
    delete[] inner_solution_2_temp;
    inner_solution_2_temp = nullptr;
  } // end if
  if (inner_solution_3_temp != nullptr) {
    delete[] inner_solution_3_temp;
    inner_solution_3_temp = nullptr;
  } // end if
  if (inner_solution_ls_temp != nullptr) {
    delete[] inner_solution_ls_temp;
    inner_solution_ls_temp = nullptr;
  } // end if
  if (inner_solution_trial_temp != nullptr) {
    delete[] inner_solution_trial_temp;
    inner_solution_trial_temp = nullptr;
  } // end if
  if (new_system_vector_temp != nullptr) {
    delete[] new_system_vector_temp;
    new_system_vector_temp = nullptr;
  } // end if
  if (right_hand_side_temp != nullptr) {
    delete[] right_hand_side_temp;
    right_hand_side_temp = nullptr;
  } // end if
  if (system_solution_temp != nullptr) {
    delete[] system_solution_temp;
    system_solution_temp = nullptr;
  } // end if
  if (system_solution_best_temp != nullptr) {
    delete[] system_solution_best_temp;
    system_solution_best_temp = nullptr;
  } // end if

  // Update lengths
  system_solution_length_ = system_solution_length_new;

} // end resizeSystemSolution

// Set augment
bool QPSolverInteriorPoint::setAugment(const Reporter* reporter,
                                       int set,
                                       int index,
                                       double system_vector[],
                                       double solution1[],
                                       double solution2[],
                                       double augmentation_value)
{

  // Initialize return value
  bool success_without_factorization_error = true;

  // Check set to augment
  if (set == 1) {

    // Print message
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "Augmenting set 1\n");

    // Add solution element
    for (int i = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size(); i > (int)omega_positive_.size(); i--) {
      system_solution_[i] = system_solution_[i - 1];
    }
    system_solution_[(int)omega_positive_.size()] = augmentation_value;

    // Add element to set
    omega_positive_.push_back(index);

    // Evaluate new system vector
    evaluateSystemVector(set, index, system_vector);

    // Augment Cholesky factor
    success_without_factorization_error = choleskyAugment(system_vector, (int)omega_positive_.size() - 1, solution1, 1, solution2, vector_[index]);

  } // end if

  else if (set == 2) {

    // Print message
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "Augmenting set 2\n");

    // Add solution element
    for (int i = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size(); i > (int)omega_positive_.size() + (int)gamma_positive_.size(); i--) {
      system_solution_[i] = system_solution_[i - 1];
    }
    system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size()] = augmentation_value;

    // Add element to set
    gamma_positive_.push_back(index);

    // Evaluate new system vector
    evaluateSystemVector(set, index, system_vector);

    // Augment Cholesky factor
    success_without_factorization_error = choleskyAugment(system_vector, (int)omega_positive_.size() + (int)gamma_positive_.size() - 1, solution1, 0, solution2, -scalar_);

  } // end else if

  else {

    // Print message
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "Augmenting set 3\n");

    // Add solution element
    system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size()] = augmentation_value;

    // Add element to set
    gamma_negative_.push_back(index);

    // Evaluate new system vector
    evaluateSystemVector(set, index, system_vector);

    // Augment Cholesky factor
    success_without_factorization_error = choleskyAugment(system_vector, (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() - 1, solution1, 0, solution2, scalar_);

  } // end else

  // Return
  return success_without_factorization_error;

} // end setAugment

// Set delete
void QPSolverInteriorPoint::setDelete(const Reporter* reporter,
                                      int set,
                                      int index,
                                      double solution1[],
                                      double solution2[])
{

  // Check set from which to delete
  if (set == 1) {

    // Remove solution element
    for (int i = index; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() - 1; i++) {
      system_solution_[i] = system_solution_[i + 1];
    }

    // Remove element from set
    omega_positive_.erase(omega_positive_.begin() + index);

    // Delete row/column from Cholesky factor
    choleskyDelete(index, solution1, solution2);

  } // end if

  else if (set == 2) {

    // Remove solution element
    for (int i = (int)omega_positive_.size() + index; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() - 1; i++) {
      system_solution_[i] = system_solution_[i + 1];
    }

    // Remove element from set
    gamma_positive_.erase(gamma_positive_.begin() + index);

    // Delete row/column from Cholesky factor
    choleskyDelete((int)omega_positive_.size() + index, solution1, solution2);

  } // end else if

  else {

    // Remove solution element
    for (int i = (int)omega_positive_.size() + (int)gamma_positive_.size() + index; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() - 1; i++) {
      system_solution_[i] = system_solution_[i + 1];
    }

    // Remove element from set
    gamma_negative_.erase(gamma_negative_.begin() + index);

    // Delete row/column from Cholesky factor
    choleskyDelete((int)omega_positive_.size() + (int)gamma_positive_.size() + index, solution1, solution2);

  } // end else

} // end setDelete

// Solve linear system
void QPSolverInteriorPoint::solveSystem(double right_hand_side[],
                                        double solution[])
{

  // Check if resizing needed
  if ((int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() > system_solution_length_) {
    resizeSystemSolution();
  }

  // Set inputs for BLASLAPACK
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  char upper_lower = 'L'; // Use "lower" since Fortran uses column-major ordering
  char transpose = 'T';   // Use "transpose" since Fortran uses column-major ordering
  char diagonal = 'N';
  int increment1 = 1;
  int incrementn = system_solution_length_;

  // Copy right_hand_side to solution
  dcopy_(&length, right_hand_side, &increment1, solution, &increment1);

  // Solve system
  dtrsv_(&upper_lower, &transpose, &diagonal, &length, factor_, &incrementn, solution, &increment1);

} // end solveSystem

// Solve triangular system with transpose
void QPSolverInteriorPoint::solveSystemTranspose(double right_hand_side[],
                                                 double solution[])
{

  // Check if resizing needed
  if ((int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() > system_solution_length_) {
    resizeSystemSolution();
  }

  // Set inputs for BLASLAPACK
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  char upper_lower = 'L'; // Use "lower" since Fortran uses column-major ordering
  char transpose = 'N';   // Use "not transpose" since Fortran uses column-major ordering
  char diagonal = 'N';
  int increment1 = 1;
  int incrementn = system_solution_length_;

  // Copy right_hand_side to solution
  dcopy_(&length, right_hand_side, &increment1, solution, &increment1);

  // Solve system
  dtrsv_(&upper_lower, &transpose, &diagonal, &length, factor_, &incrementn, solution, &increment1);

} // end solveSystemTranspose

// // Solve linear system
// void QPSolverDualActiveSet::solveLinearSystem(int size,
//                                               double matrix[],
//                                               double right_hand_side[],
//                                               double solution[])
// {

//   // Set inputs for BLASLAPACK
//   char upper_lower = 'L';
//   int length = size;
//   int length_squared = size*size;
//   int nrhs = 1;
//   int increment = 1;
//   int info;
//   int piv_info [size];
//   int lw;
//   double matrix_copy [length_squared];
//   double work;

//   // Copy right_hand_side to solution
//   dcopy_(&length_squared, matrix, &increment, matrix_copy, &increment);
//   dcopy_(&length, right_hand_side, &increment, solution, &increment);

//   // Solve system
//   dsysv_(&upper_lower, &length, &nrhs, matrix_copy, &length, piv_info, solution, &length, &work, &lw, &info);

// } // end solveLinearSystem


void QPSolverInteriorPoint::solveLinearSystem(int size,
                                              double matrix[],
                                              double right_hand_side[],
                                              double solution[]) {
    char uplo = 'L';  // Lower triangular storage
    int n = size;
    int nrhs = 1;
    int info;
    int* ipiv = new int[size];
    
    // Optimal workspace query
    double work_query;
    int lwork = -1;
    dsysv_(&uplo, &n, &nrhs, matrix, &n, ipiv, solution, &n, &work_query, &lwork, &info);
    
    // Allocate optimal workspace
    lwork = static_cast<int>(work_query);
    double* work = new double[lwork];
    
    // Copy right-hand side
    std::copy(right_hand_side, right_hand_side + size, solution);
    
    // Solve the system
    dsysv_(&uplo, &n, &nrhs, matrix, &n, ipiv, solution, &n, work, &lwork, &info);
    
    // Clean up
    delete[] ipiv;
    delete[] work;
}

void QPSolverInteriorPoint::calculate_step_sizes(const Vector& delta_x, 
                            const double& delta_y, 
                            const Vector& delta_z, 
                            const double& r_pri_, 
                            const Vector& r_dual_, 
                            const Vector& vector_x_, 
                            const Vector& vector_z_,
                            const Vector& vector_c_, 
                            SymmetricMatrixDense& matrix_Q_, 
                            const Vector& matrix_At_, 
                            double scalar_y_, 
                            double scalar_b_, 
                            double beta_, 
                            double& alpha_x, 
                            double& alpha_y, 
                            double& alpha_z) {
    int total_length = vector_x_.length();

    // Step 1: Prepare vectors
    Vector vector_r(total_length + 1, 0.0);
    vector_r.set(0, r_pri_);
    for (int i = 1; i < total_length + 1; i++) {
        vector_r.set(i, r_dual_.values()[i - 1]);
    }

    Vector vector_sf(total_length, 0.0);
    SymmetricMatrixDense matrix_for_s;
    matrix_for_s.setAsDiagonal(total_length, 0);
    for (int i = 0; i < total_length; i++) {
        for (int j = i; j < total_length; j++) {
            matrix_for_s.setElement(i, j, matrix_Q_.element(i, j));
        }
    }

    matrix_for_s.matrixVectorProduct(delta_x, vector_sf);

    Vector vector_s(total_length + 1, 0.0);
    vector_s.set(0, -matrix_At_.innerProduct(delta_x));
    for (int i = 1; i < total_length + 1; i++) {
        vector_s.set(i, vector_sf.values()[i - 1]);
    }

    Vector vector_t(total_length + 1, 0.0);
    for (int i = 1; i < total_length + 1; i++) {
        vector_t.set(i, -matrix_At_.values()[i-1]);
    }
    vector_t.scale(delta_y);

    Vector vector_u(total_length + 1, 0.0);
    for (int i = 1; i < total_length + 1; i++) {
        vector_u.set(i, -delta_z.values()[i - 1]);
    }

    // Step 2: Compute bar_alpha, bar_alpha_x, bar_alpha_z
    double bar_alpha_x, bar_alpha_z, bar_alpha;

    bar_alpha_x = 1;
    bar_alpha_z = 1;
    for (int i = 0; i < total_length; i++){
      if (delta_x.values()[i] < 0){
        bar_alpha_x = fmin( bar_alpha_x, (beta_ -1)*vector_x_.values()[i]/delta_x.values()[i]);
      }
      if (delta_z.values()[i] < 0){
        bar_alpha_z = fmin( bar_alpha_z, (beta_ - 1)*vector_z_.values()[i]/delta_z.values()[i]);
      }
    }
    bar_alpha = std::fmin(bar_alpha_x, bar_alpha_z);


    // Step 3: Compute QP coefficients...
    Vector s_plus_u(total_length + 1, 0.0);
    s_plus_u.copy(vector_s);
    s_plus_u.addScaledVector(1, vector_u);
    double Q1_11 = s_plus_u.innerProduct(s_plus_u) + delta_x.innerProduct(delta_z);
    double Q1_12 = s_plus_u.innerProduct(vector_t);
    double Q1_22 = vector_t.innerProduct(vector_t);
    double c1_1 = vector_r.innerProduct(s_plus_u) + 0.5 * (delta_x.innerProduct(vector_z_) + vector_x_.innerProduct(delta_z));
    double c1_2 = vector_r.innerProduct(vector_t);

    double alpha_x1 = (c1_2 * Q1_12 / Q1_22 - c1_1) / (Q1_11 - Q1_12 * Q1_12 / Q1_22);
    alpha_x1 = std::max(0.0, std::min(alpha_x1, bar_alpha));  // clamp alpha_x1

    double alpha_y1 = (-c1_2 - alpha_x1 * Q1_12) / Q1_22;
    double alpha_z1 = alpha_x1;

    double alpha_x2, alpha_y2, alpha_z2;
    if (bar_alpha_x <= bar_alpha_z) {
        double Q2_11 = vector_t.innerProduct(vector_t);
        double Q2_12 = vector_t.innerProduct(vector_u);
        double Q2_22 = vector_u.innerProduct(vector_u);
        double c2_1 = vector_r.innerProduct(vector_t) + bar_alpha_x * vector_s.innerProduct(vector_t);
        double c2_2 = vector_r.innerProduct(vector_u) + 0.5 * vector_x_.innerProduct(delta_z)
                      + bar_alpha_x * (vector_s.innerProduct(vector_u) + 0.5 * delta_x.innerProduct(delta_z));

        alpha_x2 = bar_alpha_x;
        alpha_z2 = (c2_1 * Q2_12 / Q2_11 - c2_2) / (Q2_22 - Q2_12 * Q2_12 / Q2_11);
        alpha_z2 = std::max(bar_alpha, std::min(alpha_z2, bar_alpha_z));  // clamp alpha_z2

        alpha_y2 = (-c2_1 - alpha_z2 * Q2_12) / Q2_11;
    } else {
        double Q2_11 = vector_s.innerProduct(vector_s);
        double Q2_12 = vector_s.innerProduct(vector_t);
        double Q2_22 = vector_t.innerProduct(vector_t);
        double c2_1 = vector_r.innerProduct(vector_s) + 0.5 * delta_x.innerProduct(vector_z_)
                      + bar_alpha_z * (vector_s.innerProduct(vector_u) + 0.5 * delta_x.innerProduct(delta_z));
        double c2_2 = vector_r.innerProduct(vector_t) + bar_alpha_z * vector_t.innerProduct(vector_u);

        alpha_x2 = (c2_2 * Q2_12 / Q2_22 - c2_1) / (Q2_11 - Q2_12 * Q2_12 / Q2_22);
        alpha_x2 = std::max(bar_alpha, std::min(alpha_x2, bar_alpha_x));  // clamp alpha_x2

        alpha_y2 = (-c2_2 - alpha_x2 * Q2_12) / Q2_22;
        alpha_z2 = bar_alpha_z;
    }

    // Step 4: Compute the final step sizes
    Vector vector_x1(total_length, 0.0), vector_z1(total_length, 0.0), vector_x2(total_length, 0.0), vector_z2(total_length, 0.0);
    double scalar_y1, scalar_y2;

    vector_x1.copy(vector_x_);
    vector_z1.copy(vector_z_);
    vector_x1.addScaledVector(alpha_x1, delta_x);
    vector_z1.addScaledVector(alpha_z1, delta_z);
    scalar_y1 = scalar_y_ + alpha_y1 * delta_y;

    vector_x2.copy(vector_x_);
    vector_z2.copy(vector_z_);
    vector_x2.addScaledVector(alpha_x2, delta_x);
    vector_z2.addScaledVector(alpha_z2, delta_z);
    scalar_y2 = scalar_y_ + alpha_y2 * delta_y;

    // Step 5: Choose the best step sizes
    Vector QX1(total_length, 0.0);
    matrix_Q_.matrixVectorProduct(vector_x1, QX1);

    Vector oper11(total_length, 0.0);
    oper11.copy(QX1);
    oper11.addScaledVector(-scalar_y1, matrix_At_);

    Vector secondnorm1(total_length, 0.0), secondnorm2(total_length, 0.0);
    secondnorm1.copy(vector_c_);
    secondnorm1.addScaledVector(-1, vector_z1);
    secondnorm1.addScaledVector(1, oper11);

    Vector QX2(total_length, 0.0);
    matrix_Q_.matrixVectorProduct(vector_x2, QX2);

    Vector oper12(total_length, 0.0);
    oper12.copy(QX2);
    oper12.addScaledVector(-scalar_y2, matrix_At_);

    secondnorm2.copy(vector_c_);
    secondnorm2.addScaledVector(-1, vector_z2);
    secondnorm2.addScaledVector(1, oper12);

    double firstnorm1 = matrix_At_.innerProduct(vector_x1) - scalar_b_;
    double firstnorm2 = matrix_At_.innerProduct(vector_x2) - scalar_b_;

    if ((firstnorm1)*(firstnorm1) + (secondnorm1.norm2())*(secondnorm1.norm2()) + vector_x1.innerProduct(vector_z1)
        < (firstnorm2)*(firstnorm2) + (secondnorm2.norm2())*(secondnorm2.norm2()) + vector_x2.innerProduct(vector_z2)) {
        alpha_x = alpha_x1;
        alpha_y = alpha_y1;
        alpha_z = alpha_z1;
    } else {
        alpha_x = alpha_x2;
        alpha_y = alpha_y2;
        alpha_z = alpha_z2;
    }
}


} // namespace NonOpt