// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Baoyu Zhou

#include <algorithm>
#include <cmath>
#include <iterator>

#include "NonOptBLAS.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptQPSolverActiveSet.hpp"

namespace NonOpt
{

// Constructor
QPSolverActiveSet::QPSolverActiveSet()
    : factor_(nullptr),
      inner_solution_1_(nullptr),
      inner_solution_2_(nullptr),
      system_solution_(nullptr),
      system_solution_best_(nullptr) {}

// Destructor
QPSolverActiveSet::~QPSolverActiveSet()
{

  // Delete arrays
  if (inner_solution_1_ != nullptr) {
    delete[] inner_solution_1_;
  }
  if (inner_solution_2_ != nullptr) {
    delete[] inner_solution_2_;
  }
  if (system_solution_ != nullptr) {
    delete[] system_solution_;
  }
  if (system_solution_best_ != nullptr) {
    delete[] system_solution_best_;
  }
  if (factor_ != nullptr) {
    delete[] factor_;
  }

}  // end destructor

// Add options
void QPSolverActiveSet::addOptions(Options* options,
                                   const Reporter* reporter)
{

  // Add bool options
  options->addBoolOption(reporter,
                         "QPAS_fail_on_factorization_error",
                         false,
                         "Indicator for whether to indicate failure on factorization error.\n"
                         "Default value: false.");
  options->addBoolOption(reporter,
                         "QPAS_allow_inexact_termination",
                         true,
                         "Indicator for whether to allow early termination.\n"
                         "Default value: true.");
  // Add double options
  options->addDoubleOption(reporter,
                           "QPAS_kkt_tolerance",
                           1e-12,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for determining optimality.  If the KKT error\n"
                           "computed by the algorithm falls below this tolerance, then\n"
                           "the algorithm terminates with a message of success.\n"
                           "Default value: 1e-12.");
  options->addDoubleOption(reporter,
                           "QPAS_cholesky_tolerance",
                           1e-12,
                           0.0,
                           1.0,
                           "Tolerance for small (or negative values) when updating the\n"
                           "Cholesky factorization of an augmented matrix.  If a diagonal\n"
                           "value in the factorization falls below this tolerance, then\n"
                           "the factorization may be re-computed from scratch and/or the\n"
                           "diagonal value may be replaced by this tolerance value.\n"
                           "Default value: 1e-12.");
  options->addDoubleOption(reporter,
                           "QPAS_inexact_termination_factor",
                           1.5,
                           0.0,
                           4.0,
                           "Factor for inexact termination, if allowed.  Factor by which\n"
                           "norm of inexact solution needs to be within true norm of true\n"
                           "(unknown) projection of origin onto convex hull of gradients.\n"
                           "Default value: 1.5.");
  options->addDoubleOption(reporter,
                           "QPAS_inexact_termination_ratio_min",
                           1e-02,
                           0.0,
                           1.0,
                           "Minimum value for ratio used in inexact termination condition.\n"
                           "Default value: 1e-02.");
  options->addDoubleOption(reporter,
                           "QPAS_linear_independence_tolerance",
                           1e-12,
                           0.0,
                           1.0,
                           "Tolerance for a linear independence check when adding a new\n"
                           "vector to an augmented matrix.  If a residual value falls\n"
                           "below this tolerance, then indices are exchanged to maintain\n"
                           "linear independence of the augmented matrix.\n"
                           "Default value: 1e-12.");

  // Add integer options
  options->addIntegerOption(reporter,
                            "QPAS_iteration_limit_minimum",
                            1e+02,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Minimum limit on the number of iterations.\n"
                            "Default value: 1e+02.");
  options->addIntegerOption(reporter,
                            "QPAS_iteration_limit_maximum",
                            1e+04,
                            0.0,
                            NONOPT_INT_INFINITY,
                            "Maximum limit on the number of iterations.\n"
                            "Default value: 1e+04.");

}  // end addOptions

// Set options
void QPSolverActiveSet::setOptions(const Options* options,
                                   const Reporter* reporter)
{

  // Read bool options
  options->valueAsBool(reporter, "QPAS_fail_on_factorization_error", fail_on_factorization_error_);
  options->valueAsBool(reporter, "QPAS_allow_inexact_termination", allow_inexact_termination_);

  // Read double options
  options->valueAsDouble(reporter, "QPAS_kkt_tolerance", kkt_tolerance_);
  options->valueAsDouble(reporter, "QPAS_cholesky_tolerance", cholesky_tolerance_);
  options->valueAsDouble(reporter, "QPAS_inexact_termination_factor", inexact_termination_factor_);
  options->valueAsDouble(reporter, "QPAS_inexact_termination_ratio_min", inexact_termination_ratio_min_);
  options->valueAsDouble(reporter, "QPAS_linear_independence_tolerance", linear_independence_tolerance_);

  // Read integer options
  options->valueAsInteger(reporter, "QPAS_iteration_limit_minimum", iteration_limit_minimum_);
  options->valueAsInteger(reporter, "QPAS_iteration_limit_maximum", iteration_limit_maximum_);

}  // end setOptions

// Initialize
void QPSolverActiveSet::initialize(const Options* options,
                                   Quantities* quantities,
                                   const Reporter* reporter)
{
  initializeData(quantities->numberOfVariables());
}

// Initialize data
void QPSolverActiveSet::initializeData(int gamma_length)
{

  // Set length parameters
  gamma_length_ = gamma_length;
  system_solution_length_ = 2 * gamma_length_;
  factor_length_ = system_solution_length_ * system_solution_length_;

  // Initialize problem data
  matrix_ = nullptr;
  vector_list_.clear();
  vector_.clear();
  scalar_ = -1.0;

  // Delete arrays (in case they exist)
  if (inner_solution_1_ != nullptr) {
    delete[] inner_solution_1_;
  }
  if (inner_solution_2_ != nullptr) {
    delete[] inner_solution_2_;
  }
  if (system_solution_ != nullptr) {
    delete[] system_solution_;
  }
  if (system_solution_best_ != nullptr) {
    delete[] system_solution_best_;
  }
  if (factor_ != nullptr) {
    delete[] factor_;
  }

  // Allocate arrays
  inner_solution_1_ = new double[system_solution_length_];
  inner_solution_2_ = new double[system_solution_length_];
  system_solution_ = new double[system_solution_length_];
  system_solution_best_ = new double[system_solution_length_];
  factor_ = new double[factor_length_];

  // Set inputs for blas
  double value = 0.0;
  int increment1 = 0;
  int increment2 = 1;

  // Initialize values
  dcopy_(&system_solution_length_, &value, &increment1, inner_solution_1_, &increment2);
  dcopy_(&system_solution_length_, &value, &increment1, inner_solution_2_, &increment2);
  dcopy_(&system_solution_length_, &value, &increment1, system_solution_, &increment2);
  dcopy_(&system_solution_length_, &value, &increment1, system_solution_best_, &increment2);
  dcopy_(&factor_length_, &value, &increment1, factor_, &increment2);

  // Set status
  setStatus(QP_UNSET);

  // Set null solution
  setNullSolution();

}  // end initializeData

// Get combination norm
double QPSolverActiveSet::combinationNormInf()
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, combination_.values(), &increment);

  // Set norm
  double combination_norm_inf = fabs(combination_.values()[i - 1]);

  // Set maximum
  if (isnan(combination_norm_inf) || combination_norm_inf > NONOPT_DOUBLE_INFINITY) {
    combination_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return combination_norm_inf;

}  // end combinationNormInf

// Get translated combination norm
double QPSolverActiveSet::combinationTranslatedNormInf()
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, combination_translated_.values(), &increment);

  // Set norm
  double combination_translated_norm_inf = fabs(combination_translated_.values()[i - 1]);

  // Set maximum
  if (isnan(combination_translated_norm_inf) || combination_translated_norm_inf > NONOPT_DOUBLE_INFINITY) {
    combination_translated_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return combination_translated_norm_inf;

}  // end combinationTranslatedNormInf

// Get dual step
void QPSolverActiveSet::dualStep(double vector[])
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, dual_step_.values(), &increment, vector, &increment);

}  // end dualStep

// Get feasible dual step
void QPSolverActiveSet::dualStepFeasible(double vector[])
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, dual_step_feasible_.values(), &increment, vector, &increment);

}  // end dualStepFeasible

// Get dual step norm
double QPSolverActiveSet::dualStepNormInf()
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, dual_step_.values(), &increment);

  // Set norm
  double dual_step_norm_inf = fabs(dual_step_.values()[i - 1]);

  // Set maximum
  if (isnan(dual_step_norm_inf) || dual_step_norm_inf > NONOPT_DOUBLE_INFINITY) {
    dual_step_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return dual_step_norm_inf;

}  // end dualStepNormInf

// Get feasible dual step norm
double QPSolverActiveSet::dualStepFeasibleNormInf()
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Find index of element with maximum absolute value
  // (returns index from 1,...,length)
  int i = idamax_(&length, dual_step_feasible_.values(), &increment);

  // Set norm
  double dual_step_norm_inf = fabs(dual_step_feasible_.values()[i - 1]);

  // Set maximum
  if (isnan(dual_step_norm_inf) || dual_step_norm_inf > NONOPT_DOUBLE_INFINITY) {
    dual_step_norm_inf = NONOPT_DOUBLE_INFINITY;
  }

  // Return inf-norm
  return dual_step_norm_inf;

}  // end dualStepFeasibleNormInf

// Get objective quadratic value
double QPSolverActiveSet::objectiveQuadraticValue()
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Set quadratic value
  double objective_quadratic_value = -ddot_(&length, dual_step_.values(), &increment, combination_translated_.values(), &increment);

  // Set maximum
  if (isnan(objective_quadratic_value) || objective_quadratic_value > NONOPT_DOUBLE_INFINITY) {
    objective_quadratic_value = NONOPT_DOUBLE_INFINITY;
  }

  // Return
  return objective_quadratic_value;

}  // end objectiveQuadraticValue

// Get objective quadratic value for feasible dual step
double QPSolverActiveSet::objectiveQuadraticValueFeasible()
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Set quadratic value
  double objective_quadratic_value = -ddot_(&length, dual_step_.values(), &increment, combination_translated_.values(), &increment);

  // Scale so value corresponds to feasible dual step
  objective_quadratic_value *= pow(dual_step_projection_scalar_, 2.0);

  // Set maximum
  if (isnan(objective_quadratic_value) || objective_quadratic_value > NONOPT_DOUBLE_INFINITY) {
    objective_quadratic_value = NONOPT_DOUBLE_INFINITY;
  }

  // Return
  return objective_quadratic_value;

}  // end objectiveQuadraticValueFeasible

// Get solution, gamma
void QPSolverActiveSet::gamma(double vector[])
{

  // Set inputs for blas
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, gamma_.values(), &increment, vector, &increment);

}  // end gamma

// Get full KKT error
double QPSolverActiveSet::KKTErrorFull()
{

  // Evaluate gradient combination
  Vector Gomega(gamma_length_);
  for (int i = 0; i < (int)vector_list_.size(); i++) {
    Gomega.addScaledVector(omega_.values()[i], *vector_list_[i].get());
  }

  // Declare dual vector
  Vector d0(gamma_length_);
  d0.addScaledVector(-1.0, Gomega);
  d0.addScaledVector(-1.0, gamma_);
  Vector d(gamma_length_);
  matrix_->matrixVectorProductOfInverse(d0, d);

  // Evaluate gradient inner products
  Vector Gtd;
  Gtd.setLength((int)vector_.size());
  for (int i = 0; i < (int)vector_list_.size(); i++) {
    Gtd.values()[i] = vector_list_[i]->innerProduct(d);
  }

  // Evaluate dual scalar
  double z = omega_.innerProduct(Gtd);
  for (int i = 0; i < (int)vector_list_.size(); i++) {
    z = z + omega_.values()[i] * vector_[i];
  }

  // Initialize KKT error
  double kkt_error = 0.0;

  // Loop through KKT conditions, update error
  for (int j = 0; j < (int)vector_.size(); j++) {
    double term = -z + vector_[j] + Gtd.values()[j];
    if (term > 0.0 && term > kkt_error) {
      kkt_error = term;
    }
  }  // end for
  for (int i = 0; i < gamma_length_; i++) {
    double term = -scalar_ + d.values()[i];
    if (term > 0.0 && term > kkt_error) {
      kkt_error = term;
    }
    term = -scalar_ - d.values()[i];
    if (term > 0.0 && term > kkt_error) {
      kkt_error = term;
    }
  }  // end for
  for (int j = 0; j < (int)vector_.size(); j++) {
    double term = -omega_.values()[j];
    if (term > 0.0 && term > kkt_error) {
      kkt_error = term;
    }
  }  // end for
  double sum = 0.0;
  for (int j = 0; j < (int)vector_.size(); j++) {
    sum = sum + omega_.values()[j];
  }
  if (fabs(sum - 1.0) > kkt_error) {
    kkt_error = fabs(sum - 1.0);
  }
  for (int i = 0; i < gamma_length_; i++) {
    double term1 = gamma_.values()[i];
    double term2 = scalar_ - d.values()[i];
    if (term1 > 0.0 && term1 * term2 > kkt_error) {
      kkt_error = term1 * term2;
    }
    term1 = -gamma_.values()[i];
    term2 = scalar_ + d.values()[i];
    if (term1 > 0.0 && term1 * term2 > kkt_error) {
      kkt_error = term1 * term2;
    }
  }  // end for

  // Set maximum
  if (isnan(kkt_error) || kkt_error > NONOPT_DOUBLE_INFINITY) {
    kkt_error = NONOPT_DOUBLE_INFINITY;
  }

  // Return
  return kkt_error;

}  // end KKTErrorFull

// Get solution, omega
void QPSolverActiveSet::omega(double vector[])
{

  // Set inputs for blas
  int length = (int)vector_.size();
  int increment = 1;

  // Copy values
  dcopy_(&length, omega_.values(), &increment, vector, &increment);

}  // end omega

// Initialize data
void QPSolverActiveSet::setNullSolution()
{

  // Algorithm parameters
  iteration_count_ = 0;
  kkt_error_ = -NONOPT_DOUBLE_INFINITY;
  dual_objective_best_ = NONOPT_DOUBLE_INFINITY;
  primal_objective_reference_ = -NONOPT_DOUBLE_INFINITY;

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
  combination_translated_.setLength(gamma_length_);
  dual_step_.setLength(gamma_length_);
  dual_step_projection_scalar_ = 0.0;
  dual_step_feasible_.setLength(gamma_length_);
  dual_step_feasible_best_.setLength(gamma_length_);

}  // end setNullSolution

// Add vectors
void QPSolverActiveSet::addData(const std::vector<std::shared_ptr<Vector> > vector_list,
                                const std::vector<double> vector)
{

  // Loop through new elements
  for (int i = 0; i < (int)vector.size(); i++) {
    vector_list_.push_back(vector_list[i]);
    vector_.push_back(vector[i]);
  }

}  // end addData

// Inexact termination condition
bool QPSolverActiveSet::inexactTerminationCondition(const Reporter* reporter)
{

  // Finalize solution (to set omega_ and gamma_)
  //  if (iteration_count_ == 0){
  //    finalizeSolution();
  //  };

  // Evaluate dual vectors
  //evaluateDualVectors();

  // Set "b"
  Vector b((int)vector_.size());
  for (int i = 0; i < (int)vector_.size(); i++) {
    b.values()[i] = vector_[i];
  }

  // Evaluate primal objective value
  double primal_objective = -0.5 * objectiveQuadraticValue();  // + omega_.innerProduct(b) - scalar_*gamma_.norm1();
  //printf("second term is= %.4e and third term =%.4e\n",omega_.innerProduct(b),scalar_*gamma_.norm1());

  // Evaluate dual objective value
  Vector bPlusGd((int)vector_.size());
  bPlusGd.copy(b);
  for (int i = 0; i < (int)vector_.size(); i++) {
    bPlusGd.values()[i] += vector_list_[i]->innerProduct(dual_step_feasible_);
  }
  double dual_objective = bPlusGd.max() + 0.5 * objectiveQuadraticValueFeasible();

  // Initialize in first iteration, else update best dual objective value
  if (iteration_count_ == 0) {
    //printf("\n");
    //printf("inexact_termination_factor_ = %e\n",inexact_termination_factor_);
    //printf("gamma^2 + 2*gamma           = %e\n",pow(inexact_termination_factor_,2.0) + 2*inexact_termination_factor_);
    //printf("%11s  %11s  %11s  %11s  %11s  %15s\n","prim-ref","prim-obj","dual-obj","dual-best","obj-ratio","term-ratio");
    dual_objective_best_ = dual_objective;
    dual_step_feasible_best_.copy(dual_step_feasible_);
    primal_objective_reference_ = primal_objective;
  }
  else {
    if (dual_objective < dual_objective_best_) {
      dual_objective_best_ = dual_objective;
      dual_step_feasible_best_.copy(dual_step_feasible_);
    }
  }

  // Set objective ratio
  double objective_ratio = (b.max() - primal_objective_reference_) / (b.max() - dual_objective_best_);

  // Set inexact termination ratio value
  double inexact_termination_ratio = 1 - (pow(inexact_termination_factor_, 2.0) + 2 * inexact_termination_factor_) / (objective_ratio - 1.0);

  if (iteration_count_ == 0) {
    reporter->printf(R_QP, R_PER_INNER_ITERATION_IN, "==========================================================================================================================================================\n");
    reporter->printf(R_QP, R_PER_INNER_ITERATION_IN, "b              b-p0          b-pk        b-qk    b-qkb       xi          alpha             LHS             RHS           LHS1        RHS1         kkt\n");
    reporter->printf(R_QP, R_PER_INNER_ITERATION_IN, "==========================================================================================================================================================\n");
  }
  reporter->printf(R_QP, R_PER_INNER_ITERATION_IN, "%+.4e  %+.4e %+.4e %+.4e  %+.4e  %+.4e  %+.8e  %+.10e  %+.10e  %+.4e  %+.4e  %+.4e\n",
                   b.max(),
                   b.max() - primal_objective_reference_,
                   b.max() - primal_objective,
                   b.max() - dual_objective,
                   b.max() - dual_objective_best_,
                   objective_ratio,
                   inexact_termination_ratio,
                   primal_objective - primal_objective_reference_,
                   fmax(inexact_termination_ratio, inexact_termination_ratio_min_) * (dual_objective_best_ - primal_objective_reference_),
                   (pow(inexact_termination_factor_, 2.0) + 2 * inexact_termination_factor_) * (b.max() - dual_objective_best_),
                   dual_objective_best_ - primal_objective,
                   kkt_error_);

  // Initialize return value
  bool condition_bool = false;

  Vector d_bar(gamma_length_);
  for (int i = 0; i < (int)vector_list_.size(); i++) {
    d_bar.addScaledVector(1.0 / (int)vector_list_.size(), *vector_list_[i].get());
  }

  Vector gTd_bar((int)vector_.size());
  for (int i = 0; i < (int)vector_.size(); i++) {
    gTd_bar.values()[i] = -1.0 * (vector_list_[i]->innerProduct(d_bar));
  }

  int length = (int)vector_.size();
  int increment = 1;
  int i_ind = idamax_(&length, gTd_bar.values(), &increment);
  int gTd_bar_opt = -1.0 * gTd_bar.values()[i_ind];
  //    int alpha_step=alpha_step1;

  double alpha_step = 0;
  Vector dual_step_simple_(gamma_length_);
  if (gTd_bar_opt > 0) {

    //      int length1 = gamma_length_;

    //      Vector combination_bar_(length1);
    //      combination_bar_.copyArray(combination_.values());
    // Set quadratic value dHd
    double quadratic_value_bar = matrix_->innerProduct(d_bar);

    // Set maximum
    if (isnan(quadratic_value_bar) || quadratic_value_bar > NONOPT_DOUBLE_INFINITY) {
      quadratic_value_bar = NONOPT_DOUBLE_INFINITY;
    }

    alpha_step = gTd_bar_opt / quadratic_value_bar;

    dual_step_simple_.addScaledVector(-1.0 * alpha_step, d_bar);

    dual_objective_simple_ = alpha_step * gTd_bar.values()[i_ind] + 0.5 * alpha_step * alpha_step * quadratic_value_bar;

    //      Vector bPlusGd_simple((int)vector_.size());
    //      bPlusGd_simple.copy(b);
    //
    //      for (int i = 0; i < (int)vector_.size(); i++) {
    //        bPlusGd_simple.values()[i] += vector_list_[i]->innerProduct(dual_step_simple_);
    //      }
    //      dual_objective_simple_ = bPlusGd_simple.max() + 0.5*quadratic_value_bar;
    //printf("condition met\n");
  }

  // dual_objective_best_ < b.max() &&     alpha_step1>0&&
  //dual_objective_best_ - primal_objective <= (pow(inexact_termination_factor_,2.0) + 2*inexact_termination_factor_)*(b.max() - dual_objective_best_)
  if (alpha_step > 0 &&
      (primal_objective - primal_objective_reference_ >=
       fmax(inexact_termination_ratio, inexact_termination_ratio_min_) * (dual_objective_simple_ - primal_objective_reference_))) {
    //printf("dual step length is %4d\n",(int)dual_step_.length());
    //printf("dual step simple length is %4d\n",(int)dual_step_simple_.length());
    dual_step_.copy(dual_step_simple_);
    condition_bool = true;
  }

  //dual_objective_best_ < b.max() &&
  // Check condition
  else if (
      (primal_objective - primal_objective_reference_ >=
           fmax(inexact_termination_ratio, inexact_termination_ratio_min_) * (dual_objective_best_ - primal_objective_reference_) ||
       dual_objective_best_ - primal_objective <= (pow(inexact_termination_factor_, 2.0) + 2 * inexact_termination_factor_) * (b.max() - dual_objective_best_))) {
    dual_step_feasible_.copy(dual_step_feasible_best_);
    condition_bool = true;
  }

  // Return
  return condition_bool;
}

// Solve
void QPSolverActiveSet::solveQP(const Options* options,
                                const Reporter* reporter)
{

  // Initialize values
  setStatus(QP_UNSET);
  setNullSolution();

  // try to solve QP, terminate on any exception
  try {

    // Check quantity compatibility
    if (!checkQuantityCompatibility()) {
      THROW_EXCEPTION(QP_INPUT_ERROR_EXCEPTION, "QP solve unsuccessful. Input error.");
    }

    // Initialize minimum index and value
    int index = -1;
    double value = NONOPT_DOUBLE_INFINITY;

    // Loop through vector list
    for (int i = 0; i < (int)vector_list_.size(); i++) {

      // Compute objective value 0.5*g_i^T*W*g_i-b_i
      double t = 0.5 * matrix_->innerProductOfInverse(*vector_list_[i]) - vector_[i];

      // Check for minimum
      if (t < value) {

        // Update minimum's index and value
        index = i;
        value = t;

      }  // end if

    }  // end for

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

    // Solve hot
    solveQPHot(options, reporter);

  }  // end try

  // catch exceptions
  catch (QP_INPUT_ERROR_EXCEPTION& exec) {
    setStatus(QP_INPUT_ERROR);
  }

}  // end solveQP

// Solve hot
void QPSolverActiveSet::solveQPHot(const Options* options,
                                   const Reporter* reporter)
{

  // Initialize values
  setStatus(QP_UNSET);
  iteration_count_ = 0;
  kkt_error_ = -NONOPT_DOUBLE_INFINITY;

  // Declare system vectors
  double* inner_solution_3 = new double[system_solution_length_];
  double* inner_solution_ls = new double[system_solution_length_];
  double* inner_solution_trial = new double[system_solution_length_];
  double* new_system_vector = new double[system_solution_length_];
  double* right_hand_side = new double[system_solution_length_];

  // try to solve QP with hot start, terminate on any exception
  try {

    // Check quantity compatibility
    if (!checkQuantityCompatibility()) {
      THROW_EXCEPTION(QP_INPUT_ERROR_EXCEPTION, "QP solve unsuccessful. Input error.");
    }

    // Print message
    reporter->printf(R_QP, R_PER_ITERATION, "\n");
    reporter->printf(R_QP, R_PER_INNER_ITERATION, "Entering main iteration loop\n");

    // Set iteration limit
    int iteration_limit = fmax(iteration_limit_minimum_, fmin(pow((int)vector_list_.size(), 2) + pow(gamma_length_, 3), iteration_limit_maximum_));

    // Iteration loop
    while (true) {

      // Print message
      if (iteration_count_ % 20 == 0) {
        reporter->printf(R_QP, R_PER_ITERATION,
                         "=======================================================\n"
                         "  Iter.    |S|     |P|     |N|    min(KKT)  Set changes\n"
                         "=======================================================\n");
      }
      reporter->printf(R_QP, R_PER_ITERATION,
                       " %6d  %6d  %6d  %6d",
                       iteration_count_,
                       (int)omega_positive_.size(),
                       (int)gamma_positive_.size(),
                       (int)gamma_negative_.size());
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "=========================\n"
                       "Starting iteration %6d\n"
                       "=========================\n",
                       iteration_count_);
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "omega_positive (%6d elements):",
                       (int)omega_positive_.size());
      for (int i = 0; i < (int)omega_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "gamma_positive (%6d elements):",
                       (int)gamma_positive_.size());
      for (int i = 0; i < (int)gamma_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "gamma_negative (%6d elements):",
                       (int)gamma_negative_.size());
      for (int i = 0; i < (int)gamma_negative_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

      // Update best solution
      bool real_solution = updateBestSolution();

      // Evaluate dual vectors
      evaluateDualVectors();

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
        kkt_residual_omega[i] = multiplier_ - vector_[i] - vector_list_[i]->innerProduct(dual_step_);
      }

      // Zero-out omega's KKT error components for positive set
      for (int i = 0; i < (int)omega_positive_.size(); i++) {
        kkt_residual_omega[omega_positive_[i]] = 0.0;
      }

      // Evaluate gamma's KKT error components for positive side
      for (int i = 0; i < gamma_length_; i++) {
        kkt_residual_gamma_positive[i] = scalar_ - dual_step_.values()[i];
      }

      // Zero-out gamma's KKT error components for positive side and positive set
      for (int i = 0; i < (int)gamma_positive_.size(); i++) {
        kkt_residual_gamma_positive[gamma_positive_[i]] = 0.0;
      }

      // Evaluate gamma's KKT error components for negative side
      for (int i = 0; i < gamma_length_; i++) {
        kkt_residual_gamma_negative[i] = scalar_ + dual_step_.values()[i];
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
      }  // end if
      else {
        if (kkt_residual_omega_minimum <= kkt_residual_gamma_positive_minimum) {
          if (kkt_residual_omega_minimum <= kkt_residual_gamma_negative_minimum) {
            kkt_error_ = kkt_residual_omega_minimum;
            kkt_residual_minimum_set = 1;
            kkt_residual_minimum_index = kkt_residual_omega_minimum_index;
          }  // end if
          else {
            kkt_error_ = kkt_residual_gamma_negative_minimum;
            kkt_residual_minimum_set = 3;
            kkt_residual_minimum_index = kkt_residual_gamma_negative_minimum_index;
          }  // end else
        }    // end if
        else {
          if (kkt_residual_gamma_positive_minimum <= kkt_residual_gamma_negative_minimum) {
            kkt_error_ = kkt_residual_gamma_positive_minimum;
            kkt_residual_minimum_set = 2;
            kkt_residual_minimum_index = kkt_residual_gamma_positive_minimum_index;
          }  // end if
          else {
            kkt_error_ = kkt_residual_gamma_negative_minimum;
            kkt_residual_minimum_set = 3;
            kkt_residual_minimum_index = kkt_residual_gamma_negative_minimum_index;
          }  // end else
        }    // end else
      }      // end else

      // Print message
      reporter->printf(R_QP, R_PER_ITERATION, "  %+.2e", kkt_error_);
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "Set of minimum KKT element is %d\n"
                       "Index of minimum KKT element is %d with value %+.16e\n",
                       kkt_residual_minimum_set,
                       kkt_residual_minimum_index,
                       kkt_error_);

      // Check for successful solve
      if (kkt_error_ >= -kkt_tolerance_ ||
          (allow_inexact_termination_ && inexactTerminationCondition(reporter))) {
        THROW_EXCEPTION(QP_SUCCESS_EXCEPTION, "QP solve successful.");
      }

      // Check for iteration limit
      if (iteration_count_ >= iteration_limit) {
        THROW_EXCEPTION(QP_ITERATION_LIMIT_EXCEPTION, "QP solve unsuccessful. Iteration limit reached.");
      }

      // Print message
      reporter->printf(R_QP, R_PER_ITERATION, "  %d", kkt_residual_minimum_set);
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "Index %d will be added to index set %d\n",
                       kkt_residual_minimum_index,
                       kkt_residual_minimum_set);

      // Increment iteration counter
      iteration_count_++;

      // Evaluate new system vector
      evaluateSystemVector(kkt_residual_minimum_set, kkt_residual_minimum_index, new_system_vector);

      // Print message
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Solving intermediate system for least squares\n");

      // Solve intermediate system for least squares
      solveSystemTranspose(new_system_vector, inner_solution_3);

      // Declare new diagonal value
      double new_diagonal_squared;

      // Check which set is being updated
      if (kkt_residual_minimum_set == 1) {
        new_diagonal_squared = 1.0 + matrix_->innerProductOfInverse(*vector_list_[kkt_residual_minimum_index]);
      }
      else {
        new_diagonal_squared = matrix_->elementOfInverse(kkt_residual_minimum_index, kkt_residual_minimum_index);
      }

      // Set inputs for blas
      int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
      int increment1 = 1;

      // Compute new diagonal for factor (squared)
      double rho2 = fmax(0.0, new_diagonal_squared - ddot_(&length, inner_solution_3, &increment1, inner_solution_3, &increment1));

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
        solveSystem(inner_solution_3, inner_solution_ls);

        // Set inputs for blas
        int length = (int)omega_positive_.size();
        int increment0 = 0;
        int increment1 = 1;
        double value = 1.0;

        // Declare residual value
        double residual = ddot_(&length, &value, &increment0, inner_solution_ls, &increment1);

        // Check for sufficiently negative value
        if (residual - 1.0 < -linear_independence_tolerance_) {
          linear_independence_flag = true;
        }

      }  // end else

      // Print message
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "Linear independence check yields %d\n",
                       linear_independence_flag);

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
            if (inner_solution_ls[i] > 0.0) {
              double temporary_scalar = system_solution_[i] / inner_solution_ls[i];
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 1;
              }  // end if
            }    // end if
          }      // end for
        }        // end if
        else {
          for (int i = 0; i < (int)gamma_positive_.size(); i++) {
            if (inner_solution_ls[(int)omega_positive_.size() + i] > 0.0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + i] / inner_solution_ls[(int)omega_positive_.size() + i];
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 2;
              }  // end if
            }    // end if
          }      // end for
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            if (inner_solution_ls[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] > 0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] / inner_solution_ls[(int)omega_positive_.size() + (int)gamma_positive_.size() + i];
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 3;
              }  // end if
            }    // end if
          }      // end for
        }        // end else

        // Update augmentation value
        augmentation_value = delete_value;

        // Print message
        reporter->printf(R_QP, R_PER_ITERATION, "  -%d", delete_set);
        if (delete_set == 1) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "Index %d will be deleted from omega's positive set\n",
                           omega_positive_[delete_index]);
        }
        else if (delete_set == 2) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "Index %d will be deleted from gamma's positive set\n",
                           gamma_positive_[delete_index]);
        }
        else if (delete_set == 3) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "Index %d will be deleted from gamma's negative set\n",
                           gamma_negative_[delete_index]);
        }
        else {
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "Uh oh!  Search for element to delete failed!\n");
        }

        // Check if search for element to delete was successful
        if (delete_set > 0) {

          // Set inputs for blas
          int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
          double value = -delete_value;
          int increment = 1;

          // Update solution
          daxpy_(&length, &value, inner_solution_ls, &increment, system_solution_, &increment);

          // Perform set deletion
          setDelete(reporter, delete_set, delete_index, inner_solution_1_, inner_solution_2_);

        }  // end if

        // Print sets
        reporter->printf(R_QP, R_PER_INNER_ITERATION,
                         "omega_positive (%6d elements):",
                         (int)omega_positive_.size());
        for (int i = 0; i < (int)omega_positive_.size(); i++) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
        }
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
        reporter->printf(R_QP, R_PER_INNER_ITERATION,
                         "gamma_positive (%6d elements):",
                         (int)gamma_positive_.size());
        for (int i = 0; i < (int)gamma_positive_.size(); i++) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
        }
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
        reporter->printf(R_QP, R_PER_INNER_ITERATION,
                         "gamma_negative (%6d elements):",
                         (int)gamma_negative_.size());
        for (int i = 0; i < (int)gamma_negative_.size(); i++) {
          reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
        }
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

      }  // end if (linear_independence_flag == false)

      // Print message
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "Performing set augmentation\n");

      // Perform set augmentation
      bool augment_success = setAugment(reporter, kkt_residual_minimum_set, kkt_residual_minimum_index, new_system_vector, inner_solution_1_, inner_solution_2_, augmentation_value);

      // Print sets
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "omega_positive (%6d elements):",
                       (int)omega_positive_.size());
      for (int i = 0; i < (int)omega_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "gamma_positive (%6d elements):",
                       (int)gamma_positive_.size());
      for (int i = 0; i < (int)gamma_positive_.size(); i++) {
        reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
      }
      reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "gamma_negative (%6d elements):",
                       (int)gamma_negative_.size());
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
        reporter->printf(R_QP, R_PER_INNER_ITERATION,
                         "Re-doing Cholesky factorization from scratch!\n");

        // Try to compute from scratch
        choleskyFromScratch(reporter);

      }  // end if

      // Compute multiplier
      evaluateDualMultiplier(inner_solution_1_, inner_solution_2_);

      // Set inner iteration limit
      int inner_iteration_limit = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();

      // Subproblem solution loop
      for (int inner_iteration_count = 0; inner_iteration_count < inner_iteration_limit; inner_iteration_count++) {

        // Print message
        reporter->printf(R_QP, R_PER_INNER_ITERATION, "Solving subproblem\n");

        // Set inputs for blas
        int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
        int increment = 1;
        double value = 1.0 - multiplier_;

        // Set right-hand side vector
        dcopy_(&length, inner_solution_2_, &increment, right_hand_side, &increment);
        daxpy_(&length, &value, inner_solution_1_, &increment, right_hand_side, &increment);

        // Solve subproblem
        solveSystem(right_hand_side, inner_solution_trial);

        // Declare boolean for correct signs
        bool correct_signs = true;

        // Check signs of subproblem solution elements
        if (correct_signs) {
          for (int i = 0; i < (int)omega_positive_.size() + (int)gamma_positive_.size(); i++) {
            if (inner_solution_trial[i] <= 0.0) {
              correct_signs = false;
              break;
            }  // end if
          }    // end for
        }      // end if
        if (correct_signs) {
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            if (inner_solution_trial[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] >= 0.0) {
              correct_signs = false;
              break;
            }  // end if
          }    // end for
        }      // end if

        // Check feasibility of subproblem solution
        if (correct_signs) {

          // Print message
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "Feasible subproblem solution! Breaking loop\n");

          // Replace current solution
          dcopy_(&length, inner_solution_trial, &increment, system_solution_, &increment);

          // Break subproblem solution loop
          break;

        }  // end if
        else {

          // Print message
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "Infeasible subproblem solution! Looking for index to delete\n");

          // Initialize minimum index and values
          int delete_index = -1;
          double delete_value = NONOPT_DOUBLE_INFINITY;
          int delete_set = -1;

          // Compute index with minimum value
          for (int i = 0; i < (int)omega_positive_.size(); i++) {
            if (inner_solution_trial[i] < 0.0) {
              double temporary_scalar = system_solution_[i] / (system_solution_[i] - inner_solution_trial[i]);
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 1;
              }  // end if
            }    // end if
          }      // end for
          for (int i = 0; i < (int)gamma_positive_.size(); i++) {
            if (inner_solution_trial[(int)omega_positive_.size() + i] < 0.0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + i] / (system_solution_[(int)omega_positive_.size() + i] - inner_solution_trial[(int)omega_positive_.size() + i]);
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 2;
              }  // end if
            }    // end if
          }      // end for
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            if (inner_solution_trial[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] > 0.0) {
              double temporary_scalar = system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] / (system_solution_[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] - inner_solution_trial[(int)omega_positive_.size() + (int)gamma_positive_.size() + i]);
              if (temporary_scalar < delete_value) {
                delete_index = i;
                delete_value = temporary_scalar;
                delete_set = 3;
              }  // end if
            }    // end if
          }      // end for

          // Print message
          reporter->printf(R_QP, R_PER_ITERATION, "  %+d", -delete_set);
          if (delete_set == 1) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION,
                             "Index %d will be deleted from omega's positive set\n",
                             omega_positive_[delete_index]);
          }
          else if (delete_set == 2) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION,
                             "Index %d will be deleted from gamma's positive set\n",
                             gamma_positive_[delete_index]);
          }
          else if (delete_set == 3) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION,
                             "Index %d will be deleted from gamma's negative set\n",
                             gamma_negative_[delete_index]);
          }
          else {
            reporter->printf(R_QP, R_PER_INNER_ITERATION,
                             "Uh oh!  Search for element to delete failed!\n");
          }

          // Check if search for element to delete was successful
          if (delete_set > 0) {

            // Update value
            delete_value = fmin(1.0, delete_value);

            // Set inputs for blas
            int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
            double value = 1.0 - delete_value;
            int increment = 1;

            // Scale solution
            dscal_(&length, &value, system_solution_, &increment);

            // Update solution
            daxpy_(&length, &delete_value, inner_solution_trial, &increment, system_solution_, &increment);

            // Perform set deletion
            setDelete(reporter, delete_set, delete_index, inner_solution_1_, inner_solution_2_);

            // Evaluate multiplier
            evaluateDualMultiplier(inner_solution_1_, inner_solution_2_);

          }  // end if

          // Print sets
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "omega_positive (%6d elements):",
                           (int)omega_positive_.size());
          for (int i = 0; i < (int)omega_positive_.size(); i++) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", omega_positive_[i]);
          }
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "gamma_positive (%6d elements):",
                           (int)gamma_positive_.size());
          for (int i = 0; i < (int)gamma_positive_.size(); i++) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_positive_[i]);
          }
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");
          reporter->printf(R_QP, R_PER_INNER_ITERATION,
                           "gamma_negative (%6d elements):",
                           (int)gamma_negative_.size());
          for (int i = 0; i < (int)gamma_negative_.size(); i++) {
            reporter->printf(R_QP, R_PER_INNER_ITERATION, " %d", gamma_negative_[i]);
          }
          reporter->printf(R_QP, R_PER_INNER_ITERATION, "\n");

        }  // end else (correct_signs)

      }  // end for (subproblem solution loop)

      // Print new line
      reporter->printf(R_QP, R_PER_ITERATION, "\n");

    }  // end of main iteration loop

  }  // end try

  // catch exceptions
  catch (QP_SUCCESS_EXCEPTION& exec) {
    setStatus(QP_SUCCESS);
  }
  catch (QP_FACTORIZATION_ERROR_EXCEPTION& exec) {
    setStatus(QP_FACTORIZATION_ERROR);
  }
  catch (QP_INPUT_ERROR_EXCEPTION& exec) {
    setStatus(QP_INPUT_ERROR);
  }
  catch (QP_ITERATION_LIMIT_EXCEPTION& exec) {
    setStatus(QP_ITERATION_LIMIT);
  }
  catch (QP_NAN_ERROR_EXCEPTION& exec) {
    setStatus(QP_NAN_ERROR);
  }

  // Print new line
  reporter->printf(R_QP, R_PER_ITERATION, "\n");

  // Print finalizing
  reporter->printf(R_QP, R_PER_INNER_ITERATION, "Finalizing solution\n");

  // Finalize solution
  finalizeSolution();

  // Delete array
  if (inner_solution_3 != nullptr) {
    delete[] inner_solution_3;
  }
  if (inner_solution_ls != nullptr) {
    delete[] inner_solution_ls;
  }
  if (inner_solution_trial != nullptr) {
    delete[] inner_solution_trial;
  }
  if (new_system_vector != nullptr) {
    delete[] new_system_vector;
  }
  if (right_hand_side != nullptr) {
    delete[] right_hand_side;
  }

}  // end solveQPHot

// Print method
void QPSolverActiveSet::printData(const Reporter* reporter)
{
  // Print matrix values
  reporter->printf(R_QP, R_BASIC, "\nMATRIX:\n");
  for (int i = 0; i < gamma_length_; i++) {
    for (int j = 0; j < gamma_length_; j++) {
      reporter->printf(R_QP, R_BASIC, " %+23.16e", matrix_->elementOfInverse(i, j));
    }
    reporter->printf(R_QP, R_BASIC, "\n");
  }  // end for

  // Print vector list
  reporter->printf(R_QP, R_BASIC, "VECTOR LIST:\n");
  for (int j = 0; j < gamma_length_; j++) {
    for (int i = 0; i < (int)vector_list_.size(); i++) {
      reporter->printf(R_QP, R_BASIC, " %+23.16e", vector_list_[i]->values()[j]);
    }
    reporter->printf(R_QP, R_BASIC, "\n");
  }  // end for

  // Print vector
  reporter->printf(R_QP, R_BASIC, "VECTOR:\n");
  for (int i = 0; i < (int)vector_.size(); i++) {
    reporter->printf(R_QP, R_BASIC, " %+23.16e\n", vector_[i]);
  }

  // Print scalar
  reporter->printf(R_QP, R_BASIC, "SCALAR:\n");
  reporter->printf(R_QP, R_BASIC, " %+23.16e\n", scalar_);

}  // end printData

// Check quantities
bool QPSolverActiveSet::checkQuantityCompatibility()
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

  }  // end for

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

}  // end checkQuantityCompatibility

// Update best solution
bool QPSolverActiveSet::updateBestSolution()
{

  // Initialize boolean
  bool real_solution = true;

  // Loop through solution
  for (int i = 0; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size(); i++) {
    if (isnan(system_solution_[i])) {
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

  }  // end if

  // Return
  return real_solution;

}  // end updateBestSolution

// Cholesky augmentation
bool QPSolverActiveSet::choleskyAugment(double system_vector[],
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

  // "Add" zero values to R by shifting values
  for (int i = 0; i < length; i++) {
    for (int j = length - 1; j > index; j--) {
      factor_[i * system_solution_length_ + j] = factor_[i * system_solution_length_ + j - 1];
    }
    factor_[i * system_solution_length_ + index] = 0.0;
  }  // end for
  for (int j = 0; j < length; j++) {
    for (int i = length - 1; i > index; i--) {
      factor_[i * system_solution_length_ + j] = factor_[(i - 1) * system_solution_length_ + j];
    }
    factor_[index * system_solution_length_ + j] = 0.0;
  }  // end for

  // "Add" zero values to solutions by shifting values
  for (int i = length - 1; i > index; i--) {
    solution1[i] = solution1[i - 1];
    solution2[i] = solution2[i - 1];
  }  // end for
  solution1[index] = 0.0;
  solution2[index] = 0.0;

  // Update Cholesky factor (first part)
  for (int i = 0; i < index; i++) {
    factor_[i * system_solution_length_ + index] = system_vector[i] / factor_[i * system_solution_length_ + i];
    for (int j = i + 1; j < index; j++) {
      system_vector[j] = system_vector[j] - factor_[i * system_solution_length_ + j] * factor_[i * system_solution_length_ + index];
    }
  }  // end for

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
  }  // end else

  // Update Cholesky factor (third part)
  for (int i = index + 1; i < length; i++) {
    temporary_scalar = 0.0;
    for (int j = 0; j < index; j++) {
      temporary_scalar = temporary_scalar + factor_[j * system_solution_length_ + i] * factor_[j * system_solution_length_ + index];
    }
    factor_[index * system_solution_length_ + i] = (system_vector[i] - temporary_scalar) / factor_[index * system_solution_length_ + index];
  }  // end for

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

  }  // end for

  // Delete temporary vector
  if (temporary_vector != nullptr) {
    delete[] temporary_vector;
  }

  // Return
  return success_without_factorization_error;

}  // end choleskyAugment

// Cholesky deletion
void QPSolverActiveSet::choleskyDelete(int index,
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
  }  // end for

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

    }  // end for

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

  }  // end for

}  // end choleskyDelete

// Cholesky factorization from scratch
void QPSolverActiveSet::choleskyFromScratch(const Reporter* reporter)
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
  std::vector<std::shared_ptr<Vector> > WG;

  // Compute WG matrix
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    std::shared_ptr<Vector> product(new Vector(gamma_length_));
    matrix_->matrixVectorProductOfInverse(*vector_list_[omega_positive_[i]], *product);
    WG.push_back(product);
  }  // end for

  // Compute (1,1)-block (upper triangle)
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    for (int j = i; j < (int)omega_positive_.size(); j++) {
      matrix[i * size + j] = 1.0 + vector_list_[omega_positive_[i]]->innerProduct(*WG[j]);
    }
  }  // end for

  // Compute (1,2)-block
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    for (int j = 0; j < (int)gamma_positive_.size(); j++) {
      matrix[i * size + (int)omega_positive_.size() + j] = WG[i]->values()[gamma_positive_[j]];
    }
  }  // end for

  // Compute (1,3)-block
  for (int i = 0; i < (int)omega_positive_.size(); i++) {
    for (int j = 0; j < (int)gamma_negative_.size(); j++) {
      matrix[i * size + (int)omega_positive_.size() + (int)gamma_positive_.size() + j] = WG[i]->values()[gamma_negative_[j]];
    }
  }  // end for

  // Compute (2,2)-block (upper triangle)
  for (int i = 0; i < (int)gamma_positive_.size(); i++) {
    for (int j = i; j < (int)gamma_positive_.size(); j++) {
      matrix[((int)omega_positive_.size() + i) * size + ((int)omega_positive_.size() + j)] = matrix_->elementOfInverse(gamma_positive_[i], gamma_positive_[j]);
    }
  }  // end for

  // Compute (2,3)-block (upper triangle)
  for (int i = 0; i < (int)gamma_positive_.size(); i++) {
    for (int j = 0; j < (int)gamma_negative_.size(); j++) {
      matrix[((int)omega_positive_.size() + i) * size + ((int)omega_positive_.size() + (int)gamma_positive_.size() + j)] = matrix_->elementOfInverse(gamma_positive_[i], gamma_negative_[j]);
    }
  }  // end for

  // Compute (3,3)-block (upper triangle)
  for (int i = 0; i < (int)gamma_negative_.size(); i++) {
    for (int j = i; j < (int)gamma_negative_.size(); j++) {
      matrix[((int)omega_positive_.size() + (int)gamma_positive_.size() + i) * size + ((int)omega_positive_.size() + (int)gamma_positive_.size() + j)] = matrix_->elementOfInverse(gamma_negative_[i], gamma_negative_[j]);
    }
  }  // end for

  // Set inputs for blas
  char upper_lower = 'L';  // Use "lower" since Fortran uses column-major ordering
  int flag = 0;

  // Compute factorization
  dpotf2_(&upper_lower, &size, matrix, &size, &flag);

  // Check tolerance on diagonal
  for (int i = 0; i < size; i++) {
    if (matrix[i * size + i] < cholesky_tolerance_) {
      reporter->printf(R_QP, R_PER_INNER_ITERATION,
                       "ARGH! Replacing Cholesky diagonal of %+23.16e with %+23.16e.\n",
                       matrix[i * size + i],
                       cholesky_tolerance_);
      matrix[i * size + i] = cholesky_tolerance_;
    }  // end if
  }    // end for

  // Reset factorization
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      factor_[i * system_solution_length_ + j] = matrix[i * size + j];
    }
  }  // end for

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
  }
  if (right_hand_side != nullptr) {
    delete[] right_hand_side;
  }

}  // end choleskyFromScratch

// Evaluate dual vectors
void QPSolverActiveSet::evaluateDualVectors()
{

  // Zero-out vectors
  combination_.scale(0.0);
  combination_translated_.scale(0.0);
  dual_step_.scale(0.0);
  dual_step_feasible_.scale(0.0);

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
  matrix_->matrixVectorProductOfInverse(combination_translated_, dual_step_);
  dual_step_.scale(-1.0);

  // Compute feasible dual step by projection
  dual_step_feasible_.copy(dual_step_);
  dual_step_projection_scalar_ = fmin(1.0, scalar_ / dual_step_.normInf());
  dual_step_feasible_.scale(dual_step_projection_scalar_);

}  // end evaluateDualVectors

// Evaluate dual multiplier
void QPSolverActiveSet::evaluateDualMultiplier(double solution1[],
                                               double solution2[])
{

  // Set inputs for blas
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  int increment = 1;

  // Evaluate inner products
  double solution1_norm_squared = ddot_(&length, solution1, &increment, solution1, &increment);
  double solution1_solution2 = ddot_(&length, solution1, &increment, solution2, &increment);

  // Evaluate multiplier
  multiplier_ = (solution1_norm_squared + solution1_solution2 - 1.0) / solution1_norm_squared;

}  // end evaluateDualMultiplier

// Evaluate system vector
void QPSolverActiveSet::evaluateSystemVector(int set,
                                             int index,
                                             double system_vector[])
{

  // Set inputs for blas
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  double value = 0.0;
  int increment0 = 0;
  int increment1 = 1;

  // Initialize values
  dcopy_(&length, &value, &increment0, system_vector, &increment1);

  // Check set
  if (set == 1) {

    // Declare temporary vector
    double* temporary_vector = new double[gamma_length_];

    // Evaluate temporary vector
    for (int i = 0; i < gamma_length_; i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(i, col);
      temporary_vector[i] = vector_list_[index]->innerProduct(col);
    }  // end for

    // Set "omega" values, i.e.,
    for (int i = 0; i < (int)omega_positive_.size(); i++) {
      system_vector[i] = system_vector[i] + ddot_(&gamma_length_, vector_list_[omega_positive_[i]]->values(), &increment1, temporary_vector, &increment1) + 1.0;
    }

    // Set "gamma positive" values, i.e.,
    for (int i = 0; i < (int)gamma_positive_.size(); i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(gamma_positive_[i], col);
      system_vector[(int)omega_positive_.size() + i] = vector_list_[index]->innerProduct(col);
    }  // end for

    // Set "gamma negative" values, i.e.,
    for (int i = 0; i < (int)gamma_negative_.size(); i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(gamma_negative_[i], col);
      system_vector[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] = vector_list_[index]->innerProduct(col);
    }  // end for

    // Delete temporary vector
    if (temporary_vector != nullptr) {
      delete[] temporary_vector;
    }

  }  // end if

  else {

    // Set "omega" values, i.e., ith value = G(:,omega_positive_[i])'*W(:,kkt_residual_minimum_index)
    for (int i = 0; i < (int)omega_positive_.size(); i++) {
      Vector col(gamma_length_, 0.0);
      matrix_->columnOfInverse(index, col);
      system_vector[i] = vector_list_[omega_positive_[i]]->innerProduct(col);
    }  // end for

    // Set "gamma positive" values, i.e., W(gamma_positive_[i],kkt_residual_minimum_index)
    for (int i = 0; i < (int)gamma_positive_.size(); i++) {
      system_vector[(int)omega_positive_.size() + i] = matrix_->elementOfInverse(gamma_positive_[i], index);
    }

    // Set "gamma negative" values, i.e., W(gamma_positive_[i],kkt_residual_minimum_index)
    for (int i = 0; i < (int)gamma_negative_.size(); i++) {
      system_vector[(int)omega_positive_.size() + (int)gamma_positive_.size() + i] = matrix_->elementOfInverse(gamma_negative_[i], index);
    }

  }  // end else

}  // end evaluateSystemVector

// Finalize solution
void QPSolverActiveSet::finalizeSolution()
{

  // Set sizes
  omega_.setLength((int)vector_.size());
  gamma_.setLength(gamma_length_);

  // Set final omega values
  for (int i = 0; i < (int)omega_positive_best_.size(); i++) {
    omega_.set(omega_positive_best_[i], system_solution_best_[i]);
  }

  // Set final gamma values
  for (int i = 0; i < (int)gamma_positive_best_.size(); i++) {
    gamma_.set(gamma_positive_best_[i], system_solution_best_[(int)omega_positive_best_.size() + i]);
  }
  for (int i = 0; i < (int)gamma_negative_best_.size(); i++) {
    gamma_.set(gamma_negative_best_[i], system_solution_best_[(int)omega_positive_best_.size() + (int)gamma_positive_best_.size() + i]);
  }

  // Set solution to best
  omega_positive_ = omega_positive_best_;
  gamma_positive_ = gamma_positive_best_;
  gamma_negative_ = gamma_negative_best_;

  // Update solution
  for (int i = 0; i < (int)omega_positive_best_.size() + (int)gamma_positive_best_.size() + (int)gamma_negative_best_.size(); i++) {
    system_solution_[i] = system_solution_best_[i];
  }

}  // end finalizeSolution

// Set augment
bool QPSolverActiveSet::setAugment(const Reporter* reporter,
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

  }  // end if

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

  }  // end else if

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

  }  // end else

  // Return
  return success_without_factorization_error;

}  // end setAugment

// Set delete
void QPSolverActiveSet::setDelete(const Reporter* reporter,
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

  }  // end if

  else if (set == 2) {

    // Remove solution element
    for (int i = (int)omega_positive_.size() + index; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() - 1; i++) {
      system_solution_[i] = system_solution_[i + 1];
    }

    // Remove element from set
    gamma_positive_.erase(gamma_positive_.begin() + index);

    // Delete row/column from Cholesky factor
    choleskyDelete((int)omega_positive_.size() + index, solution1, solution2);

  }  // end else if

  else {

    // Remove solution element
    for (int i = (int)omega_positive_.size() + (int)gamma_positive_.size() + index; i < (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size() - 1; i++) {
      system_solution_[i] = system_solution_[i + 1];
    }

    // Remove element from set
    gamma_negative_.erase(gamma_negative_.begin() + index);

    // Delete row/column from Cholesky factor
    choleskyDelete((int)omega_positive_.size() + (int)gamma_positive_.size() + index, solution1, solution2);

  }  // end else

}  // end setDelete

// Solve linear system
void QPSolverActiveSet::solveSystem(double right_hand_side[],
                                    double solution[])
{

  // Set inputs for blas
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  char upper_lower = 'L';  // Use "lower" since Fortran uses column-major ordering
  char transpose = 'T';    // Use "transpose" since Fortran uses column-major ordering
  char diagonal = 'N';
  int increment1 = 1;
  int incrementn = system_solution_length_;

  // Copy right_hand_side to solution
  dcopy_(&length, right_hand_side, &increment1, solution, &increment1);

  // Solve system
  dtrsv_(&upper_lower, &transpose, &diagonal, &length, factor_, &incrementn, solution, &increment1);

}  // end solveSystem

// Solve triangular system with transpose
void QPSolverActiveSet::solveSystemTranspose(double right_hand_side[],
                                             double solution[])
{

  // Set inputs for blas
  int length = (int)omega_positive_.size() + (int)gamma_positive_.size() + (int)gamma_negative_.size();
  char upper_lower = 'L';  // Use "lower" since Fortran uses column-major ordering
  char transpose = 'N';    // Use "not transpose" since Fortran uses column-major ordering
  char diagonal = 'N';
  int increment1 = 1;
  int incrementn = system_solution_length_;

  // Copy right_hand_side to solution
  dcopy_(&length, right_hand_side, &increment1, solution, &increment1);

  // Solve system
  dtrsv_(&upper_lower, &transpose, &diagonal, &length, factor_, &incrementn, solution, &increment1);

}  // end solveSystemTranspose

}  // namespace NonOpt
