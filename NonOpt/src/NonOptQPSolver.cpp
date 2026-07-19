// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptBLASLAPACK.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptQPSolver.hpp"
#include "NonOptSymmetricMatrix.hpp"

namespace NonOpt
{

// Get combination norm
double QPSolver::combinationNormInf()
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
double QPSolver::combinationTranslatedNormInf()
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
double QPSolver::combinationTranslatedNorm2Squared()
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
double QPSolver::dualObjectiveQuadraticValue()
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
double QPSolver::dualObjectiveQuadraticValueScaled()
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Set quadratic value
  double dual_objective_quadratic_value_scaled = -ddot_(&length, primal_solution_.values(), &increment, combination_translated_.values(), &increment);

  // Scale so value corresponds to feasible dual step
  dual_objective_quadratic_value_scaled *= primal_solution_projection_scalar_ * primal_solution_projection_scalar_;

  // Set maximum
  if (std::isnan(dual_objective_quadratic_value_scaled) || dual_objective_quadratic_value_scaled > NONOPT_DOUBLE_INFINITY) {
    dual_objective_quadratic_value_scaled = NONOPT_DOUBLE_INFINITY;
  }

  // Return
  return dual_objective_quadratic_value_scaled;

} // end dualObjectiveQuadraticValueScaled

// Get dual solution
void QPSolver::dualSolution(double omega[], double gamma[])
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
void QPSolver::dualSolutionOmega(double omega[])
{

  // Set inputs for BLASLAPACK
  int length = (int)vector_.size();
  int increment = 1;

  // Copy values
  dcopy_(&length, omega_.values(), &increment, omega, &increment);

} // end dualSolutionOmega

// Get primal solution
void QPSolver::primalSolution(double d[])
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, primal_solution_.values(), &increment, d, &increment);

} // end primalSolution

// Get primal solution inf-norm
double QPSolver::primalSolutionNormInf()
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
double QPSolver::primalSolutionNorm2Squared()
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

// Get feasible primal solution
void QPSolver::primalSolutionFeasible(double d_feasible[])
{

  // Set inputs for BLASLAPACK
  int length = gamma_length_;
  int increment = 1;

  // Copy values
  dcopy_(&length, primal_solution_feasible_.values(), &increment, d_feasible, &increment);

} // end primalSolutionFeasible

// Get feasible primal solution norm
double QPSolver::primalSolutionFeasibleNormInf()
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

// Sanity check
bool QPSolver::checkQuantityCompatibility()
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

} // namespace NonOpt
