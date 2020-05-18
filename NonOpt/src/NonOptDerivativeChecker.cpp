// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptDerivativeChecker.hpp"
#include "NonOptEnumerations.hpp"

namespace NonOpt
{

// Check derivatives
void DerivativeChecker::checkDerivatives(const Reporter *reporter,
                                         const std::shared_ptr<Point> point) const
{

  // Print introductory message
  reporter->printf(R_NL, R_PER_ITERATION, "Checking first-order derivatives (increment = %+23.16e):", increment_);

  // Copy point's vector
  std::shared_ptr<Vector> perturbation_plus = point->vector()->makeNewCopy();
  std::shared_ptr<Vector> perturbation_minus = point->vector()->makeNewCopy();

  // Initialize counter of poor derivatives
  int poor_count = 0;

  // Loop over entries
  for (int i = 0; i < point->vector()->length(); i++) {

    // Add perturbations to element
    perturbation_plus->valuesModifiable()[i] += increment_;
    perturbation_minus->valuesModifiable()[i] -= increment_;

    // Declare values
    double objective_plus, objective_minus;

    // Evaluate objective at perturbed points
    bool evaluated_plus = point->problem()->evaluateObjective(point->vector()->length(), perturbation_plus->values(), objective_plus);
    bool evaluated_minus = point->problem()->evaluateObjective(point->vector()->length(), perturbation_minus->values(), objective_minus);

    // Apply objective scaling
    objective_plus *= point->scale();
    objective_minus *= point->scale();

    // Check evaluation errors
    if (!evaluated_plus || !evaluated_minus) {

      // Increment poor counter
      poor_count++;

      // Print message
      reporter->printf(R_NL, R_PER_ITERATION, "\n  g[%8d] = evaluation error", i);

    }  // end if
    else {

      // Evaluate errors
      double finite_difference_derivative_forward = (objective_plus - point->objective())/increment_;
      double finite_difference_derivative_backward = (point->objective() - objective_minus)/increment_;
      double finite_difference_derivative_central = (objective_plus - objective_minus)/(2*increment_);
      double error_forward = fabs(finite_difference_derivative_forward - point->gradient()->values()[i]) / fmax(finite_difference_derivative_central,1.0);
      double error_backward = fabs(finite_difference_derivative_backward - point->gradient()->values()[i]) / fmax(finite_difference_derivative_central,1.0);
      double error_central = fabs(finite_difference_derivative_central - point->gradient()->values()[i]) / fmax(finite_difference_derivative_central,1.0);

      // Set best error
      double finite_difference_derivative, error;
      if (error_forward <= error_backward && error_forward <= error_central) {
        finite_difference_derivative = finite_difference_derivative_forward;
        error = error_forward;
      }
      else if (error_backward <= error_forward && error_backward <= error_central) {
        finite_difference_derivative = finite_difference_derivative_backward;
        error = error_backward;
      }
      else {
        finite_difference_derivative = finite_difference_derivative_central;
        error = error_central;
      }

      // Check error
      if (error > tolerance_) {

        // Increment poor counter
        poor_count++;

        // Print message
        reporter->printf(R_NL, R_PER_ITERATION, "\n  g[%8d]=%+23.16e != %+23.16e", i, point->gradient()->values()[i], finite_difference_derivative);

      }  // end if
    }    // end else

    // Remove perturbations from element
    perturbation_plus->valuesModifiable()[i] -= increment_;
    perturbation_minus->valuesModifiable()[i] += increment_;

  }  // end for

  // Check poor count
  if (poor_count == 0) {
    reporter->printf(R_NL, R_PER_ITERATION, " ACCURATE!\n");
  }
  else {
    reporter->printf(R_NL, R_PER_ITERATION, "\nErrors above due to nonsmoothness at point or bug in gradient computation.\n");
  }

}  // end checkDerivativesAtPoint

}  // namespace NonOpt
