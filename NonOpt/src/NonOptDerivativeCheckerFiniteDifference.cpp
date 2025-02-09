// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptDefinitions.hpp"
#include "NonOptDerivativeCheckerFiniteDifference.hpp"

namespace NonOpt
{

// Add options
void DerivativeCheckerFiniteDifference::addOptions(Options* options)
{

  // Add bool options
  options->addBoolOption("DEFD_check_derivatives",
                         false,
                         "Determines whether to check derivatives at iterates.\n"
                         "Default     : false.");

  // Add double options
  options->addDoubleOption("DEFD_increment",
                           1e-08,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Increment for derivative checker.\n"
                           "Default     : 1e-08.");
  options->addDoubleOption("DEFD_tolerance",
                           1e-04,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for derivative checker.\n"
                           "Default     : 1e-04.");

  // Add integer options

} // end addOptions

// Set options
void DerivativeCheckerFiniteDifference::setOptions(Options* options)
{

  // Read bool options
  options->valueAsBool("DEFD_check_derivatives", check_derivatives_);

  // Read double options
  options->valueAsDouble("DEFD_increment", increment_);
  options->valueAsDouble("DEFD_tolerance", tolerance_);

  // Read integer options

} // end setOptions

// Initialize
void DerivativeCheckerFiniteDifference::initialize(const Options* options,
                                                   Quantities* quantities,
                                                   const Reporter* reporter)
{
} // end initialize

// Check derivatives
void DerivativeCheckerFiniteDifference::checkDerivatives(const Options* options,
                                                         Quantities* quantities,
                                                         const Reporter* reporter)
{

  // Initialize status
  setStatus(DE_UNSET);

  // Check whether to check!
  if (check_derivatives_) {

    // Print message
    reporter->printf(R_NL, R_BASIC, "Checking first-order derivatives (increment = %+23.16e):", increment_);

    // Copy current iterate's vector
    std::shared_ptr<Vector> perturbation_plus = quantities->currentIterate()->vector()->makeNewCopy();
    std::shared_ptr<Vector> perturbation_minus = quantities->currentIterate()->vector()->makeNewCopy();

    // Initialize counter of poor derivatives
    int poor_count = 0;

    // Loop over entries
    for (int i = 0; i < quantities->currentIterate()->vector()->length(); i++) {

      // Add perturbations to element
      perturbation_plus->valuesModifiable()[i] += increment_;
      perturbation_minus->valuesModifiable()[i] -= increment_;

      // Declare values
      double objective_plus, objective_minus;

      // Evaluate objective at perturbed points
      bool evaluated_plus;
      bool evaluated_minus;
      if (quantities->evaluateFunctionWithGradient()) {
        evaluated_plus = quantities->currentIterate()->problem()->evaluateObjective(quantities->currentIterate()->vector()->length(), perturbation_plus->values(), objective_plus);
        evaluated_minus = quantities->currentIterate()->problem()->evaluateObjective(quantities->currentIterate()->vector()->length(), perturbation_minus->values(), objective_minus);
      }
      else {
        double* dummy = new double[quantities->currentIterate()->vector()->length()];
        evaluated_plus = quantities->currentIterate()->problem()->evaluateObjectiveAndGradient(quantities->currentIterate()->vector()->length(), perturbation_plus->values(), objective_plus, dummy);
        evaluated_minus = quantities->currentIterate()->problem()->evaluateObjectiveAndGradient(quantities->currentIterate()->vector()->length(), perturbation_minus->values(), objective_minus, dummy);
      }

      // Apply objective scaling
      objective_plus *= quantities->currentIterate()->scale();
      objective_minus *= quantities->currentIterate()->scale();

      // Check evaluation errors
      if (!evaluated_plus || !evaluated_minus) {

        // Increment poor counter
        poor_count++;

        // Print message
        reporter->printf(R_NL, R_BASIC, "\n  g[%8d] = evaluation error", i);

      } // end if
      else {

        // Evaluate errors
        double finite_difference_derivative_forward = (objective_plus - quantities->currentIterate()->objective()) / increment_;
        double finite_difference_derivative_backward = (quantities->currentIterate()->objective() - objective_minus) / increment_;
        double finite_difference_derivative_central = (objective_plus - objective_minus) / (2 * increment_);
        double error_forward = fabs(finite_difference_derivative_forward - quantities->currentIterate()->gradient()->values()[i]) / fmax(finite_difference_derivative_central, 1.0);
        double error_backward = fabs(finite_difference_derivative_backward - quantities->currentIterate()->gradient()->values()[i]) / fmax(finite_difference_derivative_central, 1.0);
        double error_central = fabs(finite_difference_derivative_central - quantities->currentIterate()->gradient()->values()[i]) / fmax(finite_difference_derivative_central, 1.0);

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
          reporter->printf(R_NL, R_BASIC, "\n  g[%8d]=%+23.16e != %+23.16e", i, quantities->currentIterate()->gradient()->values()[i], finite_difference_derivative);

        } // end if
      }   // end else

      // Remove perturbations from element
      perturbation_plus->valuesModifiable()[i] -= increment_;
      perturbation_minus->valuesModifiable()[i] += increment_;

    } // end for

    // Check poor count
    if (poor_count == 0) {
      reporter->printf(R_NL, R_BASIC, " ACCURATE!\n");
    }
    else {
      reporter->printf(R_NL, R_BASIC, "\nErrors above due to nonsmoothness at point or bug in gradient computation.\n");
    }

  } // end if

  // Update status
  setStatus(DE_SUCCESS);

} // end checkDerivatives

} // namespace NonOpt
