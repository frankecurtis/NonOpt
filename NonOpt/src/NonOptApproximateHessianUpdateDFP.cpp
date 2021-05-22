// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptApproximateHessianUpdateDFP.hpp"
#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptVector.hpp"

namespace NonOpt
{

// Add options
void ApproximateHessianUpdateDFP::addOptions(Options* options)
{

  // Add bool options
  options->addBoolOption("DFP_fail_on_tolerance_violation",
                         false,
                         "Indicator for whether to indicate failure on violated tolerance.\n"
                         "Default     : false.");

  // Add double options
  options->addDoubleOption("DFP_correction_threshold_1",
                           1e-20,
                           0.0,
                           1.0,
                           "DFP correction threshold.  If DFP update is corrected, then\n"
                           "              gradient displacement v is set so that inner product with step,\n"
                           "              call it s, is such that <s,v>/<s,s> is at least this threshold.\n"
                           "Default     : 1e-20.");
  options->addDoubleOption("DFP_correction_threshold_2",
                           1e+02,
                           1.0,
                           NONOPT_DOUBLE_INFINITY,
                           "DFP correction threshold.  If DFP update is corrected, then\n"
                           "              gradient displacement v is set so that inner product with step,\n"
                           "              call it s, is such that <v,v>/<s,v> is at most this threshold.\n"
                           "Default     : 1e+02.");
  options->addDoubleOption("DFP_norm_tolerance",
                           1e-20,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for allowing a DFP update to occur.  If the norm\n"
                           "              of the iterate displacement or the gradient displacement is\n"
                           "              greater than this tolerance, then the DFP update may occur;\n"
                           "              else, it is skipped.\n"
                           "Default     : 1e-20.");
  options->addDoubleOption("DFP_product_tolerance",
                           1e-20,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Tolerance for allowing a DFP update to occur.  If the inner\n"
                           "              product between the iterate and gradient displacements is at\n"
                           "              least this tolerance times the product of the 2-norms of the\n"
                           "              displacements, then a DFP update occurs; else, it is skipped.\n"
                           "Default     : 1e-20.");

} // end addOptions

// Set options
void ApproximateHessianUpdateDFP::setOptions(Options* options)
{

  // Read bool options
  options->valueAsBool("DFP_fail_on_tolerance_violation", fail_on_tolerance_violation_);

  // Read double options
  options->valueAsDouble("DFP_correction_threshold_1", correction_threshold_1_);
  options->valueAsDouble("DFP_correction_threshold_2", correction_threshold_2_);
  options->valueAsDouble("DFP_norm_tolerance", norm_tolerance_);
  options->valueAsDouble("DFP_product_tolerance", product_tolerance_);

} // end setOptions

// Initialize
void ApproximateHessianUpdateDFP::initialize(const Options* options,
                                             Quantities* quantities,
                                             const Reporter* reporter)
{

  // Indicate that initial update has not yet been performed
  initial_update_performed_ = false;

} // end initialize

// Update approximate Hessian
void ApproximateHessianUpdateDFP::updateApproximateHessian(const Options* options,
                                                           Quantities* quantities,
                                                           const Reporter* reporter,
                                                           Strategies* strategies)
{

  // Initialize values
  setStatus(AH_UNSET);
  double correction_scalar = 0.0;
  bool perform_update = false;

  // try update, terminate on any exception
  try {

    // Declare iterate displacement
    Vector iterate_displacement(quantities->numberOfVariables());

    // Set iterate displacement
    iterate_displacement.linearCombination(1.0,
                                           *quantities->trialIterate()->vector(),
                                           -1.0,
                                           *quantities->currentIterate()->vector());

    // Check for iterate displacement norm tolerance violation
    if (iterate_displacement.norm2() <= norm_tolerance_) {
      THROW_EXCEPTION(AH_NORM_TOLERANCE_VIOLATION_EXCEPTION, "Approximate Hessian update unsuccessful. Norm tolerance violated.");
    }

    // Declare gradient displacement
    Vector gradient_displacement(quantities->numberOfVariables());

    // Evaluate current iterate gradient
    bool evaluation_success = quantities->currentIterate()->evaluateGradient(*quantities);

    // Check for successful evaluation
    if (!evaluation_success) {
      THROW_EXCEPTION(AH_EVALUATION_FAILURE_EXCEPTION, "Approximate Hessian update unsuccessful. Evaluation failed.");
    }

    // Evaluate trial iterate gradient
    evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

    // Check for successful evaluation
    if (!evaluation_success) {
      THROW_EXCEPTION(AH_EVALUATION_FAILURE_EXCEPTION, "Approximate Hessian update unsuccessful. Evaluation failed.");
    }

    // Set gradient displacement
    gradient_displacement.linearCombination(1.0,
                                            *quantities->trialIterate()->gradient(),
                                            -1.0,
                                            *quantities->currentIterate()->gradient());

    // Check for gradient displacement norm tolerance violation
    if (gradient_displacement.norm2() <= norm_tolerance_) {
      THROW_EXCEPTION(AH_NORM_TOLERANCE_VIOLATION_EXCEPTION, "Approximate Hessian update unsuccessful. Norm tolerance violated.");
    }

    // Evaluate correction scalar
    evaluateSelfCorrectingScalar(iterate_displacement, gradient_displacement, correction_scalar);

    // Update gradient displacement
    if (correction_scalar > 0.0) {

      // Scale it
      gradient_displacement.scale(1.0 - correction_scalar);

      // Add scaled vector
      gradient_displacement.addScaledVector(correction_scalar, iterate_displacement);

    } // end if

    // Determine whether approximate Hessian should be updated
    perform_update = (iterate_displacement.innerProduct(gradient_displacement) >= product_tolerance_ * iterate_displacement.norm2() * gradient_displacement.norm2());

    // Check for inner product violation
    if (!perform_update) {
      THROW_EXCEPTION(AH_PRODUCT_TOLERANCE_VIOLATION_EXCEPTION, "Approximate Hessian update unsuccessful. Product tolerance violated.");
    }

    // Check for initial scaling
    if (!initial_update_performed_ && quantities->approximateHessianInitialScaling()) {
      strategies->symmetricMatrix()->setAsDiagonal(quantities->numberOfVariables(), iterate_displacement.innerProduct(gradient_displacement) / pow(gradient_displacement.norm2(), 2.0));
      initial_update_performed_ = true;
    } // end if

    // Update approximate Hessian
    strategies->symmetricMatrix()->update(iterate_displacement, gradient_displacement);

    // Terminate
    THROW_EXCEPTION(AH_SUCCESS_EXCEPTION, "Approximate Hessian update successful.");

  } // end try

  // catch exceptions
  catch (AH_SUCCESS_EXCEPTION& exec) {
    setStatus(AH_SUCCESS);
  } catch (AH_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(AH_EVALUATION_FAILURE);
  } catch (AH_NORM_TOLERANCE_VIOLATION_EXCEPTION& exec) {
    if (fail_on_tolerance_violation_) {
      setStatus(AH_NORM_TOLERANCE_VIOLATION);
    }
    else {
      setStatus(AH_SUCCESS);
    }
  } // end catch
  catch (AH_PRODUCT_TOLERANCE_VIOLATION_EXCEPTION& exec) {
    if (fail_on_tolerance_violation_) {
      setStatus(AH_PRODUCT_TOLERANCE_VIOLATION);
    }
    else {
      setStatus(AH_SUCCESS);
    }
  } // end catch

  // Print messages
  reporter->printf(R_NL, R_PER_ITERATION, " %+.2e %2d", correction_scalar, perform_update);

} // end updateApproximateHessian

// Evaluate self-correcting DFP scalar
void ApproximateHessianUpdateDFP::evaluateSelfCorrectingScalar(Vector& s,
                                                               Vector& y,
                                                               double& scalar)
{

  // Compute products
  double u = pow(s.norm2(), 2.0);
  double v = s.innerProduct(y);
  double w = pow(y.norm2(), 2.0);

  // Set scalar for lower bound
  double scalar1 = 0.0;
  if (u <= 0.0) {
    scalar1 = 1.0;
  }
  else {
    if (v / u < correction_threshold_1_) {
      if (correction_threshold_1_ * u - v > 0.0 &&
          u - v > 0.0) {
        scalar1 = (correction_threshold_1_ * u - v) / (u - v);
      }
      else {
        scalar1 = 0.0;
      }
    }
    if (scalar1 > 0.0 &&
        pow(scalar1, 2.0) * u + scalar1 * (1.0 - scalar1) * v + pow(1.0 - scalar1, 2.0) * w > correction_threshold_1_ * (scalar1 * u + (1 - scalar1) * v)) {
      scalar1 = 1.0;
    }
  } // end else

  // Set scalar for upper bound
  double scalar2 = 0.0;
  if (v <= 0.0) {
    scalar2 = 1.0;
  }
  else if (w / v > correction_threshold_2_) {
    double temporary1 = u - 2.0 * v + w;
    double temporary2 = 2.0 * (v - w) - correction_threshold_2_ * (u - v);
    double temporary3 = w - correction_threshold_2_ * v;
    if (temporary3 > 0.0 &&
        -temporary2 + sqrt(pow(temporary2, 2.0) - 4.0 * temporary1 * temporary3) > 0.0) {
      scalar2 = 2.0 * temporary3 / (-temporary2 + sqrt(pow(temporary2, 2.0) - 4.0 * temporary1 * temporary3));
    }
    else {
      scalar2 = 1.0;
    }
  } // end else if

  // Set scalar
  scalar = fmax(scalar1, scalar2);

} // end evaluateSelfCorrectingScalar

} // namespace NonOpt
