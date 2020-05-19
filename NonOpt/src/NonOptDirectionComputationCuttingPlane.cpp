// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptDirectionComputationCuttingPlane.hpp"

namespace NonOpt
{

// Add options
void DirectionComputationCuttingPlane::addOptions(Options* options,
                                                  const Reporter* reporter)
{

  // Add bool options
  options->addBoolOption(reporter,
                         "DCCP_add_far_points",
                         false,
                         "Determines whether to add points far outside stationarity\n"
                         "radius to point set during subproblem solve.\n"
                         "Default value: false.");
  options->addBoolOption(reporter,
                         "DCCP_fail_on_iteration_limit",
                         false,
                         "Determines whether to fail if iteration limit exceeded.\n"
                         "Default value: false.");
  options->addBoolOption(reporter,
                         "DCCP_fail_on_QP_failure",
                         false,
                         "Determines whether to fail if QP solver ever fails.\n"
                         "Default value: false.");
  options->addBoolOption(reporter,
                         "DCCP_try_shortened_step",
                         true,
                         "Determines whether to consider shortened step if subproblem\n"
                         "solver does not terminate after considering full QP step.\n"
                         "Shortened stepsize set by DCCP_shortened_stepsize parameter.\n"
                         "Default value: true.");

  // Add double options
  options->addDoubleOption(reporter,
                           "DCCP_downshift_constant",
                           1e-02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Downshifting constant.  The linear term corresponding to an\n"
                           "added cut is set as the minimum of the linearization value\n"
                           "corresponding to the bundle point and the negative of this value\n"
                           "times the squared norm difference between the bundle point and\n"
                           "the current iterate.\n"
                           "Default value: 1e-02.");
  options->addDoubleOption(reporter,
                           "DCCP_shortened_stepsize",
                           1e+00,
                           0.0,
                           1.0,
                           "Shortened stepsize.  If full QP step does not offer desired\n"
                           "objective reduction, then a shortened step corresponding to\n"
                           "this stepsize is considered if DCCP_try_shortened_step == true.\n"
                           "In particular, the shortened stepsize that is considered is\n"
                           "DCCP_shortened_stepsize*min(stat. rad.,||qp_step||_inf)/||qp_step||_inf.\n"
                           "Default value: 1e+00.");
  options->addDoubleOption(reporter,
                           "DCCP_step_acceptance_tolerance",
                           1e-12,
                           0.0,
                           1.0,
                           "Tolerance for step acceptance.\n"
                           "Default value: 1e-12.");

  // Add integer options
  options->addIntegerOption(reporter,
                            "DCCP_inner_iteration_limit",
                            2e+01,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of inner iterations that will be performed.\n"
                            "Default value: 2e+01.");

} // end addOptions

// Set options
void DirectionComputationCuttingPlane::setOptions(const Options* options,
                                                  const Reporter* reporter)
{

  // Read integer options
  options->valueAsBool(reporter, "DCCP_add_far_points", add_far_points_);
  options->valueAsBool(reporter, "DCCP_fail_on_iteration_limit", fail_on_iteration_limit_);
  options->valueAsBool(reporter, "DCCP_fail_on_QP_failure", fail_on_QP_failure_);
  options->valueAsBool(reporter, "DCCP_try_shortened_step", try_shortened_step_);

  // Read double options
  options->valueAsDouble(reporter, "DCCP_downshift_constant", downshift_constant_);
  options->valueAsDouble(reporter, "DCCP_shortened_stepsize", shortened_stepsize_);
  options->valueAsDouble(reporter, "DCCP_step_acceptance_tolerance", step_acceptance_tolerance_);

  // Read integer options
  options->valueAsInteger(reporter, "DCCP_inner_iteration_limit", inner_iteration_limit_);

} // end setOptions

// Initialize
void DirectionComputationCuttingPlane::initialize(const Options* options,
                                                  Quantities* quantities,
                                                  const Reporter* reporter) {}

// Iteration header
std::string DirectionComputationCuttingPlane::iterationHeader()
{
  return "In.-Iters.   QP-Iters.  QP   KKT Error   |G. Combo.|     |Step|      |Step|_H ";
}

// Iteration null values string
std::string DirectionComputationCuttingPlane::iterationNullValues()
{
  return "----------  ----------  --  -----------  -----------  -----------  -----------";
}

// Compute direction
void DirectionComputationCuttingPlane::computeDirection(const Options* options,
                                                        Quantities* quantities,
                                                        const Reporter* reporter,
                                                        Strategies* strategies)
{

  // Initialize values
  setStatus(DC_UNSET);
  strategies->qpSolver()->setPrimalSolutionToZero();
  quantities->resetInnerIterationCounter();
  quantities->resetQPIterationCounter();
  quantities->setTrialIterateToCurrentIterate();

  // try direction computation, terminate on any exception
  try {

    // Initialize boolean for evaluation
    bool evaluation_success = false;

    // Evaluate current objective
    evaluation_success = quantities->currentIterate()->evaluateObjective(*quantities);

    // Check for successful evaluation
    if (!evaluation_success) {
      THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION, "Direction computation unsuccessful. Evaluation failed.");
    }

    // Evaluate current gradient
    evaluation_success = quantities->currentIterate()->evaluateGradient(*quantities);

    // Check for successful evaluation
    if (!evaluation_success) {
      THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION, "Direction computation unsuccessful. Evaluation failed.");
    }

    // Declare QP quantities
    std::vector<std::shared_ptr<Vector>> QP_gradient_list;
    std::vector<double> QP_vector;

    // Add pointer to current gradient to list
    QP_gradient_list.push_back(quantities->currentIterate()->gradient());

    // Add linear term value
    QP_vector.push_back(quantities->currentIterate()->objective());

    // Loop through point set
    for (int point_count = 0; point_count < (int)quantities->pointSet()->size(); point_count++) {

      // Create difference vector
      std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *(*quantities->pointSet())[point_count]->vector());

      // Check distance between point and current iterate
      if (difference->normInf() <= quantities->stationarityRadius()) {

        // Evaluate gradient
        evaluation_success = (*quantities->pointSet())[point_count]->evaluateGradient(*quantities);

        // Check for evaluation success
        if (evaluation_success) {

          // Add pointer to gradient in point set to list
          QP_gradient_list.push_back((*quantities->pointSet())[point_count]->gradient());

          // Evaluate linearization and downshifting values
          double linearization_value = (*quantities->pointSet())[point_count]->objective() + (*quantities->pointSet())[point_count]->gradient()->innerProduct(*quantities->currentIterate()->vector()) - (*quantities->pointSet())[point_count]->gradient()->innerProduct(*((*quantities->pointSet())[point_count])->vector());
          double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

          // Add linear term value based on cutting plane
          QP_vector.push_back(fmin(linearization_value, downshifting_value));

        } // end if

      } // end if

    } // end for

    // Set QP data
    strategies->qpSolver()->setVectorList(QP_gradient_list);
    strategies->qpSolver()->setVector(QP_vector);
    strategies->qpSolver()->setScalar(quantities->trustRegionRadius());
    strategies->qpSolver()->setInexactSolutionTolerance(quantities->stationarityRadius());

    // Solve QP
    strategies->qpSolver()->solveQP(options, reporter, quantities);

    // Convert QP solution to step
    convertQPSolutionToStep(quantities, strategies);

    // Check for termination on QP failure
    if (strategies->qpSolver()->status() != QP_SUCCESS && fail_on_QP_failure_) {
      THROW_EXCEPTION(DC_QP_FAILURE_EXCEPTION, "Direction computation unsuccessful. QP solver failed.");
    }

    // Check for QP failure
    if (strategies->qpSolver()->status() != QP_SUCCESS) {

      // Clear data
      QP_gradient_list.clear();
      QP_vector.clear();

      // Add pointer to current gradient to list
      QP_gradient_list.push_back(quantities->currentIterate()->gradient());

      // Add linear term value
      QP_vector.push_back(quantities->currentIterate()->objective());

      // Set QP data
      strategies->qpSolver()->setVectorList(QP_gradient_list);
      strategies->qpSolver()->setVector(QP_vector);

      // Solve QP
      strategies->qpSolver()->solveQP(options, reporter, quantities);

      // Convert QP solution to step
      convertQPSolutionToStep(quantities, strategies);

    } // end if

    // Inner loop
    while (true) {

      // Flush buffer
      reporter->flushBuffer();

      // Evaluate trial iterate objective
      evaluation_success = quantities->trialIterate()->evaluateObjective(*quantities);

      // Check for sufficient decrease
      if (evaluation_success &&
          (quantities->trialIterate()->objective() - quantities->currentIterate()->objective() < -step_acceptance_tolerance_ * std::max(strategies->qpSolver()->combinationTranslatedNorm2Square(), strategies->qpSolver()->primalSolutionNorm2Square()) ||
           (strategies->qpSolver()->primalSolutionNormInf() <= quantities->stationarityRadius() &&
            strategies->qpSolver()->combinationNormInf() <= quantities->stationarityRadius() &&
            strategies->qpSolver()->combinationTranslatedNormInf() <= quantities->stationarityRadius()))) {
        THROW_EXCEPTION(DC_SUCCESS_EXCEPTION, "Direction computation successful.");
      }

      // Check for inner iteration limit
      if (quantities->innerIterationCounter() > inner_iteration_limit_) {
        if (fail_on_iteration_limit_) {
          THROW_EXCEPTION(DC_ITERATION_LIMIT_EXCEPTION, "Direction computation unsuccessful. Iterate limit exceeded.");
        }
        else {
          THROW_EXCEPTION(DC_SUCCESS_EXCEPTION, "Direction computation successful.");
        }
      } // end if

      // Declare new QP data
      std::vector<std::shared_ptr<Vector>> QP_gradient_list_new;
      std::vector<double> QP_vector_new;

      // Check if adding far points
      if (add_far_points_) {

        // Check for objective evaluation success
        if (evaluation_success) {

          // Evaluate trial iterate gradient
          evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

          // Check for gradient evaluation success
          if (evaluation_success) {

            // Add trial iterate to point set
            quantities->pointSet()->push_back(quantities->trialIterate());

            // Add pointer to gradient in point set to list
            QP_gradient_list_new.push_back(quantities->trialIterate()->gradient());

            // Create difference vector
            std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *quantities->trialIterate()->vector());

            // Evaluate linearization and downshifting values
            double linearization_value = quantities->trialIterate()->objective() + quantities->trialIterate()->gradient()->innerProduct(*quantities->currentIterate()->vector()) - quantities->trialIterate()->gradient()->innerProduct(*quantities->trialIterate()->vector());
            double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

            // Add linear term value based on cutting plane
            QP_vector_new.push_back(fmin(linearization_value, downshifting_value));

          } // end if

        } // end if

      } // end if

      // Try shortened step?
      if (try_shortened_step_) {

        // Set shortened stepsize
        double shortened_stepsize = shortened_stepsize_ * fmin(quantities->stationarityRadius(), strategies->qpSolver()->primalSolutionNormInf()) / strategies->qpSolver()->primalSolutionNormInf();

        // Compute shortened trial iterate
        quantities->setTrialIterate(quantities->currentIterate()->makeNewLinearCombination(1.0, shortened_stepsize, *quantities->direction()));

        // Evaluate trial objective
        evaluation_success = quantities->trialIterate()->evaluateObjective(*quantities);

        // Check for sufficient decrease
        if (evaluation_success &&
            (quantities->trialIterate()->objective() - quantities->currentIterate()->objective() < -step_acceptance_tolerance_ * shortened_stepsize * std::max(strategies->qpSolver()->combinationTranslatedNorm2Square(), strategies->qpSolver()->primalSolutionNorm2Square()) ||
             (strategies->qpSolver()->primalSolutionNormInf() <= quantities->stationarityRadius() &&
              strategies->qpSolver()->combinationNormInf() <= quantities->stationarityRadius() &&
              strategies->qpSolver()->combinationTranslatedNormInf() <= quantities->stationarityRadius()))) {
          THROW_EXCEPTION(DC_SUCCESS_EXCEPTION, "Direction computation successful.");
        }

        // Check for objective evaluation success
        if (evaluation_success) {

          // Evaluate trial gradient
          evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

          // Check for gradient evaluation success
          if (evaluation_success) {

            // Add trial iterate to point set
            quantities->pointSet()->push_back(quantities->trialIterate());

            // Add pointer to gradient in point set to list
            QP_gradient_list_new.push_back(quantities->trialIterate()->gradient());

            // Create difference vector
            std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *quantities->trialIterate()->vector());

            // Evaluate linearization and downshifting values
            double linearization_value = quantities->trialIterate()->objective() + quantities->trialIterate()->gradient()->innerProduct(*quantities->currentIterate()->vector()) - quantities->trialIterate()->gradient()->innerProduct(*quantities->trialIterate()->vector());
            double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

            // Add linear term value based on cutting plane
            QP_vector_new.push_back(fmin(linearization_value, downshifting_value));

          } // end if

        } // end if

      } // end if (try_shortened_step_)

      // Print QP solve / step information
      reporter->printf(R_NL, R_PER_INNER_ITERATION, "  %10d  %10d  %2d  %+.4e  %+.4e  %+.4e  %+.4e", quantities->innerIterationCounter(), quantities->QPIterationCounter(), strategies->qpSolver()->status(), strategies->qpSolver()->KKTErrorDual(), strategies->qpSolver()->combinationNormInf(), strategies->qpSolver()->primalSolutionNormInf(), strategies->qpSolver()->dualObjectiveQuadraticValue());

      // Set blank solve string
      std::string blank_solve = "";
      if (strategies->lineSearch()->iterationNullValues().length() > 0) {
        blank_solve += "  ";
        blank_solve += strategies->lineSearch()->iterationNullValues();
      }
      if (strategies->inverseHessianUpdate()->iterationNullValues().length() > 0) {
        blank_solve += "  ";
        blank_solve += strategies->inverseHessianUpdate()->iterationNullValues();
      }
      if (strategies->pointSetUpdate()->iterationNullValues().length() > 0) {
        blank_solve += "  ";
        blank_solve += strategies->pointSetUpdate()->iterationNullValues();
      }

      // Print blank solve information
      reporter->printf(R_NL, R_PER_INNER_ITERATION, "%s\n%s", blank_solve.c_str(), quantities->iterationNullValues().c_str());

      // Add QP data
      strategies->qpSolver()->addData(QP_gradient_list_new, QP_vector_new);

      // Solve QP hot
      strategies->qpSolver()->solveQPHot(options, reporter, quantities);

      // Convert QP solution to step
      convertQPSolutionToStep(quantities, strategies);

      // Check for termination on QP failure
      if (strategies->qpSolver()->status() != QP_SUCCESS && fail_on_QP_failure_) {
        THROW_EXCEPTION(DC_QP_FAILURE_EXCEPTION, "Direction computation unsuccessful. QP solver failed.");
      }

      // Check for QP failure
      if (strategies->qpSolver()->status() != QP_SUCCESS) {

        // Clear data
        QP_gradient_list.clear();
        QP_vector.clear();

        // Add pointer to current gradient to list
        QP_gradient_list.push_back(quantities->currentIterate()->gradient());

        // Add linear term value
        QP_vector.push_back(quantities->currentIterate()->objective());

        // Set QP data
        strategies->qpSolver()->setVectorList(QP_gradient_list);
        strategies->qpSolver()->setVector(QP_vector);

        // Solve QP
        strategies->qpSolver()->solveQP(options, reporter, quantities);

        // Convert QP solution to step
        convertQPSolutionToStep(quantities, strategies);

      } // end if

    } // end while

  } // end try

  // catch exceptions
  catch (DC_SUCCESS_EXCEPTION& exec) {
    setStatus(DC_SUCCESS);
  } catch (DC_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(DC_EVALUATION_FAILURE);
  } catch (DC_ITERATION_LIMIT_EXCEPTION& exec) {
    setStatus(DC_ITERATION_LIMIT);
  } catch (DC_QP_FAILURE_EXCEPTION& exec) {
    setStatus(DC_QP_FAILURE);
  }

  // Print iteration information
  reporter->printf(R_NL, R_PER_ITERATION, "  %10d  %10d  %2d  %+.4e  %+.4e  %+.4e  %+.4e", quantities->innerIterationCounter(), quantities->QPIterationCounter(), strategies->qpSolver()->status(), strategies->qpSolver()->KKTErrorDual(), strategies->qpSolver()->combinationNormInf(), strategies->qpSolver()->primalSolutionNormInf(), strategies->qpSolver()->dualObjectiveQuadraticValue());

  // Increment total inner iteration counter
  quantities->incrementTotalInnerIterationCounter();

  // Increment total QP iteration counter
  quantities->incrementTotalQPIterationCounter();

} // end computeDirection

// Convert QP solution to step
void DirectionComputationCuttingPlane::convertQPSolutionToStep(Quantities* quantities,
                                                               Strategies* strategies)
{

  // Increment QP iteration counter
  quantities->incrementQPIterationCounter(strategies->qpSolver()->numberOfIterations());

  // Increment inner iteration counter
  quantities->incrementInnerIterationCounter(1);

  // Get primal solution
  strategies->qpSolver()->primalSolution(quantities->direction()->valuesModifiable());

  // Set trial iterate
  quantities->setTrialIterate(quantities->currentIterate()->makeNewLinearCombination(1.0, 1.0, *quantities->direction()));

} // end convertQPSolutionToStep

} // namespace NonOpt
