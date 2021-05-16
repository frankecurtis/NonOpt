// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <ctime>
#include <vector>

#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptDirectionComputationGradient.hpp"

namespace NonOpt
{

// Add options
void DirectionComputationGradient::addOptions(Options* options,
                                              const Reporter* reporter)
{

  // Add bool options
  options->addBoolOption(reporter,
                         "DCG_fail_on_QP_failure",
                         false,
                         "Determines whether to fail if QP solver ever fails.\n"
                         "Default     : false.");

} // end addOptions

// Set options
void DirectionComputationGradient::setOptions(const Options* options,
                                              const Reporter* reporter)
{

  // Read bool options
  options->valueAsBool(reporter, "DCG_fail_on_QP_failure", fail_on_QP_failure_);

} // end setOptions

// Initialize
void DirectionComputationGradient::initialize(const Options* options,
                                              Quantities* quantities,
                                              const Reporter* reporter) {}

// Iteration header
std::string DirectionComputationGradient::iterationHeader()
{
  return "In. Its.  QP Its. QP   QP KKT    |Step|   |Step|_H";
}

// Iteration null values string
std::string DirectionComputationGradient::iterationNullValues()
{
  return "-------- -------- -- --------- --------- ---------";
}

// Compute direction
void DirectionComputationGradient::computeDirection(const Options* options,
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
  clock_t start_time = clock();

  // try direction computation, terminate on any exception
  try {

    // Declare bool for evaluations
    bool evaluation_success;

    // Check whether to evaluate function with gradient
    if (quantities->evaluateFunctionWithGradient()) {

      // Evaluate current objective
      evaluation_success = quantities->currentIterate()->evaluateObjectiveAndGradient(*quantities);

      // Check for successful evaluation
      if (!evaluation_success) {
        THROW_EXCEPTION(DC_EVALUATION_FAILURE_EXCEPTION, "Direction computation unsuccessful. Evaluation failed.");
      }

    }
    else {

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

    } // end else

    // Declare QP quantities
    std::vector<std::shared_ptr<Vector>> QP_gradient_list;
    std::vector<double> QP_vector;

    // Add pointer to current gradient to list
    QP_gradient_list.push_back(quantities->currentIterate()->gradient());

    // Add linear term value
    QP_vector.push_back(quantities->currentIterate()->objective());

    // Set QP data
    strategies->qpSolver()->setVectorList(QP_gradient_list);
    strategies->qpSolver()->setVector(QP_vector);
    strategies->qpSolver()->setScalar(quantities->trustRegionRadius());
    strategies->qpSolver()->setInexactSolutionTolerance(quantities->stationarityRadius());

    // Solve QP
    strategies->qpSolver()->solveQP(options, reporter, quantities);

    // Convert QP solution to step
    convertQPSolutionToStep(quantities, strategies);

    // Check for QP failure
    if (strategies->qpSolver()->status() != QP_SUCCESS && fail_on_QP_failure_) {
      THROW_EXCEPTION(DC_QP_FAILURE_EXCEPTION, "Direction computation unsuccessful. QP solver failed.");
    }
    else {
      THROW_EXCEPTION(DC_SUCCESS_EXCEPTION, "Direction computation successful.")
    }

  } // end try

  // catch exceptions
  catch (DC_SUCCESS_EXCEPTION& exec) {
    setStatus(DC_SUCCESS);
  } catch (DC_CPU_TIME_LIMIT_EXCEPTION& exec) {
    setStatus(DC_CPU_TIME_LIMIT);
    quantities->incrementDirectionComputationTime(clock() - start_time);
    THROW_EXCEPTION(NONOPT_CPU_TIME_LIMIT_EXCEPTION, "CPU time limit has been reached.");
  } catch (DC_EVALUATION_FAILURE_EXCEPTION& exec) {
    setStatus(DC_EVALUATION_FAILURE);
  } catch (DC_QP_FAILURE_EXCEPTION& exec) {
    setStatus(DC_QP_FAILURE);
  }

  // Print iteration information
  reporter->printf(R_NL, R_PER_ITERATION, " %8d %8d %2d %+.2e %+.2e %+.2e", quantities->innerIterationCounter(), quantities->QPIterationCounter(), strategies->qpSolver()->status(), strategies->qpSolver()->KKTErrorDual(), strategies->qpSolver()->primalSolutionNormInf(), strategies->qpSolver()->dualObjectiveQuadraticValue());

  // Increment total inner iteration counter
  quantities->incrementTotalInnerIterationCounter();

  // Increment total QP iteration counter
  quantities->incrementTotalQPIterationCounter();

  // Increment direction computation time
  quantities->incrementDirectionComputationTime(clock() - start_time);

} // end computeDirection

// Convert QP solution to step
void DirectionComputationGradient::convertQPSolutionToStep(Quantities* quantities,
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
