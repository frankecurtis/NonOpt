// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include <cmath>

#include "NonOptDeclarations.hpp"
#include "NonOptDefinitions.hpp"
#include "NonOptDirectionComputationAggregation.hpp"
#include <iostream>
namespace NonOpt
{

// Add options
void DirectionComputationAggregation::addOptions(Options* options,
                                                 const Reporter* reporter)
{

  // Add bool options
  options->addBoolOption(reporter,
                         "DCAgg_fail_on_iteration_limit",
                         false,
                         "Determines whether to fail if iteration limit exceeded.\n"
                         "Default value: false.");
  options->addBoolOption(reporter,
                         "DCAgg_fail_on_QP_failure",
                         false,
                         "Determines whether to fail if QP solver ever fails.\n"
                         "Default value: false.");
  options->addBoolOption(reporter,
                         "DCAgg_try_shortened_step",
                         true,
                         "Determines whether to consider shortened step if subproblem\n"
                         "solver does not fail after considering full QP step.\n"
                         "Shortened stepsize set by DCAgg_shortened_stepsize parameter.\n"
                         "Default value: true.");

  // Add double options
  options->addDoubleOption(reporter,
                           "DCAgg_downshift_constant",
                           1e-02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Downshifting constant.  The linear term corresponding to an\n"
                           "added cut is set as the objective value at the current iterate\n"
                           "minus this value times the squared norm difference between the\n"
                           "sampled point and the current iterate.\n"
                           "Default value: 1e-02.");
  options->addDoubleOption(reporter,
                           "DCAgg_random_sample_fraction",
                           2.5e-02,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Fraction of number of variables used to determine number of\n"
                           "points to sample in adaptive gradient sampling scheme.\n"
                           "Default value: 1e-01.");
  options->addDoubleOption(reporter,
                           "DCAgg_shortened_stepsize",
                           1e+00,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Shortened stepsize.  If full QP step does not offer desired\n"
                           "objective reduction, then a shortened step corresponding to\n"
                           "this stepsize is considered if DCAgg_try_shortened_step == true.\n"
                           "In particular, the shortened stepsize that is considered is\n"
                           "DCAgg_shortened_stepsize*min(stat. rad.,||qp_step||_inf)/||qp_step||_inf.\n"
                           "Default value: 1e+00.");
  options->addDoubleOption(reporter,
                           "DCAgg_step_acceptance_tolerance",
                           1e-12,
                           0.0,
                           1.0,
                           "Tolerance for step acceptance.\n"
                           "Default value: 1e-12.");
  options->addDoubleOption(reporter,
                           "DCAgg_size_factor",
                           1e+01,
                           0.0,
                           NONOPT_DOUBLE_INFINITY,
                           "Size factor for removing points from point set.  If size of\n"
                           "point set exceeds this factor times the number of variables,\n"
                           "then the oldest members are removed.\n"
                           "Default value: 1e+02.");

  // Add integer options
  options->addIntegerOption(reporter,
                            "DCAgg_inner_iteration_limit",
                            2e+01,
                            0,
                            NONOPT_INT_INFINITY,
                            "Limit on the number of inner iterations that will be performed.\n"
                            "Default value: 2e+01.");

} // end addOptions

// Set options
void DirectionComputationAggregation::setOptions(const Options* options,
                                                 const Reporter* reporter)
{

  // Read bool options
  options->valueAsBool(reporter, "DCAgg_fail_on_iteration_limit", fail_on_iteration_limit_);
  options->valueAsBool(reporter, "DCAgg_fail_on_QP_failure", fail_on_QP_failure_);
  options->valueAsBool(reporter, "DCAgg_try_shortened_step", try_shortened_step_);

  // Read double options
  options->valueAsDouble(reporter, "DCAgg_downshift_constant", downshift_constant_);
  options->valueAsDouble(reporter, "DCAgg_random_sample_fraction", random_sample_fraction_);
  options->valueAsDouble(reporter, "DCAgg_shortened_stepsize", shortened_stepsize_);
  options->valueAsDouble(reporter, "DCAgg_step_acceptance_tolerance", step_acceptance_tolerance_);
  options->valueAsDouble(reporter, "DCAgg_size_factor", full_size_factor_);

  // Read integer options
  options->valueAsInteger(reporter, "DCAgg_inner_iteration_limit", inner_iteration_limit_);

} // end setOptions

// Initialize
void DirectionComputationAggregation::initialize(const Options* options,
                                                 Quantities* quantities,
                                                 const Reporter* reporter) {}

// Iteration header
std::string DirectionComputationAggregation::iterationHeader()
{
  return "In.-Iters.   QP-Iters.  QP   KKT Error   |G. Combo.|     |Step|      |Step|_H ";
}

// Iteration null values string
std::string DirectionComputationAggregation::iterationNullValues()
{
  return "----------  ----------  --  -----------  -----------  -----------  -----------";
}

// Compute direction
void DirectionComputationAggregation::computeDirection(const Options* options,
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
    std::vector<std::shared_ptr<Vector>> QP_gradient_list_agg;
    std::vector<double> QP_vector_agg;

    // Declare new QP data
    std::vector<std::shared_ptr<Vector>> QP_gradient_list;
    std::vector<double> QP_vector;

    // Add pointer to current gradient to list
    QP_gradient_list_agg.push_back(quantities->currentIterate()->gradient());
    QP_gradient_list.push_back(quantities->currentIterate()->gradient());

    // Add linear term value
    QP_vector_agg.push_back(quantities->currentIterate()->objective());
    QP_vector.push_back(quantities->currentIterate()->objective());

    // Loop through point set
    for (int point_count = 0; point_count < (int)quantities->pointSet()->size(); point_count++) {

      // Create difference vector
      std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *(*quantities->pointSet())[point_count]->vector());

      // Check distance between point and current iterate
      if (difference->normInf() <= quantities->stationarityRadius()) {

        // Evaluate gradient
        evaluation_success = (*quantities->pointSet())[point_count]->evaluateGradient(*quantities);

        // Check for successful evaluation
        if (evaluation_success) {

          // Add pointer to gradient in point set to list
          QP_gradient_list_agg.push_back((*quantities->pointSet())[point_count]->gradient());
          QP_gradient_list.push_back((*quantities->pointSet())[point_count]->gradient());

          // Evaluate downshifting value
          double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

          // Add linear term value
          QP_vector_agg.push_back(downshifting_value);
          QP_vector.push_back(downshifting_value);

        } // end if

      } // end if

    } // end for

    // Set QP data
    strategies->qpSolver()->setVectorList(QP_gradient_list_agg);
    strategies->qpSolver()->setVector(QP_vector_agg);
    //    used_full_last_=true;
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

    // Clear data on QP failure
    if (strategies->qpSolver()->status() != QP_SUCCESS) {

      // Clear data
      QP_gradient_list_agg.clear();
      QP_vector_agg.clear();

      QP_gradient_list.clear();
      QP_vector.clear();

      // Add pointer to current gradient to list
      QP_gradient_list_agg.push_back(quantities->currentIterate()->gradient());
      QP_gradient_list.push_back(quantities->currentIterate()->gradient());

      // Add linear term value
      QP_vector_agg.push_back(quantities->currentIterate()->objective());
      QP_vector.push_back(quantities->currentIterate()->objective());

      // Set QP data
      strategies->qpSolver()->setVectorList(QP_gradient_list_agg);
      strategies->qpSolver()->setVector(QP_vector_agg);

      //      used_full_last_=true;

      // Solve QP
      strategies->qpSolver()->solveQP(options, reporter, quantities);

      // Convert QP solution to step
      convertQPSolutionToStep(quantities, strategies);

    } // end if

    do_agg_next_ = true;

    // Inner loop

    while (true) {

      // Flush reporter buffer
      reporter->flushBuffer();

      // Evaluate trial iterate objective
      evaluation_success = quantities->trialIterate()->evaluateObjective(*quantities);

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

      if (do_agg_next_) {

        yk = new double[strategies->qpSolver()->omega_length()];
        gamm = new double[quantities->numberOfVariables()];
        strategies->qpSolver()->dualSolution(yk, gamm);

        std::shared_ptr<Vector> combination_temp(new Vector(quantities->numberOfVariables()));
        double downshifting_obj = 0.0;

        for (int i = 0; i < strategies->qpSolver()->omega_length(); i++) {
          combination_temp->addScaledVector(yk[i], *QP_gradient_list_agg[i]);
          downshifting_obj += yk[i] * (QP_vector[i]);
        }

        delete[] yk;
        delete[] gamm;

        // Clear data
        QP_gradient_list_agg.clear();
        QP_vector_agg.clear();

        // Add pointer to current gradient to list
        QP_gradient_list_agg.push_back(quantities->currentIterate()->gradient());
        QP_vector_agg.push_back(quantities->currentIterate()->objective());

        // Add combination term
        QP_gradient_list_agg.push_back(combination_temp);
        QP_vector_agg.push_back(downshifting_obj);
      }

      // Add trial point, if inside stationarity radius
      if (strategies->qpSolver()->primalSolutionNormInf() <= quantities->stationarityRadius()) {

        // Evaluate trial point gradient
        evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

        // Check for successful evaluation
        if (evaluation_success) {

          // Add trial iterate to point set
          quantities->pointSet()->push_back(quantities->trialIterate());

          // Add pointer to gradient in point set to list
          QP_gradient_list_agg.push_back(quantities->trialIterate()->gradient());
          QP_gradient_list.push_back(quantities->trialIterate()->gradient());

          // Create difference vector
          std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *quantities->trialIterate()->vector());

          // Evaluate downshifting value
          double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

          // Add linear term value
          QP_vector_agg.push_back(downshifting_value);
          QP_vector.push_back(downshifting_value);

          if (!do_agg_next_) {
            QP_gradient_list_new.push_back(quantities->trialIterate()->gradient());
            QP_vector_new.push_back(downshifting_value);
          }

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

        // Evaluate trial gradient
        evaluation_success = quantities->trialIterate()->evaluateGradient(*quantities);

        // Check for objective evaluation success
        if (evaluation_success) {

          // Add trial iterate to point set
          quantities->pointSet()->push_back(quantities->trialIterate());

          // Add pointer to gradient in point set to list
          QP_gradient_list_agg.push_back(quantities->trialIterate()->gradient());
          QP_gradient_list.push_back(quantities->trialIterate()->gradient());

          // Create difference vector
          std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *quantities->trialIterate()->vector());

          // Evaluate downshifting value
          double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

          // Add linear term value
          QP_vector_agg.push_back(downshifting_value);
          QP_vector.push_back(downshifting_value);

          if (!do_agg_next_) {
            QP_gradient_list_new.push_back(quantities->trialIterate()->gradient());
            QP_vector_new.push_back(downshifting_value);
          }

        } // end if

      } // end if (try_shortened_step_)

      // Loop over sample size

      int added_points = 0;
      if (random_sample_fraction_ > 1.0) {
        added_points = (int)random_sample_fraction_;
      }
      else {
        added_points = (int)(random_sample_fraction_ * quantities->numberOfVariables());
      }

      //int added_points=(int) random_sample_fraction_ *quantities->pointSet()->size();
      // Loop over sample size
      for (int point_count = 0; point_count < std::max(1, added_points); point_count++) {
        // Randomly generate new point
        std::shared_ptr<Point> random_point = quantities->currentIterate()->makeNewRandom(quantities->stationarityRadius(), &random_number_generator_);

        // Evaluate gradient at random point
        evaluation_success = random_point->evaluateGradient(*quantities);

        //        std::cout<<point_count<<"\n";
        // Check for gradient evaluation success
        if (evaluation_success) {

          // Add random point to point set
          quantities->pointSet()->push_back(random_point);

          // Add pointer to gradient to list
          QP_gradient_list.push_back(random_point->gradient());
          QP_gradient_list_agg.push_back(random_point->gradient());

          // Create difference vector
          std::shared_ptr<Vector> difference = quantities->currentIterate()->vector()->makeNewLinearCombination(1.0, -1.0, *random_point->vector());

          // Evaluate downshifting value
          double downshifting_value = quantities->currentIterate()->objective() - downshift_constant_ * pow(difference->norm2(), 2.0);

          // Add linear term value (simply current objective value)
          QP_vector.push_back(downshifting_value);
          QP_vector_agg.push_back(downshifting_value);

          if (!do_agg_next_) {
            QP_gradient_list_new.push_back(random_point->gradient());
            QP_vector_new.push_back(downshifting_value);
          }

        } // end if

      } // end for

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

      if (do_agg_next_) {
        if ((double)quantities->pointSet()->size() > std::min(5000.0, full_size_factor_ * (double)quantities->numberOfVariables())) {
          // Add QP data
          strategies->qpSolver()->setVectorList(QP_gradient_list);
          strategies->qpSolver()->setVector(QP_vector);
          do_agg_next_ = false;
        }
        else {
          // Add QP data
          strategies->qpSolver()->setVectorList(QP_gradient_list_agg);
          strategies->qpSolver()->setVector(QP_vector_agg);
        }

        // Solve QP
        strategies->qpSolver()->solveQP(options, reporter, quantities);
      }
      else {
        // Add QP data
        strategies->qpSolver()->addData(QP_gradient_list_new, QP_vector_new);

        // Solve QP hot
        strategies->qpSolver()->solveQPHot(options, reporter, quantities);
      }

      // Convert QP solution to step
      convertQPSolutionToStep(quantities, strategies);

      // Check for termination on QP failure
      if (strategies->qpSolver()->status() != QP_SUCCESS && fail_on_QP_failure_) {
        THROW_EXCEPTION(DC_QP_FAILURE_EXCEPTION, "Direction computation unsuccessful. QP solver failed.");
      }

      // Check for QP failure
      if (strategies->qpSolver()->status() != QP_SUCCESS) {

        // Clear data
        QP_gradient_list_agg.clear();
        QP_vector_agg.clear();

        QP_gradient_list.clear();
        QP_vector.clear();

        // Add pointer to current gradient to list
        QP_gradient_list_agg.push_back(quantities->currentIterate()->gradient());
        QP_gradient_list.push_back(quantities->currentIterate()->gradient());

        // Add linear term value
        QP_vector_agg.push_back(quantities->currentIterate()->objective());
        QP_vector.push_back(quantities->currentIterate()->objective());

        // Set QP data
        strategies->qpSolver()->setVectorList(QP_gradient_list);
        strategies->qpSolver()->setVector(QP_vector);

        //        used_full_last_=true;
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
void DirectionComputationAggregation::convertQPSolutionToStep(Quantities* quantities,
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