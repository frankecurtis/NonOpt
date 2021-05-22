// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptStrategies.hpp"
#include "NonOptApproximateHessianUpdateBFGS.hpp"
#include "NonOptApproximateHessianUpdateDFP.hpp"
#include "NonOptDerivativeCheckerFiniteDifference.hpp"
#include "NonOptDirectionComputationCuttingPlane.hpp"
#include "NonOptDirectionComputationGradient.hpp"
#include "NonOptDirectionComputationGradientCombination.hpp"
#include "NonOptLineSearchBacktracking.hpp"
#include "NonOptLineSearchWeakWolfe.hpp"
#include "NonOptPointSetUpdateProximity.hpp"
#include "NonOptQPSolverDualActiveSet.hpp"
#include "NonOptSymmetricMatrixDense.hpp"
#include "NonOptSymmetricMatrixLimitedMemory.hpp"
#include "NonOptTerminationBasic.hpp"
#include "NonOptTerminationSecondQP.hpp"

namespace NonOpt
{

// Add options
void Strategies::addOptions(Options* options)
{

  // Add string options
  options->addStringOption("approximate_hessian_update",
                           "BFGS",
                           "Approximate Hessian update strategy to use.\n"
                           "Default     : BFGS.");
  options->addStringOption("derivative_checker",
                           "FiniteDifference",
                           "Derivative checker strategy to use.\n"
                           "Default     : FiniteDifference.");
  options->addStringOption("direction_computation",
                           "CuttingPlane",
                           "Direction computation strategy to use.\n"
                           "Default     : CuttingPlane.");
  options->addStringOption("line_search",
                           "WeakWolfe",
                           "Line search strategy to use.\n"
                           "Default     : WeakWolfe.");
  options->addStringOption("point_set_update",
                           "Proximity",
                           "Point set update strategy to use.\n"
                           "Default     : Proximity.");
  options->addStringOption("qp_solver",
                           "DualActiveSet",
                           "QP solver strategy to use.\n"
                           "Default     : DualActiveSet.");
  options->addStringOption("symmetric_matrix",
                           "Dense",
                           "Symmetric matrix strategy to use.\n"
                           "Default     : Dense.");
  options->addStringOption("termination",
                           "Basic",
                           "Termination strategy to use.\n"
                           "Default     : Basic.");

  // Add options for approximate Hessian update strategies
  std::shared_ptr<ApproximateHessianUpdate> approximate_hessian_update;
  approximate_hessian_update = std::make_shared<ApproximateHessianUpdateBFGS>();
  approximate_hessian_update->addOptions(options);
  approximate_hessian_update = std::make_shared<ApproximateHessianUpdateDFP>();
  approximate_hessian_update->addOptions(options);
  // ADD NEW APPROXIMATE HESSIAN UPDATE STRATEGIES HERE //

  // Add options for derivative checker strategies
  std::shared_ptr<DerivativeChecker> derivative_checker;
  derivative_checker = std::make_shared<DerivativeCheckerFiniteDifference>();
  derivative_checker->addOptions(options);
  // ADD NEW POINT SET UPDATE STRATEGIES HERE //

  // Add options for direction computation strategies
  std::shared_ptr<DirectionComputation> direction_computation;
  direction_computation = std::make_shared<DirectionComputationCuttingPlane>();
  direction_computation->addOptions(options);
  direction_computation = std::make_shared<DirectionComputationGradient>();
  direction_computation->addOptions(options);
  direction_computation = std::make_shared<DirectionComputationGradientCombination>();
  direction_computation->addOptions(options);
  // ADD NEW DIRECTION COMPUTATION STRATEGIES HERE //

  // Add options for line search strategies
  std::shared_ptr<LineSearch> line_search;
  line_search = std::make_shared<LineSearchBacktracking>();
  line_search->addOptions(options);
  line_search = std::make_shared<LineSearchWeakWolfe>();
  line_search->addOptions(options);
  // ADD NEW LINE SEARCH STRATEGIES HERE //

  // Add options for point set update strategies
  std::shared_ptr<PointSetUpdate> point_set_update;
  point_set_update = std::make_shared<PointSetUpdateProximity>();
  point_set_update->addOptions(options);
  // ADD NEW POINT SET UPDATE STRATEGIES HERE //

  // Add options for QP solver strategies
  std::shared_ptr<QPSolver> qp_solver;
  qp_solver = std::make_shared<QPSolverDualActiveSet>();
  qp_solver->addOptions(options);
  // ADD NEW QP SOLVER STRATEGIES HERE //

  // Add options for symmetric matrix strategies
  std::shared_ptr<SymmetricMatrix> symmetric_matrix;
  symmetric_matrix = std::make_shared<SymmetricMatrixDense>();
  symmetric_matrix->addOptions(options);
  symmetric_matrix = std::make_shared<SymmetricMatrixLimitedMemory>();
  symmetric_matrix->addOptions(options);
  // ADD NEW SYMMETRIC MATRIX STRATEGIES HERE //

  // Add options for termination strategies
  std::shared_ptr<Termination> termination;
  termination = std::make_shared<TerminationBasic>();
  termination->addOptions(options);
  termination = std::make_shared<TerminationSecondQP>();
  termination->addOptions(options);
  // ADD NEW TERMINATION STRATEGIES HERE //

} // end addOptions

// Set options
void Strategies::setOptions(Options* options)
{

  // Declare strategy names
  std::string approximate_hessian_update_name;
  std::string derivative_checker_name;
  std::string direction_computation_name;
  std::string line_search_name;
  std::string point_set_update_name;
  std::string qp_solver_name;
  std::string symmetric_matrix_name;
  std::string termination_name;

  // Read integer options
  options->valueAsString("approximate_hessian_update", approximate_hessian_update_name);
  options->valueAsString("derivative_checker", derivative_checker_name);
  options->valueAsString("direction_computation", direction_computation_name);
  options->valueAsString("line_search", line_search_name);
  options->valueAsString("point_set_update", point_set_update_name);
  options->valueAsString("qp_solver", qp_solver_name);
  options->valueAsString("symmetric_matrix", symmetric_matrix_name);
  options->valueAsString("termination", termination_name);

  // Set approximate Hessian update strategy
  if (approximate_hessian_update_name.compare("BFGS") == 0) {
    approximate_hessian_update_ = std::make_shared<ApproximateHessianUpdateBFGS>();
  }
  else if (approximate_hessian_update_name.compare("DFP") == 0) {
    approximate_hessian_update_ = std::make_shared<ApproximateHessianUpdateDFP>();
  }
  else {
    approximate_hessian_update_ = std::make_shared<ApproximateHessianUpdateBFGS>();
  }

  // Set derivative checker strategy
  if (derivative_checker_name.compare("FiniteDifference") == 0) {
    derivative_checker_ = std::make_shared<DerivativeCheckerFiniteDifference>();
  }
  else {
    derivative_checker_ = std::make_shared<DerivativeCheckerFiniteDifference>();
  }

  // Set direction computation strategy
  if (direction_computation_name.compare("CuttingPlane") == 0) {
    direction_computation_ = std::make_shared<DirectionComputationCuttingPlane>();
  }
  else if (direction_computation_name.compare("Gradient") == 0) {
    direction_computation_ = std::make_shared<DirectionComputationGradient>();
  }
  else if (direction_computation_name.compare("GradientCombination") == 0) {
    direction_computation_ = std::make_shared<DirectionComputationGradientCombination>();
  }
  else {
    direction_computation_ = std::make_shared<DirectionComputationCuttingPlane>();
  }

  // Set line search strategy
  if (line_search_name.compare("Backtracking") == 0) {
    line_search_ = std::make_shared<LineSearchBacktracking>();
  }
  else if (line_search_name.compare("WeakWolfe") == 0) {
    line_search_ = std::make_shared<LineSearchWeakWolfe>();
  }
  else {
    line_search_ = std::make_shared<LineSearchWeakWolfe>();
  }

  // Set point set update strategy
  if (point_set_update_name.compare("Proximity") == 0) {
    point_set_update_ = std::make_shared<PointSetUpdateProximity>();
  }
  else {
    point_set_update_ = std::make_shared<PointSetUpdateProximity>();
  }

  // Set QP solver strategy
  if (qp_solver_name.compare("DualActiveSet") == 0) {
    qp_solver_ = std::make_shared<QPSolverDualActiveSet>();
    qp_solver_termination_ = std::make_shared<QPSolverDualActiveSet>();
  }
  else {
    qp_solver_ = std::make_shared<QPSolverDualActiveSet>();
    qp_solver_termination_ = std::make_shared<QPSolverDualActiveSet>();
  }

  // Set symmetric matrix strategy
  if (symmetric_matrix_name.compare("Dense") == 0) {
    symmetric_matrix_ = std::make_shared<SymmetricMatrixDense>();
  }
  else if (symmetric_matrix_name.compare("LimitedMemory") == 0) {
    symmetric_matrix_ = std::make_shared<SymmetricMatrixLimitedMemory>();
  }
  else {
    symmetric_matrix_ = std::make_shared<SymmetricMatrixDense>();
  }
  // ALWAYS USE LIMITED MEMORY FOR TERMINATION QP
  symmetric_matrix_termination_ = std::make_shared<SymmetricMatrixLimitedMemory>();

  // Set termination strategy
  if (termination_name.compare("Basic") == 0) {
    termination_ = std::make_shared<TerminationBasic>();
  }
  else if (termination_name.compare("SecondQP") == 0) {
    termination_ = std::make_shared<TerminationSecondQP>();
  }
  else {
    termination_ = std::make_shared<TerminationBasic>();
  }

  // Set approximate Hessian update options
  approximate_hessian_update_->setOptions(options);

  // Set derivative checker options
  derivative_checker_->setOptions(options);

  // Set direction computation options
  direction_computation_->setOptions(options);

  // Set line search options
  line_search_->setOptions(options);

  // Set point set update options
  point_set_update_->setOptions(options);

  // Set QP solver options
  qp_solver_->setOptions(options);
  qp_solver_termination_->setOptions(options);

  // Set symmetric matrix options
  symmetric_matrix_->setOptions(options);
  symmetric_matrix_termination_->setOptions(options);

  // Set termination options
  termination_->setOptions(options);

} // end setOptions

// Initialize
void Strategies::initialize(const Options* options,
                            Quantities* quantities,
                            const Reporter* reporter)
{

  // Initialize approximate Hessian update
  approximate_hessian_update_->initialize(options, quantities, reporter);

  // Initialize derivative checker
  derivative_checker_->initialize(options, quantities, reporter);

  // Initialize direction computation
  direction_computation_->initialize(options, quantities, reporter);

  // Initialize line search
  line_search_->initialize(options, quantities, reporter);

  // Initialize point set update
  point_set_update_->initialize(options, quantities, reporter);

  // Initialize QP data
  qp_solver_->initialize(options, quantities, reporter);
  qp_solver_termination_->initialize(options, quantities, reporter);

  // Initialize symmetric matrix
  symmetric_matrix_->initialize(options, quantities, reporter);
  symmetric_matrix_termination_->initialize(options, quantities, reporter);

  // Initialize termination
  termination_->initialize(options, quantities, reporter);

  // Set QP matrix as pointer to approximate Hessian
  qp_solver_->setMatrix(symmetric_matrix_);
  qp_solver_termination_->setMatrix(symmetric_matrix_termination_);

} // end initialize

// Set iteration header
void Strategies::setIterationHeader()
{

  // Set iteration header string based on strategy objects
  iteration_header_ = "";
  if (direction_computation_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += direction_computation_->iterationHeader();
  } // end if
  if (termination_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += termination_->iterationHeader();
  } // end if
  if (line_search_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += line_search_->iterationHeader();
  } // end if
  if (approximate_hessian_update_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += approximate_hessian_update_->iterationHeader();
  } // end if
  if (point_set_update_->iterationHeader().length() > 0) {
    iteration_header_ += " ";
    iteration_header_ += point_set_update_->iterationHeader();
  } // end if

} // end setIterationHeader

// Print header
void Strategies::printHeader(const Reporter* reporter)
{

  // Print header
  reporter->printf(R_NL, R_BASIC, "Approximate Hessian update strategy.. : %s\n"
                                  "Derivative checker strategy.......... : %s\n"
                                  "Direction computation strategy....... : %s\n"
                                  "Line search strategy................. : %s\n"
                                  "Point set update strategy............ : %s\n"
                                  "QP solver strategy................... : %s\n"
                                  "Symmetric matrix strategy............ : %s\n"
                                  "Termination strategy................. : %s\n",
                   approximate_hessian_update_->name().c_str(),
                   derivative_checker_->name().c_str(),
                   direction_computation_->name().c_str(),
                   line_search_->name().c_str(),
                   point_set_update_->name().c_str(),
                   qp_solver_->name().c_str(),
                   symmetric_matrix_->name().c_str(),
                   termination_->name().c_str());

} // end printHeader

// Print footer
void Strategies::printFooter(const Reporter* reporter) {}

} // namespace NonOpt
