// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#include "NonOptStrategies.hpp"
#include "NonOptApproximateHessianUpdateBFGS.hpp"
#include "NonOptApproximateHessianUpdateDFP.hpp"
#include "NonOptDirectionComputationCuttingPlane.hpp"
#include "NonOptDirectionComputationGradient.hpp"
#include "NonOptDirectionComputationGradientCombination.hpp"
#include "NonOptLineSearchBacktracking.hpp"
#include "NonOptLineSearchWeakWolfe.hpp"
#include "NonOptPointSetUpdateProximity.hpp"
#include "NonOptQPSolverDualActiveSet.hpp"
#include "NonOptSymmetricMatrixDense.hpp"
#include "NonOptSymmetricMatrixLimitedMemory.hpp"
#include "NonOptTerminationGradientCombination.hpp"

namespace NonOpt
{

// Add options
void Strategies::addOptions(Options* options,
                            const Reporter* reporter)
{

  // Add string options
  options->addStringOption(reporter,
                           "direction_computation",
                           "CuttingPlane",
                           "Direction computation strategy to use.\n"
                           "Default     : CuttingPlane.");
  options->addStringOption(reporter,
                           "approximate_hessian_update",
                           "BFGS",
                           "Approximate Hessian update strategy to use.\n"
                           "Default     : BFGS.");
  options->addStringOption(reporter,
                           "line_search",
                           "WeakWolfe",
                           "Line search strategy to use.\n"
                           "Default     : WeakWolfe.");
  options->addStringOption(reporter,
                           "point_set_update",
                           "Proximity",
                           "Point set update strategy to use.\n"
                           "Default     : Proximity.");
  options->addStringOption(reporter,
                           "qp_solver",
                           "DualActiveSet",
                           "QP solver strategy to use.\n"
                           "Default     : DualActiveSet.");
  options->addStringOption(reporter,
                           "symmetric_matrix",
                           "Dense",
                           "Symmetric matrix strategy to use.\n"
                           "Default     : Dense.");
  options->addStringOption(reporter,
                           "termination",
                           "GradientCombination",
                           "Termination strategy to use.\n"
                           "Default     : GradientCombination.");

  // Add options for direction computation strategies
  std::shared_ptr<DirectionComputation> direction_computation;
  direction_computation = std::make_shared<DirectionComputationCuttingPlane>();
  direction_computation->addOptions(options, reporter);
  direction_computation = std::make_shared<DirectionComputationGradient>();
  direction_computation->addOptions(options, reporter);
  direction_computation = std::make_shared<DirectionComputationGradientCombination>();
  direction_computation->addOptions(options, reporter);
  // ADD NEW DIRECTION COMPUTATION STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for approximate Hessian update strategies
  std::shared_ptr<ApproximateHessianUpdate> approximate_hessian_update;
  approximate_hessian_update = std::make_shared<ApproximateHessianUpdateBFGS>();
  approximate_hessian_update->addOptions(options, reporter);
  approximate_hessian_update = std::make_shared<ApproximateHessianUpdateDFP>();
  approximate_hessian_update->addOptions(options, reporter);
  // ADD NEW APPROXIMATE HESSIAN UPDATE STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for line search strategies
  std::shared_ptr<LineSearch> line_search;
  line_search = std::make_shared<LineSearchBacktracking>();
  line_search->addOptions(options, reporter);
  line_search = std::make_shared<LineSearchWeakWolfe>();
  line_search->addOptions(options, reporter);
  // ADD NEW LINE SEARCH STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for point set update strategies
  std::shared_ptr<PointSetUpdate> point_set_update;
  point_set_update = std::make_shared<PointSetUpdateProximity>();
  point_set_update->addOptions(options, reporter);
  // ADD NEW POINT SET UPDATE STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for QP solver strategies
  std::shared_ptr<QPSolver> qp_solver;
  qp_solver = std::make_shared<QPSolverDualActiveSet>();
  qp_solver->addOptions(options, reporter);
  // ADD NEW QP SOLVER STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for symmetric matrix strategies
  std::shared_ptr<SymmetricMatrix> symmetric_matrix;
  symmetric_matrix = std::make_shared<SymmetricMatrixDense>();
  symmetric_matrix->addOptions(options, reporter);
  symmetric_matrix = std::make_shared<SymmetricMatrixLimitedMemory>();
  symmetric_matrix->addOptions(options, reporter);
  // ADD NEW SYMMETRIC MATRIX STRATEGIES HERE AND IN SWITCH BELOW //

  // Add options for termination strategies
  std::shared_ptr<Termination> termination;
  termination = std::make_shared<TerminationGradientCombination>();
  termination->addOptions(options, reporter);
  // ADD NEW TERMINATION STRATEGIES HERE AND IN SWITCH BELOW //

} // end addOptions

// Set options
void Strategies::setOptions(const Options* options,
                            const Reporter* reporter)
{

  // Declare strategy names
  std::string direction_computation_name;
  std::string approximate_hessian_update_name;
  std::string line_search_name;
  std::string point_set_update_name;
  std::string qp_solver_name;
  std::string symmetric_matrix_name;
  std::string termination_name;

  // Read integer options
  options->valueAsString(reporter, "direction_computation", direction_computation_name);
  options->valueAsString(reporter, "approximate_hessian_update", approximate_hessian_update_name);
  options->valueAsString(reporter, "line_search", line_search_name);
  options->valueAsString(reporter, "point_set_update", point_set_update_name);
  options->valueAsString(reporter, "qp_solver", qp_solver_name);
  options->valueAsString(reporter, "symmetric_matrix", symmetric_matrix_name);
  options->valueAsString(reporter, "termination", termination_name);

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
  }
  else {
    qp_solver_ = std::make_shared<QPSolverDualActiveSet>();
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

  // Set termination strategy
  if (termination_name.compare("GradientCombination") == 0) {
    termination_ = std::make_shared<TerminationGradientCombination>();
  }
  else {
    termination_ = std::make_shared<TerminationGradientCombination>();
  }

  // Set direction computation options
  direction_computation_->setOptions(options, reporter);

  // Set approximate Hessian update options
  approximate_hessian_update_->setOptions(options, reporter);

  // Set line search options
  line_search_->setOptions(options, reporter);

  // Set point set update options
  point_set_update_->setOptions(options, reporter);

  // Set QP solver options
  qp_solver_->setOptions(options, reporter);

  // Set symmetric matrix options
  symmetric_matrix_->setOptions(options, reporter);

  // Set termination options
  termination_->setOptions(options, reporter);

} // end setOptions

// Initialize
void Strategies::initialize(const Options* options,
                            Quantities* quantities,
                            const Reporter* reporter)
{

  // Initialize direction computation
  direction_computation_->initialize(options, quantities, reporter);

  // Initialize approximate Hessian update
  approximate_hessian_update_->initialize(options, quantities, reporter);

  // Initialize line search
  line_search_->initialize(options, quantities, reporter);

  // Initialize point set update
  point_set_update_->initialize(options, quantities, reporter);

  // Initialize QP data
  qp_solver_->initialize(options, quantities, reporter);

  // Initialize symmetric matrix
  symmetric_matrix_->initialize(options, quantities, reporter);

  // Initialize termination
  termination_->initialize(options, quantities, reporter);

  // Set QP matrix as pointer to approximate Hessian
  qp_solver_->setMatrix(symmetric_matrix_);

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
  reporter->printf(R_NL, R_BASIC, "Direction computation strategy....... : %s\n"
                                  "Approximate Hessian update strategy.. : %s\n"
                                  "Line search strategy................. : %s\n"
                                  "Point set update strategy............ : %s\n"
                                  "QP solver strategy................... : %s\n"
                                  "Symmetric matrix strategy............ : %s\n"
                                  "Termination strategy................. : %s\n",
                   direction_computation_->name().c_str(),
                   approximate_hessian_update_->name().c_str(),
                   line_search_->name().c_str(),
                   point_set_update_->name().c_str(),
                   qp_solver_->name().c_str(),
                   symmetric_matrix_->name().c_str(),
                   termination_->name().c_str());

} // end printHeader

// Print footer
void Strategies::printFooter(const Reporter* reporter) {}

} // namespace NonOpt
