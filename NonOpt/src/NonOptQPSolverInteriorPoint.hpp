// Copyright (C) 2025 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Lara Zebiane

#ifndef NONOPTQPSOLVERINTERIORPOINT_HPP_
#define NONOPTQPSOLVERINTERIORPOINT_HPP_

#include <deque>

#include "NonOptQPSolver.hpp"
#include "NonOptSymmetricMatrixDense.hpp"

namespace NonOpt
{

/**
 * QPSolverInteriorPoint class
 */
class QPSolverInteriorPoint : public QPSolver
{

public:
  /** @name Constructor */
  //@{
  /**
   * Constructor
   */
  QPSolverInteriorPoint();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~QPSolverInteriorPoint();
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   */
  void addOptions(Options* options);
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   */
  void setOptions(Options* options);
  //@}

  /** @name Initialization method */
  //@{
  /**
   * Initialize strategy
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in,out] quantities is pointer to Quantities object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void initialize(const Options* options,
                  Quantities* quantities,
                  const Reporter* reporter);
  /**
   * Initialize data
   * \param[in] gamma_length is length of gamma solution vector
   */
  void initializeData(int gamma_length);
  //@}

  /** @name Get methods */
  //@{
  /**
   * Get KKT error
   * \return solver's KKT error
   */
  double KKTError() { return kkt_error_; };
  /**
   * Get KKT error full
   * \return full KKT error corresponding to dual solution
   */
  double KKTErrorDual() { return kkt_error_; };  
  /**
   * Get iteration count
   * \return number of iterations performed
   */
  int numberOfIterations() { return iteration_count_; };  //change to the outer iteration NO sum of iteration :)
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "InteriorPoint"; };
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set inexact solution tolerance
   */
  void setInexactSolutionTolerance(double tolerance) { inexact_solution_tolerance_ = tolerance; };
  /**
   * Set "d" to zero
   */
  void setPrimalSolutionToZero()
  {
    primal_solution_.scale(0.0);
    primal_solution_feasible_.scale(0.0);
    primal_solution_feasible_best_.scale(0.0);
  };
  /**
   * Set matrix
   * \param[in] matrix is pointer to SymmetricMatrix, for which "W" is the "Inverse"
   */
  void setMatrix(const std::shared_ptr<SymmetricMatrix> matrix) { matrix_ = matrix; };
  /**
   * Set null solution
   */
  void setNullSolution();
  /**
   * Set vector list
   * \param[in] vector_list is vector of pointers to Vectors to be set as QP "G" data
   */
  void setVectorList(const std::vector<std::shared_ptr<Vector>> vector_list) { vector_list_ = vector_list; };
  /**
   * Set vector
   * \param[in] vector is vector of double values to be set as QP "b" data
   */
  void setVector(const std::vector<double> vector) { vector_ = vector; };
  /**
   * Set scalar
   * \param[in] scalar is double value to be set as QP "r" data
   */
  void setScalar(double scalar) { scalar_ = scalar; };
  //@}

  /** @name Add data methods */
  //@{
  /**
   * Add vectors
   * \param[in] vector_list is vector of pointers to Vectors to be added to QP "G" data
   * \param[in] vector is vector of double values to be added to QP "b" data
   */
  void addData(const std::vector<std::shared_ptr<Vector>> vector_list,
               const std::vector<double> vector);
  //@}

  /** @name Solve methods */
  //@{
  /**
   * Solve QP
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void solveQP(const Options* options,
               const Reporter* reporter,
               Quantities* quantities);
  /**
   * Solve QP hot, after new data added, re-using previous solution, factorization, etc.
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void solveQPHot(const Options* options,
                  const Reporter* reporter,
                  Quantities* quantities);
  //@}

  /** @name Print methods */
  //@{
  /**
   * Print data
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void printData(const Reporter* reporter);
  //@}

private:
  /** @name Default compiler generated methods
   * (Hidden to avoid implicit creation/calling.)
   */
  //@{
  /**
   * Copy constructor
   */
  QPSolverInteriorPoint(const QPSolverInteriorPoint&) = delete;
  /**
   * Overloaded equals operator
   */
  void operator=(const QPSolverInteriorPoint&) = delete;
  //@}

  /** @name Private members */
  //@{
  /**
   * Input parameters
   */
  bool allow_inexact_termination_;
  double barrier_parameter_factor_;
  double barrier_parameter_initial_;
  double barrier_parameter_maximum_;
  double barrier_parameter_minimum_;
  double fraction_to_boundary_tolerance_;
  double kkt_tolerance_;
  double inexact_solution_tolerance_;
  double inexact_termination_descent_tolerance_;
  double inexact_termination_initialization_factor_;
  double inexact_termination_ratio_minimum_;
  double solution_initialization_factor_;
  int inexact_termination_check_interval_;
  int iteration_limit_;

  /**
   * Interior-point system quantities
   */
  Vector At_;
  double b_;
  SymmetricMatrixDense Q_;
  Vector c_;
  Vector r_dual_;
  Vector r_comp_;
  double r_prim_;
  Vector theta_;
  double u_; 
  Vector v_;
  double mu_;

  /**
   * Cached LAPACK workspace for solveLinearSystem, re-queried/reallocated
   * only when the system size changes instead of on every solve
   */
  double* linear_system_workspace_;
  int linear_system_lwork_;
  int linear_system_workspace_size_;

  /**
   * Algorithm parameters
   */
  bool predict_boundary_;
  int iteration_count_;
  double kkt_error_;
  double dual_objective_reference_;
  double primal_directional_derivative_feasible_best_;
  double primal_objective_feasible_best_;
  double primal_objective_reference_;
  double primal_objective_simple_;
  double primal_quadratic_feasible_best_;
  double primal_solution_feasible_best_norm_inf_;

  /**
   * Solution quantities
   */
  Vector Gomega_;
  Vector primal_solution_feasible_best_;
  Vector primal_solution_simple_;
  //@}

  /** @name Private methods */
  //@{
  /**
   * Termination condition for inexact solution
   */
  bool inexactTerminationCondition(const Quantities* quantities,
                                   const Reporter* reporter);
  /**
   * Evaluate primal vectors
   */
  void evaluatePrimalVectors();
  /**
   * Finalize solution
   */
  void finalizeSolution();
  /**
   * Solve linear system
   */
  void solveLinearSystem(int size,
                         double matrix[],
                         double right_hand_side[],
                         int ipiv[],
                         double solution[]);
  void solveLinearSystemReuseFactorization(int size,
                                           double matrix[],
                                           double right_hand_side[],
                                           int ipiv[],
                                           double solution[]);
  /**
   * Compute step sizes
   */
  void computeStepSizes(const Vector& dtheta,
                        const double& du,
                        const Vector& dv,
                        double& atheta,
                        double& au,
                        double& av);
  //@}

}; // end QPSolverInteriorPoint

} // namespace NonOpt

#endif /* NONOPTQPSOLVERINTERIORPOINT_HPP_ */
