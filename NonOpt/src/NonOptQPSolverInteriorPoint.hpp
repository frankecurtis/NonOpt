// Copyright (C) 2022 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Baoyu Zhou

#ifndef __NONOPTQPSOLVERINTERIORPOINT_HPP__
#define __NONOPTQPSOLVERINTERIORPOINT_HPP__

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
   * Get combination of vectors' infinity norm
   * \return "||G*omega||_inf"
   */
  double combinationNormInf();
  /**
   * Get translated combination of vectors' infinity norm
   * \return "||G*omega + gamma||_inf"
   */
  double combinationTranslatedNormInf();
  /**
   * Get translated combination of vectors' infinity norm
   * \return "||G*omega + gamma||_2^2"
   */
  double combinationTranslatedNorm2Squared();
  /**
   * Get dual objective quadratic value
   * \return "(G*omega + gamma)'*W*(G*omega + gamma)"
   */
  double dualObjectiveQuadraticValue();
  /**
   * Get dual solution
   * \param[out] omega is dual solution, omega part
   * \param[out] gamma is dual solution, gamma part
   */
  void dualSolution(double omega[], double gamma[]);
  /**
   * Get dual solution, omega part
   * \param[out] omega is dual solution, omega part
   */
  void dualSolutionOmega(double omega[]);
  /**
   * Get dual solution, omega part, length
   * \return dual solution, omega part, length
   */
  int dualSolutionOmegaLength() { return (int)vector_.size(); };
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
   * Get primal solution
   * \param[out] "d" (equal to "-W*(G*omega + gamma)" if solution is exact)
   */
  void primalSolution(double d[]);
  /**
   * Get primal solution that is feasible
   * \param[out] feasible primal solution
   */
  void primalSolutionFeasible(double d_feasible[]);
  /**
   * Get primal solution infinity norm
   * \return "||d||_inf"
   */
  double primalSolutionNormInf();
  /**
   * Get primal solution 2-norm square
   * \return "||d||_2^2"
   */
  double primalSolutionNorm2Squared();
  /**
   * Get feasible primal solution infinity norm
   * \return inf-norm of feasible primal solution
   */
  double primalSolutionFeasibleNormInf();
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "InteriorPoint"; };
  /**
   * Vector list length
   */
  inline int const vectorListLength() const { return (int)vector_list_.size(); };
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
  QPSolverInteriorPoint(const QPSolverInteriorPoint&);
  /**
   * Overloaded equals operator
   */
  void operator=(const QPSolverInteriorPoint&);
  //@}

  /** @name Private members */
  //@{
  /**
   * Length parameters
   */
  int gamma_length_;
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
   * QP data quantities
   */
  double scalar_;                                    /**< "r" (delta) */ 
  std::shared_ptr<SymmetricMatrix> matrix_;          /**< "H" */
  std::vector<std::shared_ptr<Vector>> vector_list_; /**< "G" */
  std::vector<double> vector_;                       /**< "b_k" */

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
   * Algorithm parameters
   */
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
  double primal_solution_projection_scalar_;
  Vector combination_;
  Vector combination_translated_;
  Vector gamma_;
  Vector Gomega_;
  Vector omega_;
  Vector primal_solution_;
  Vector primal_solution_feasible_;
  Vector primal_solution_feasible_best_;
  Vector primal_solution_simple_;
  //@}

  /** @name Private methods */
  //@{
  /**
   * Sanity check
   */
  bool checkQuantityCompatibility();
  /**
   * Dual objective quadratic value, scaled to correspond to feasible "d"
   */
  double dualObjectiveQuadraticValueScaled();
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

#endif /* __NONOPTQPSOLVERINTERIORPOINT_HPP__ */
