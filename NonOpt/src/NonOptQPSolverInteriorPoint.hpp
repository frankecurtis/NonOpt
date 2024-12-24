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
  double KKTError() { return kkt_error_; };  //residual r
  /**
   * Get KKT error full
   * \return full KKT error corresponding to dual solution
   */
  double KKTErrorDual() { return kkt_error_; };  
  /**
   * Get iteration count
   * \return number of iterations performed
   */
  int numberOfIterations() { return sum_inner_iter_count_; };  //change to the outer iteration NO sum of iteration :)
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
                  Quantities* quantities){
                    solveQP(options,reporter,quantities
                    );
                  }
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
  int factor_length_;
  int gamma_length_; /*n
  */
  int system_solution_length_;
  /**
   * Input parameters
   */
  bool fail_on_factorization_error_;
  bool allow_inexact_termination_;
  double cholesky_tolerance_;
  double kkt_tolerance_;
  double inexact_solution_tolerance_;
  double inexact_termination_descent_tolerance_;
  double inexact_termination_initialization_factor_;
  double inexact_termination_ratio_minimum_;
  double linear_independence_tolerance_;
  int inexact_termination_check_interval_;
  int iteration_limit_minimum_;
  int iteration_limit_maximum_;
  /**
   * Lara
   */
  double beta_;
  double sigma_;
  double eps_;
  double tol_in_;
  double tol_out_;
  double mu_factor_;
  int max_iter_in_;
  int max_iter_out_;

  /**
   * QP data quantities
   */
  double scalar_;                                    /**< "r" (delta) */ 
  std::shared_ptr<SymmetricMatrix> matrix_;          /**< "H" */
  std::vector<std::shared_ptr<Vector>> vector_list_; /**< "G" */
  std::vector<double> vector_;                       /**< "b_k" */
  /**
   * Lara
   */
  Vector matrix_At_;  /**< "A transpose" */
  double scalar_b_;
  SymmetricMatrixDense matrix_Q_;
  Vector vector_c_;
  SymmetricMatrixDense matrix_J_;
  Vector r_dual_;
  Vector r_cent_;
  double r_pri_;
  Vector vector_x_;
  double scalar_y_; 
  Vector vector_z_;
  double mu_;
  Vector primal_solution_; 
  double scalar_tau_;


  int inner_iter_count_;
  int outer_iter_count_;
  int sum_inner_iter_count_;

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
   * Algorithm quantities
   */
  double* factor_;
  std::deque<int> gamma_negative_;
  std::deque<int> gamma_negative_best_;
  std::deque<int> gamma_positive_;
  std::deque<int> gamma_positive_best_;
  double* inner_solution_1_;
  double* inner_solution_2_;
  double* inner_solution_3_;
  double* inner_solution_ls_;
  double* inner_solution_trial_;
  double* new_system_vector_;
  double* right_hand_side_;
  std::deque<int> omega_positive_;
  std::deque<int> omega_positive_best_;
  double* system_solution_;
  double* system_solution_best_;

  /**
   * Solution quantities
   */
  double multiplier_;
  double primal_solution_projection_scalar_;
  Vector combination_;
  Vector combination_translated_;
  Vector gamma_;
  Vector Gomega_;
  Vector omega_;
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
   * Update best solution
   */
  bool updateBestSolution();
  /**
   * Termination condition for inexact solution
   */
  bool inexactTerminationCondition(const Quantities* quantities,
                                   const Reporter* reporter);
  /**
   * Internal solve methods
   */
  bool choleskyAugment(double system_vector[],
                       int index,
                       double solution1[],
                       double value1,
                       double solution2[],
                       double value2);
  void choleskyDelete(int index,
                      double solution1[],
                      double solution2[]);
  void choleskyFromScratch(const Reporter* reporter);
  void evaluatePrimalVectors();
  void evaluatePrimalMultiplier(double solution1[],
                                double solution2[]);
  void evaluateSystemVector(int set,
                            int index,
                            double system_vector[]);
  void finalizeSolution(const Reporter* reporter);
  void resizeSystemSolution();
  bool setAugment(const Reporter* reporter,
                  int set,
                  int index,
                  double system_vector[],
                  double solution1[],
                  double solution2[],
                  double augmentation_value);
  void setDelete(const Reporter* reporter,
                 int set,
                 int index,
                 double solution1[],
                 double solution2[]);
  void solveSystem(double right_hand_side[],
                   double solution[]);
  void solveSystemTranspose(double right_hand_side[],
                            double solution[]);
  void solveLinearSystem(int size,
                         double matrix[],
                         double right_hand_side[],
                         double solution[]);
  void calculate_step_sizes(const Vector& delta_x, 
                            const double& delta_y, 
                            const Vector& delta_z, 
                            const double& r_pri_, 
                            const Vector& r_dual_, 
                            const Vector& vector_x_, 
                            const Vector& vector_z_,
                            const Vector& vector_c_, 
                            SymmetricMatrixDense& matrix_Q_, 
                            const Vector& matrix_At_, 
                            double scalar_y_, 
                            double scalar_b_, 
                            double beta_,  // Assuming beta is another parameter needed
                            double& alpha_x, 
                            double& alpha_y, 
                            double& alpha_z);
  //@}

}; // end QPSolverInteriorPoint

} // namespace NonOpt

#endif /* __NONOPTQPSOLVERINTERIORPOINT_HPP__ */
