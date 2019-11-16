// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis, Baoyu Zhou

#ifndef __NONOPTQPSOLVERACTIVESET_HPP__
#define __NONOPTQPSOLVERACTIVESET_HPP__

#include <deque>

#include "NonOptQPSolver.hpp"

namespace NonOpt
{

/**
 * QPSolverActiveSet class
 */
class QPSolverActiveSet : public QPSolver
{

 public:
  /** @name Constructor */
  //@{
  /**
   * Constructor
   */
  QPSolverActiveSet();
  //@}

  /** @name Destructor */
  //@{
  /**
   * Destructor
   */
  ~QPSolverActiveSet();
  //@}

  /** @name Options handling methods */
  //@{
  /**
   * Add options
   * \param[in,out] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void addOptions(Options* options,
                  const Reporter* reporter);
  /**
   * Set options
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void setOptions(const Options* options,
                  const Reporter* reporter);
  //@}

  /** @name Initialize method */
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
   * Get dual step
   * \param[out] vector is dual step given by "-W*(G*omega + gamma)"
   */
  void dualStep(double vector[]);
  /**
   * Get feasible dual step
   * \param[out] vector is dual step given by "-W*(G*omega + gamma)" projected onto feasible region
   */
  void dualStepFeasible(double vector[]);
  /**
   * Get dual step's infinity norm
   * \return "||W*(G*omega + gamma)||_inf"
   */
  double dualStepNormInf();
  /**
   * Get feasible dual step's infinity norm
   * \return "||d||_inf" where d is projection of -W*(G*omega + gamma) onto feasible region
   */
  double dualStepFeasibleNormInf();
  /**
   * Get objective quadratic value
   * \return (G*omega + gamma)'*W*(G*omega + gamma)
   */
  double objectiveQuadraticValue();
  /**
   * Get objective quadratic value
   * \return "d'*H*d" where d is projection of -W*(G*omega + gamma) onto feasible region
   */
  double objectiveQuadraticValueFeasible();
  /**
   * Get gamma
   * \param[out] vector is "gamma" solution value
   */
  void gamma(double vector[]);
  /**
   * Get KKT error
   * \return solver's KKT error corresponding to current solution (ignores certain conditions assumed to be satisfied during solve)
   */
  double KKTError() { return kkt_error_; };
  /**
   * Get KKT error full
   * \return full KKT error corresponding to current solution
   */
  double KKTErrorFull();
  /**
   * Get name of strategy
   * \return string with name of strategy
   */
  std::string name() { return "ActiveSet"; };
  /**
   * Get iteration count
   * \return number of iterations performed
   */
  int numberOfIterations() { return iteration_count_; };
  /**
   * Get omega
   * \param[out] vector is "omega" solution value
   */
  void omega(double vector[]);
  //@}

  /** @name Set methods */
  //@{
  /**
   * Set dual step to zero
   */
  void setDualStepToZero()
  {
    dual_step_.scale(0.0);
    dual_step_feasible_.scale(0.0);
    dual_step_feasible_best_.scale(0.0);
  };
  /**
   * Set matrix
   * \param[in] matrix is pointer to SymmetricMatrix to be set as QP "W" data
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
  void setVectorList(const std::vector<std::shared_ptr<Vector> > vector_list) { vector_list_ = vector_list; };
  /**
   * Set vector
   * \param[in] vector is vector of double values to be set as QP "b" data
   */
  //void setVector(const std::vector<double> vector){vector_=vector; };
  void setVector(const std::vector<double> vector)
  {
    std::vector<double> vector2(vector.size(), 0.0);
    vector_ = vector2;
  };
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
  void addData(const std::vector<std::shared_ptr<Vector> > vector_list,
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
               const Reporter* reporter);
  /**
   * Solve QP hot, after new data added, re-using previous solution, factorization, etc.
   * \param[in] options is pointer to Options object from NonOpt
   * \param[in] reporter is pointer to Reporter object from NonOpt
   */
  void solveQPHot(const Options* options,
                  const Reporter* reporter);
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
  QPSolverActiveSet(const QPSolverActiveSet&);
  /**
   * Overloaded equals operator
   */
  void operator=(const QPSolverActiveSet&);
  //@}

  /** @name Private members */
  //@{
  /**
   * Length parameters
   */
  int factor_length_;
  int gamma_length_;
  int system_solution_length_;
  /**
   * Input parameters
   */
  bool fail_on_factorization_error_;
  bool allow_inexact_termination_;
  double cholesky_tolerance_;
  double kkt_tolerance_;
  double inexact_termination_factor_;
  double inexact_termination_ratio_min_;
  double linear_independence_tolerance_;
  int iteration_limit_minimum_;
  int iteration_limit_maximum_;
  /**
   * QP data quantities
   */
  double scalar_;                                     /**< "r" */
  std::shared_ptr<SymmetricMatrix> matrix_;           /**< "W" */
  std::vector<std::shared_ptr<Vector> > vector_list_; /**< "G" */
  std::vector<double> vector_;                        /**< "b" */
  /**
   * Algorithm parameters
   */
  int iteration_count_;
  double kkt_error_;
  double dual_objective_best_;
  double primal_objective_reference_;
  double dual_objective_simple_;
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
  std::deque<int> omega_positive_;
  std::deque<int> omega_positive_best_;
  double* system_solution_;
  double* system_solution_best_;
  /**
   * Solution quantities
   */
  Vector combination_;
  Vector combination_translated_;
  Vector dual_step_;
  double dual_step_projection_scalar_;
  Vector dual_step_feasible_;
  Vector dual_step_feasible_best_;
  Vector dual_step_simple_;
  Vector gamma_;
  double multiplier_;
  Vector omega_;
  //@}

  /** @name Private methods */
  //@{
  /**
   * Sanity check
   */
  bool checkQuantityCompatibility();
  /**
   * Update best solution
   */
  bool updateBestSolution();
  /**
   * Termination condition for inexact solution
   */
  bool inexactTerminationCondition(const Reporter* reporter);
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
  void evaluateDualVectors();
  void evaluateDualMultiplier(double solution1[],
                              double solution2[]);
  void evaluateSystemVector(int set,
                            int index,
                            double system_vector[]);
  void finalizeSolution();
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
  //@}

};  // end QPSolverActiveSet

}  // namespace NonOpt

#endif /* __NONOPTQPSOLVERACTIVESET_HPP__ */
