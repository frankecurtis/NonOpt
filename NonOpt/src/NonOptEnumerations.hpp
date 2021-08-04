// Copyright (C) 2019 Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Frank E. Curtis

#ifndef __NONOPTENUMERATIONS_HPP__
#define __NONOPTENUMERATIONS_HPP__

namespace NonOpt
{

/** @name Enumerations */
//@{
/**
 * NonOpt enumerations
 */
enum NonOpt_Status
{
  NONOPT_UNSET = -1,
  NONOPT_SUCCESS,
  NONOPT_OBJECTIVE_SIMILARITY,
  NONOPT_OBJECTIVE_TOLERANCE,
  NONOPT_CPU_TIME_LIMIT,
  NONOPT_ITERATE_NORM_LIMIT,
  NONOPT_ITERATION_LIMIT,
  NONOPT_FUNCTION_EVALUATION_LIMIT,
  NONOPT_GRADIENT_EVALUATION_LIMIT,
  NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE,
  NONOPT_DERIVATIVE_CHECKER_FAILURE,
  NONOPT_DIRECTION_COMPUTATION_FAILURE,
  NONOPT_FUNCTION_EVALUATION_FAILURE,
  NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE,
  NONOPT_GRADIENT_EVALUATION_FAILURE,
  NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE,
  NONOPT_LINE_SEARCH_FAILURE,
  NONOPT_POINT_SET_UPDATE_FAILURE,
  NONOPT_PROBLEM_DATA_FAILURE,
  NONOPT_SYMMETRIC_MATRIX_ASSERT_FAILURE,
  NONOPT_TERMINATION_FAILURE,
  NONOPT_VECTOR_ASSERT_FAILURE
};
/**
 * Approximate Hessian update enumerations
 */
enum AH_Status
{
  AH_UNSET = -1,
  AH_SUCCESS,
  AH_EVALUATION_FAILURE,
  AH_NORM_TOLERANCE_VIOLATION,
  AH_PRODUCT_TOLERANCE_VIOLATION
};
/**
 * Derivative checker enumerations
 */
enum DE_Status
{
  DE_UNSET = -1,
  DE_SUCCESS,
  DE_FAILURE
};
/**
 * Direction computation enumerations
 */
enum DC_Status
{
  DC_UNSET = -1,
  DC_SUCCESS,
  DC_CPU_TIME_LIMIT,
  DC_ITERATION_LIMIT,
  DC_EVALUATION_FAILURE,
  DC_QP_FAILURE
};
/**
 * Line search enumerations
 */
enum LS_Status
{
  LS_UNSET = -1,
  LS_SUCCESS,
  LS_ITERATION_LIMIT,
  LS_EVALUATION_FAILURE,
  LS_STEPSIZE_TOO_SMALL,
  LS_INTERVAL_TOO_SMALL
};
/**
 * Point set update enumerations
 */
enum PS_Status
{
  PS_UNSET = -1,
  PS_SUCCESS,
  PS_FAILURE
};
/**
 * QP solver enumerations
 */
enum QP_Status
{
  QP_UNSET = -1,
  QP_SUCCESS,
  QP_CPU_TIME_LIMIT,
  QP_ITERATION_LIMIT,
  QP_FACTORIZATION_ERROR,
  QP_INPUT_ERROR,
  QP_NAN_ERROR
};
/**
 * Report type enumerations
 */
enum ReportType
{
  R_NL = 0,
  R_QP
};
/**
 * Report level enumerations
 */
enum ReportLevel
{
  R_NONE = 0,
  R_BASIC,
  R_PER_ITERATION,
  R_PER_INNER_ITERATION
};
/**
 * Symmetric matrix enumerations
 */
enum SM_Status
{
  SM_UNSET = -1,
  SM_SUCCESS,
  SM_FAILURE
};
/**
 * Termination enumerations
 */
enum TE_Status
{
  TE_UNSET = -1,
  TE_SUCCESS,
  TE_EVALUATION_FAILURE
};
//@}

} // namespace NonOpt

#endif /* __NONOPTENUMERATIONS_HPP__ */
