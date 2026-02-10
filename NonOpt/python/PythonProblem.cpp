// Copyright (C) 2026 Lara Zebiane, Frank E. Curtis
//
// This code is published under the MIT License.
//
// Author(s) : Lara Zebiane

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "NonOptEnumerations.hpp"
#include "NonOptProblem.hpp"
#include "NonOptSolver.hpp"

namespace py = pybind11;

namespace NonOpt
{

/**
 * PyProblem class
 */
class PyProblem : public Problem
{

public:
  using Problem::Problem;

  /** @name Get methods */
  //@{
  /**
   * Number of variables
   * \param[out] n is the number of variables, an integer (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool numberOfVariables(int& n) override
  {
    auto func = py::get_override(this, "numberOfVariables");
    bool evaluation_success = false;
    if (func) {
      n = func().cast<int>(); 
      evaluation_success = true;
    }
    return evaluation_success;
  } // end numberOfVariables
  /**
   * Initial point
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[out] x is the initial point/iterate, a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool initialPoint(int n, double* x) override
  { 
    auto func = py::get_override(this, "initialPoint");
    bool evaluation_success = false;
    if (func) {
      auto x_arr = py::array_t<double>({n},{sizeof(double)},x,py::cast(this));
      evaluation_success = func(n, x_arr).cast<bool>();
    }
    return evaluation_success;
  }
  //@}

  /** @name Evaluate methods */
  //@{
  /**
   * Evaluates objective
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateObjective(int n, const double* x, double& f) override
  {
    auto func = py::get_override(this, "evaluateObjective");
    bool evaluation_success = false;
    if (func) {
      auto x_arr = py::array_t<double>({n},{sizeof(double)},const_cast<double*>(x));
      f = func(x_arr).cast<double>();
      evaluation_success = true;
    }
    return evaluation_success;
  }
  /**
   * Evaluates objective and gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] f is the objective value at "x", a double (return value)
   * \param[out] g is the gradient value at "x", a double array (return value)
   */
  bool evaluateObjectiveAndGradient(int n, const double* x, double& f, double* g) override
  {
    auto func = py::get_override(this, "evaluateObjectiveAndGradient");
    bool evaluation_success = false;
    if (func) {
      auto x_arr = py::array_t<double>({n},{sizeof(double)},const_cast<double*>(x));
      py::list f_val;
      f_val.append(0.0);
      auto g_arr = py::array_t<double>({n},{sizeof(double)},g,py::cast(this));
      evaluation_success = func(n, x_arr, f_val, g_arr).cast<bool>();
      f = f_val[0].cast<double>();
    }
    return evaluation_success;
  }
  /**
   * Evaluates gradient
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is a given point/iterate, a constant double array
   * \param[out] g is the gradient value at "x", a double array (return value)
   * \return indicator of success (true) or failure (false)
   */
  bool evaluateGradient(int n, const double* x, double* g) override
  {
    auto func = py::get_override(this, "evaluateGradient");
    bool evaluation_success = false;
    if (func) {
      auto x_arr = py::array_t<double>({n},{sizeof(double)},const_cast<double*>(x));
      auto g_arr = py::array_t<double>({n},{sizeof(double)},g,py::cast(this));
      evaluation_success = func(n, x_arr, g_arr).cast<bool>();
    }
    return evaluation_success;
  }
  //@}

  /** @name Finalize methods */
  //@{
  /**
   * Finalizes solution
   * \param[in] n is the number of variables, the size of "x", a constant integer
   * \param[in] x is the final point/iterate, a constant double array
   * \param[in] f is the objective value at "x", a constant double
   * \param[in] g is the gradient value at "x", a constant double array
   * \return indicator of success (true) or failure (false)
   */
  // Helper to convert raw double* to numpy without copying
  py::array_t<double> ptr_to_array(int n, double* ptr)
  {
    return py::array_t<double>({n}, {sizeof(double)}, ptr, py::handle());
  }
  bool finalizeSolution(int n, const double* x, double f, const double* g) override
  {
    auto x_arr = ptr_to_array(n, const_cast<double*>(x));
    auto g_arr = ptr_to_array(n, const_cast<double*>(g));
    PYBIND11_OVERRIDE_PURE(bool, Problem, finalizeSolution, n, x_arr, f, g_arr);
  }
  //@}

}; // end PyProblem

// pybind module
PYBIND11_MODULE(nonoptpy, m)
{

  py::enum_<NonOpt_Status>(m, "NonOpt_Status")
    .value("NONOPT_UNSET", NonOpt_Status::NONOPT_UNSET)
    .value("NONOPT_SUCCESS", NonOpt_Status::NONOPT_SUCCESS)
    .value("NONOPT_OBJECTIVE_SIMILARITY", NonOpt_Status::NONOPT_OBJECTIVE_SIMILARITY)
    .value("NONOPT_OBJECTIVE_TOLERANCE", NonOpt_Status::NONOPT_OBJECTIVE_TOLERANCE)
    .value("NONOPT_CPU_TIME_LIMIT", NonOpt_Status::NONOPT_CPU_TIME_LIMIT)
    .value("NONOPT_ITERATE_NORM_LIMIT", NonOpt_Status::NONOPT_ITERATE_NORM_LIMIT)
    .value("NONOPT_ITERATION_LIMIT", NonOpt_Status::NONOPT_ITERATION_LIMIT)
    .value("NONOPT_FUNCTION_EVALUATION_LIMIT", NonOpt_Status::NONOPT_FUNCTION_EVALUATION_LIMIT)
    .value("NONOPT_GRADIENT_EVALUATION_LIMIT", NonOpt_Status::NONOPT_GRADIENT_EVALUATION_LIMIT)
    .value("NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE", NonOpt_Status::NONOPT_APPROXIMATE_HESSIAN_UPDATE_FAILURE)
    .value("NONOPT_DERIVATIVE_CHECKER_FAILURE", NonOpt_Status::NONOPT_DERIVATIVE_CHECKER_FAILURE)
    .value("NONOPT_DIRECTION_COMPUTATION_FAILURE", NonOpt_Status::NONOPT_DIRECTION_COMPUTATION_FAILURE)
    .value("NONOPT_FUNCTION_EVALUATION_FAILURE", NonOpt_Status::NONOPT_FUNCTION_EVALUATION_FAILURE)
    .value("NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE", NonOpt_Status::NONOPT_FUNCTION_EVALUATION_ASSERT_FAILURE)
    .value("NONOPT_GRADIENT_EVALUATION_FAILURE", NonOpt_Status::NONOPT_GRADIENT_EVALUATION_FAILURE)
    .value("NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE", NonOpt_Status::NONOPT_GRADIENT_EVALUATION_ASSERT_FAILURE)
    .value("NONOPT_LINE_SEARCH_FAILURE", NonOpt_Status::NONOPT_LINE_SEARCH_FAILURE)
    .value("NONOPT_POINT_SET_UPDATE_FAILURE", NonOpt_Status::NONOPT_POINT_SET_UPDATE_FAILURE)
    .value("NONOPT_PROBLEM_DATA_FAILURE", NonOpt_Status::NONOPT_PROBLEM_DATA_FAILURE)
    .value("NONOPT_SYMMETRIC_MATRIX_ASSERT_FAILURE", NonOpt_Status::NONOPT_SYMMETRIC_MATRIX_ASSERT_FAILURE)
    .value("NONOPT_TERMINATION_FAILURE", NonOpt_Status::NONOPT_TERMINATION_FAILURE)
    .value("NONOPT_VECTOR_ASSERT_FAILURE", NonOpt_Status::NONOPT_VECTOR_ASSERT_FAILURE)
    .export_values();

  py::enum_<ReportType>(m, "ReportType")
    .value("R_NL", ReportType::R_NL)
    .value("R_QP", ReportType::R_QP);

  py::enum_<ReportLevel>(m, "ReportLevel")
    .value("R_NONE", ReportLevel::R_NONE)
    .value("R_BASIC", ReportLevel::R_BASIC)
    .value("R_PER_ITERATION", ReportLevel::R_PER_ITERATION)
    .value("R_PER_INNER_ITERATION", ReportLevel::R_PER_INNER_ITERATION);

  py::class_<Problem, PyProblem, std::shared_ptr<Problem>>(m, "Problem")
    .def(py::init<>());

  py::class_<NonOptSolver>(m, "Solver")
    .def(py::init<>())
    .def("functionEvaluations", &NonOptSolver::functionEvaluations)
    .def("gradientEvaluations", &NonOptSolver::gradientEvaluations)
    .def("iterations", &NonOptSolver::iterations)
    .def("numberOfVariables", &NonOptSolver::numberOfVariables)
    .def("totalInnerIterations", &NonOptSolver::totalInnerIterations)
    .def("totalQPIterations", &NonOptSolver::totalQPIterations)
    .def("objective", &NonOptSolver::objective)
    .def("stationarityRadius", &NonOptSolver::stationarityRadius)
    .def("time", &NonOptSolver::time)
    .def("timeEvaluations", &NonOptSolver::timeEvaluations)
    .def("timeNonOpt", &NonOptSolver::timeNonOpt)
    .def("status", &NonOptSolver::status)
    .def("optimize", &NonOptSolver::optimize)
    .def("options", &NonOptSolver::options, py::return_value_policy::reference)
    .def("reporter", &NonOptSolver::reporter, py::return_value_policy::reference)
    .def("solution", [](NonOptSolver &self, py::array_t<double> input) {
      py::buffer_info buff = input.request();
      if (buff.ndim != 1) {
        throw std::runtime_error("Input must be a 1D array");
      }
      self.solution(static_cast<double*>(buff.ptr));
    });

  py::class_<Options>(m, "Options")
    .def("modifyOptionsFromFile", &Options::modifyOptionsFromFile)
    .def("modifyBoolValue", &Options::modifyBoolValue)
    .def("modifyDoubleValue", &Options::modifyDoubleValue)
    .def("modifyIntegerValue", &Options::modifyIntegerValue)
    .def("modifyStringValue", &Options::modifyStringValue);

  py::class_<Reporter>(m, "Reporter")
    .def("addFileReport", &Reporter::addFileReport);

} // end pybind module

} // namespace NonOpt