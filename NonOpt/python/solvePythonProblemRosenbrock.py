# Copyright (C) 2026 Lara Zebiane, Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Lara Zebiane

import nonoptpy as nonopt
import numpy as np

class Rosenbrock(nonopt.Problem):

  def numberOfVariables(self):
    return 2
    
  def initialPoint(self, n, x):
    # set initial point
    x[0] = -1.2
    x[1] = 1.0
    return True

  def evaluateObjective(self, x):
    # evaluate objective function value
    return (1.0 - x[0])**2 + 100.0 * (x[1] - x[0]**2)**2

  def evaluateObjectiveAndGradient(self, n, x, f, g):
    # evaluate objective function and gradient values
    f[0] = (1.0 - x[0])**2 + 100.0 * (x[1] - x[0]**2)**2
    g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0]**2)
    g[1] = 200.0 * (x[1] - x[0]**2)
    return True

  def evaluateGradient(self, n, x, g):
    # evaluate objective gradient value
    g[0] = -2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0]**2)
    g[1] = 200.0 * (x[1] - x[0]**2)
    return True

  def finalizeSolution(self, n, x, f, g):
    return True

# Declare algorithm
solver = nonopt.Solver()

# Modify options from file
solver.options().modifyOptionsFromFile("nonoptpy.opt")

# Optimize
print("Running NonOpt from Python...")
solver.optimize(Rosenbrock())

# Print solver statistics
print("\nSolver statistics in Python...")
print(f"Number of function evaluations....... : {solver.functionEvaluations()}")
print(f"Number of gradient evaluations....... : {solver.gradientEvaluations()}")
print(f"Number of iterations................. : {solver.iterations()}")
print(f"Number of variables.................. : {solver.numberOfVariables()}")
print(f"Number of inner iterations........... : {solver.totalInnerIterations()}")
print(f"Number of QP iterations.............. : {solver.totalQPIterations()}")
print(f"Final objective function value....... : {solver.objective()}")
print(f"Final stationarity radius............ : {solver.stationarityRadius()}")
print(f"CPU seconds.......................... : {solver.time()}")
print(f"CPU seconds in evaluations........... : {solver.timeEvaluations()}")
print(f"CPU seconds in NonOpt................ : {solver.timeNonOpt()}")
print(f"Solver status........................ : {solver.status()}")

# Print solution
x = np.zeros(solver.numberOfVariables())
solver.solution(x)
print(x)