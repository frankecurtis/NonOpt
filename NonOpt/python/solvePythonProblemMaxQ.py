# Copyright (C) 2026 Lara Zebiane, Frank E. Curtis
#
# This code is published under the MIT License.
#
# Author(s) : Lara Zebiane, Frank E. Curtis

import nonoptpy as nonopt
import numpy as np

class MaxQ(nonopt.Problem):

  def numberOfVariables(self):
    return 1000
    
  def initialPoint(self, n, x):
    # set initial point
    for i in range(n//2):
      x[i] = i + 1
    for i in range(n//2,n):
      x[i] = -i - 1
    return True

  def evaluateObjective(self, x):
    # evaluate objective function value
    return np.max(np.abs(x))**2

  def evaluateObjectiveAndGradient(self, n, x, f, g):
    # evaluate objective function and gradient values
    f[0] = np.max(np.abs(x))**2
    index = np.argmax(np.abs(x))
    g.fill(0)
    g[index] = 2 * x[index]
    return True

  def evaluateGradient(self, n, x, g):
    # evaluate objective gradient value
    index = np.argmax(np.abs(x))
    g.fill(0)
    g[index] = 2 * x[index]
    return True

  def finalizeSolution(self, n, x, f, g):
    return True

# Declare algorithm
solver = nonopt.Solver()

# Modify options from file
solver.options().modifyOptionsFromFile("nonoptpy.opt")

# Add file report
solver.reporter().addFileReport("report","nonoptpy.out",nonopt.ReportType.R_NL,nonopt.ReportLevel.R_PER_ITERATION)

# Optimize
print("Running NonOpt from Python...")
solver.optimize(MaxQ())

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