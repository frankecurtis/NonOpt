import nonoptpy as nonopt
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# ------------------ Problem Definition ------------------
class NoisyNonsmooth(nonopt.Problem):
    def __init__(self, eps_f, rng, eps_g=None):
        super().__init__()
        if eps_g is None:
            eps_g = 10*np.sqrt(eps_f)
        self.eps_f = eps_f
        self.eps_g = eps_g
        self.rng = rng

        self.true_f_history = []
        self.x_history = [] 
        self.last_recorded_iter = -1
        self.noise_cache = {}
        self.grad_noise_cache = {}
        self.solver = None

    def _get_true_values(self, x):
        residual = x[1] - 2.0 * x[0]**2 + 1.0
        f = (1.0 - x[0])**2 + 10.0 * np.abs(residual)

        if residual > 0:
            sign = 1.0
        elif residual < 0:
            sign = -1.0
        else:
            sign = 0.0

        g = np.array([
            -2.0 * (1.0 - x[0]) + 10.0 * sign * (-4.0 * x[0]),
            10.0 * sign
        ])
        return f, g

    def numberOfVariables(self):
        return 2

    def initialPoint(self, n, x):
        x[0], x[1] = -1.2, 1.0
        return True

    def evaluateObjective(self, x):
        f_true, _ = self._get_true_values(x)
        key = tuple(np.round(x, 12))
        if key not in self.noise_cache:
            self.noise_cache[key] = self.rng.uniform(-self.eps_f, self.eps_f)
        return f_true + self.noise_cache[key]

    def evaluateGradient(self, n, x, g):
        f_true, g_true = self._get_true_values(x)
        key = tuple(np.round(x, 12))

        if key not in self.grad_noise_cache:
            noise = self.rng.normal(size=2)
            norm = np.linalg.norm(noise)
            if norm > 0:
                noise = (noise / norm) * self.rng.uniform(0, self.eps_g)
            else:
                noise = np.zeros(2)
            self.grad_noise_cache[key] = noise

        noise = self.grad_noise_cache[key]
        g[0] = g_true[0] + noise[0]
        g[1] = g_true[1] + noise[1]

        if self.solver is not None:
            current_iter = self.solver.iterations()
            if current_iter > self.last_recorded_iter:
                self.true_f_history.append(f_true)
                self.x_history.append(x.copy()) 
                self.last_recorded_iter = current_iter
        return True

    def finalizeSolution(self, n, x, f, g):
        f_true, _ = self._get_true_values(x)
        self.true_f_history.append(f_true)
        self.x_history.append(x.copy()) 
        return True


# ------------------ Execution ------------------
eps_f_values = [0.1, 0.01, 0.001, 0.0001]
base_seed = 42
histories = []

for idx, eps_f in enumerate(eps_f_values):
    rng = np.random.default_rng(base_seed + idx * 10)
    prob = NoisyNonsmooth(eps_f, rng)
    solver = nonopt.Solver()
    prob.solver = solver

    opts = solver.options()
    opts.modifyIntegerValue("print_level", 2)
    opts.modifyDoubleValue("LSB_stepsize_sufficient_decrease_fudge_factor", 2.1 * eps_f)
    opts.modifyDoubleValue("stationarity_tolerance", 1e-6)
    opts.modifyStringValue("direction_computation", "GradientCombination")
    # opts.modifyBoolValue("DCGC_try_gradient_step", 0)
    # opts.modifyBoolValue("DCGC_try_shortened_step", 0)
    # opts.modifyIntegerValue("TB_objective_similarity_limit", 100)
    # opts.modifyDoubleValue("DCGC_random_sample_factor", 100)
    # opts.modifyDoubleValue("DCGC_gradient_stepsize", 1.0)
    # opts.modifyDoubleValue("DCGC_shortened_stepsize", 1.0)

    try:
        solver.optimize(prob)
    except Exception as e:
        print(f"Solver error for eps_f={eps_f}:", e)
    
    histories.append((eps_f, prob.true_f_history, np.array(prob.x_history)))

# --- FIGURE 1: CONVERGENCE (TRUE OBJECTIVE) ---
all_f_values = [val for _, f_hist, _ in histories for val in f_hist]
all_eps_values = [eps for eps, _, _ in histories]
y_min = min(min(all_f_values), min(all_eps_values)) * 0.5 
y_max = max(all_f_values) * 2.0

fig1, axs1 = plt.subplots(2, 2, figsize=(12, 8))
# fig1.suptitle("TObjective Function Convergence (Shared Y-Axis)", fontsize=14)
axs1 = axs1.flatten()

for idx, (eps_f, f_hist, _) in enumerate(histories):
    ax = axs1[idx]
    ax.plot(f_hist, color="#3182bd", linewidth=2, label="True $f(x_k)$")
    ax.axhline(y=eps_f, color="black", linestyle=":", alpha=0.8, label=rf"$\epsilon_f = {eps_f}$")
    ax.set_yscale("log")
    ax.set_ylim(y_min, y_max)
    ax.set_title(rf"$\epsilon_f={eps_f}$")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("True $f(x)$")
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    ax.legend(loc='upper right', fontsize=8)

fig1.tight_layout()

# --- FIGURE 2: X-TRAJECTORIES (x0 vs x1) ---
all_x0 = np.concatenate([h[2][:, 0] for h in histories])
all_x1 = np.concatenate([h[2][:, 1] for h in histories])
x0_min, x0_max = all_x0.min(), all_x0.max()
x1_min, x1_max = all_x1.min(), all_x1.max()

pad0 = (x0_max - x0_min) * 0.1
pad1 = (x1_max - x1_min) * 0.1

fig2, axs2 = plt.subplots(2, 2, figsize=(12, 10))
# fig2.suptitle("Optimization Trajectories (Shared X/Y-Axes)", fontsize=14)
axs2 = axs2.flatten()

for idx, (eps_f, _, x_hist) in enumerate(histories):
    ax = axs2[idx]
    ax.plot(x_hist[:, 0], x_hist[:, 1], '-o', markersize=3, color="#e6550d", alpha=0.6, label="Path")
    ax.plot(x_hist[0, 0], x_hist[0, 1], 'ks', markersize=8, label="Start")
    ax.plot(x_hist[-1, 0], x_hist[-1, 1], 'r*', markersize=12, label="End")
    
    ax.set_xlim(x0_min - pad0, x0_max + pad0)
    ax.set_ylim(x1_min - pad1, x1_max + pad1)
    ax.set_title(rf"$\epsilon_f={eps_f}$")
    ax.set_xlabel("$x_0$")
    ax.set_ylabel("$x_1$")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(fontsize=8, loc='upper right')

fig2.tight_layout()

# This displays both windows
plt.show()