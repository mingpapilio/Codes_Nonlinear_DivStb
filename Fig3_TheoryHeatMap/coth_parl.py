#!/usr/bin/env python3
"""
Overview
--------
This script explores how our general model behaves, which is a blended
generalised Lotka-Volterra (gLV) and θ-logistic model:
  1) the nonlinear crowding exponent alpha (theta-logistic curvature), and
  2) an interaction-similarity weight w that mixes two interaction channels
     (A_alpha, A_beta) derived from the same random matrix.

For each parameter triple (alpha, S, w), the code:
  • builds sparse interaction matrices A_alpha and A_beta from a base matrix,
  • integrates the ODEs to a pre-perturbation state,
  • applies a pulse dilution (multiplicative knock-down), then integrates again,
  • rejects runs with any extinctions, and
  • computes the leading real part of the Jacobian's eigenvalues at the terminal
    post-perturbation state as a local-stability proxy (more negative ⇒ faster
    return to equilibrium).

Model form:
  dN_i/dt = N_i * [ (r_alpha[i]/alpha) * (1 - (Σ_j A_alpha[i,j] N_j)^alpha)
                  + (r_beta[i]/beta) * (1 - (Σ_j A_beta[i,j] N_j)^beta) ]

Matrix construction:
  • Draw a sparse random base matrix A ~ Normal(mean=a_alpha_mean, sd=a_alpha_std),
    with off-diagonals present with probability Connectivity.
  • Set A_alpha = |A| * w with diag(A_alpha)=1 (self-regulation through the alpha-channel).
  • Set A_beta = |A| * (1-w) with diag(A_beta)=0 (no direct self-term through beta-channel).
    Thus w∈[0,1] shifts “how much” regulation flows through the alpha vs beta channel.

Outputs:
  For each (alpha, S, w) we save a row to result_array.txt:
    [alpha, N_species, w, leading_eigenvalue, extinction_counter, attempt_counter]
  where leading_eigenvalue from the first full-survival success (or NaN if none).

Example execution code:
  # from the script directory
  python3 coth_parl.py &
"""

# =========================
# Import Python packages
# =========================
import os
# Pin BLAS/OpenMP threads inside each worker to avoid oversubscription in parallel runs.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np
from scipy import integrate
import threading
from itertools import product
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# =========================
# Global options (sweep I/O)
# =========================
N_workers = 32          # Number of worker processes for ProcessPoolExecutor
precision = 1e-6        # atol/rtol for solve_ivp
max_attempts = 5        # Max retries to get a full-survival post-perturbation run
max_time = 200          # Wall-clock timeout (seconds) for each IVP solve

# =============
# ODE time grid
# =============
T_preturb = 500.0
T_posturb = 500.0
N_points_pre = int(T_preturb * 100)
N_points_post = int(T_posturb * 100)
t_pre = np.linspace(0., T_preturb, N_points_pre)
t_post = np.linspace(T_preturb, T_preturb + T_posturb, N_points_post)

# ===================
# Parameter settings
# ===================
Connectivity = 1.0      # Probability that an off-diagonal a[i,j] is present and sampled
dilution_coeff = 0.1    # Multiplicative pulse knock-down after the pre-perturbation phase
N_ini = 0.1              # Initial abundance for all species
beta = 1.0              # θ-exponent for the beta channel (held fixed here)
a_alpha_mean = 0.0      # Mean of the base random matrix A (Zero keeps |A| as half-normal; using a non-zero mean makes |A| folded normal. For non-zero means, use the supplementary code script.)
a_alpha_std = 0.05      # SD of the base random matrix A

# =============================
# Parameter grids for the sweep
# =============================
alpha_values = np.linspace(-2.375, 2.375, 20)   # θ-exponent grid (include sub- and super-linear)
N_species_values = [10, 25, 50, 100]            # Community sizes
w_values = np.linspace(0.0, 1.0, 21)            # Mix between A_alpha and A_beta

# ======================
# ODE: community dynamics
# ======================

def derivative(t, N, r_alpha, r_beta, a_alpha, a_beta, alpha, beta):
    """ Compute dN/dt for the blended θ-logistic (alpha-channel) + beta-channel model.

    Notes
    -----
    • We clip interaction sums to ε>0 to keep x^alpha well-defined for non-integer alpha.
    • Signs: because a_alpha and a_beta are made non-negative below, the model here
      behaves like purely competitive crowding (no mutualism). Drop abs() in the
      matrix construction to allow mixed signs.
    """
    N = np.array(N, dtype=float)  # Ensure a NumPy array (copy is cheap for small S)

    # Interaction pressures through the two channels
    interaction_alpha = np.dot(a_alpha, N)  # shape (S,)
    interaction_beta  = np.dot(a_beta,  N)  # shape (S,)

    # To avoid invalid powers / gradients when alpha is not an integer
    eps = 1e-12
    base_alpha = np.clip(interaction_alpha, eps, None)
    base_beta  = np.clip(interaction_beta,  eps, None)

    # Core θ-logistic style crowding in two channels
    growth = ((r_alpha / alpha) * (1.0 - np.power(base_alpha, alpha)) +
              (r_beta  / beta)  * (1.0 - np.power(base_beta,  beta  )))

    # Mass-action factor N_i (standard in gLV/θ-logistic formulations)
    dotN = N * growth
    return dotN

# ==============================
# Random sparse interaction matrix
# ==============================

def generate_interaction_matrix(num_species, connectivity, mean, sd):
    """ Generate a random interaction matrix with unit diagonal.

    Off-diagonals are realised with probability `connectivity`; values are Normal(mean, sd).
    The diagonal is set to 1.0 here, but downstream code will overwrite it depending on
    the alpha/beta channel design.
    """
    a = np.zeros((num_species, num_species), dtype=float)
    np.fill_diagonal(a, 1.0)
    for i in range(num_species):
        for j in range(num_species):
            if i != j and np.random.random() < connectivity:
                a[i, j] = np.random.normal(loc=mean, scale=sd)
    return a

# =============================================
# Timeout wrapper to guard stiff/slow ODE solves
# =============================================

class ODESolverWithTimeout:
    """ Run solve_ivp in a thread and join with a wall-time limit.

    If the solve exceeds max_time seconds, return None.
    """
    def __init__(self, fun, t_span, y0, args=(), max_time=5, precision=1e-6, max_step=0.1):
        self.fun = fun
        self.t_span = t_span
        self.y0 = y0
        self.args = args
        self.max_time = max_time
        self.precision = precision
        self.max_step = max_step
        self.result = None
        self.thread = None

    def solve(self):
        self.thread = threading.Thread(target=self._solve_ode)
        self.thread.start()
        self.thread.join(timeout=self.max_time)
        if self.thread.is_alive():
            print("Solver timed out")
            self.result = None
        return self.result

    def _solve_ode(self):
        try:
            sol = integrate.solve_ivp(
                self.fun,
                self.t_span,
                self.y0,
                args=self.args,
                method='BDF',           # Stiff solver; alternatives: 'Radau', 'LSODA' (Fortran via scipy)
                atol=self.precision,
                rtol=self.precision,
                max_step=self.max_step
            )
            self.result = sol
        except Exception as e:
            print(f"Solver failed: {e}")
            self.result = None

# ============================
# Linearisation / eigenvalues
# ============================

def compute_eigenvalues(N, r_alpha, r_beta, a_alpha, a_beta, alpha, beta):
    """ Approximate Jacobian eigenvalues at state N.

    We linearise f_i(N) = N_i * g_i(N), where
      g_i(N) = (r_alpha[i]/alpha) * (1 - (Σ_j A_alpha[i,j] N_j)^alpha)
       + (r_beta[i]/beta) * (1 - (Σ_j A_beta[i,j] N_j)^beta).

    The Jacobian is: J[i,j] = δ_ij * g_i(N) + N_i * ∂g_i/∂N_j.
    Using ∂/∂N_j (1 - x^alpha) = -alpha x^{alpha-1} ∂x/∂N_j with x = Σ_k A[i,k] N_k, the alpha,beta cancel:
      ∂g_i/∂N_j = - r_alpha[i] * A_alpha[i,j] * (Σ_k A_alpha[i,k] N_k)^{alpha-1}
                  - r_beta[i] * A_beta[i,j] * (Σ_k A_beta[i,k] N_k)^{beta-1}.

    Below we omit the δ_ij * g_i(N) term under the assumption that we are at (or
    extremely near) equilibrium so g_i≈0. If not guaranteed, consider adding it.
    """
    S = len(N)
    if np.isscalar(r_alpha):
        r_alpha = np.full(S, r_alpha, dtype=float)
    if np.isscalar(r_beta):
        r_beta = np.full(S, r_beta, dtype=float)

    sum_alpha = np.dot(a_alpha, N)  # shape (S,)
    sum_beta  = np.dot(a_beta,  N)  # shape (S,)

    # Optional safety to avoid undefined powers when alpha<1. Uncomment if needed.
    # eps = 1e-12
    # sum_alpha = np.clip(sum_alpha, eps, None)
    # sum_beta  = np.clip(sum_beta,  eps, None)

    J = np.zeros((S, S), dtype=float)
    for i in range(S):
        for j in range(S):
            term_alpha = r_alpha[i] * a_alpha[i, j] * (sum_alpha[i] ** (alpha - 1))
            term_beta  = r_beta[i]  * a_beta[i, j]  * (sum_beta[i]  ** (beta  - 1))
            # Omit δ_ij * g_i(N) term (assumes g_i≈0 at the evaluated state)
            J[i, j] = -N[i] * (term_alpha + term_beta)
    eigenvalues = np.linalg.eigvals(J)
    return eigenvalues

# ========================================
# Single (alpha, S, w) simulation with retries
# ========================================

def simulate_combination(params):
    """ Run one parameter triple with up to `max_attempts` tries to achieve zero extinctions.

    Returns a record: [alpha, N_species, w, leading_eigenvalue, extinction_counter, attempt_counter]
    leading_eigenvalue is the leading real part from the first full-survival success (or NaN if none).
    """
    alpha, N_species, w = params
    leading_eig = np.nan
    extinction_counter = 0

    for attempt in range(max_attempts):
        # Base interaction matrix (will be split into alpha and beta channels)
        a_base = generate_interaction_matrix(N_species, connectivity=Connectivity,
                                             mean=a_alpha_mean, sd=a_alpha_std)

        # Mix into alpha and beta channels; absolute value enforces non-negative entries.
        # CLARIFY: If mixed-sign interactions are desired, replace abs(a_base) with a_base.
        a_alpha = np.abs(a_base) * w
        a_beta  = np.abs(a_base) * (1 - w)

        # Channel-specific diagonals
        np.fill_diagonal(a_alpha, 1.0)   # self-regulation routed through alpha-channel
        np.fill_diagonal(a_beta,  0.0)   # no beta self-term

        # Species-specific intrinsic rates (can be vectors; here constants)
        r_alpha = np.repeat(1.0,  N_species)
        r_beta  = np.repeat(0.25, N_species)

        # ---------- Stage 1: pre-perturbation integration ----------
        N_initial = np.repeat(N_ini, N_species)
        solver1 = ODESolverWithTimeout(
            derivative,
            [0.0, T_preturb],
            y0=N_initial,
            args=(r_alpha, r_beta, a_alpha, a_beta, alpha, beta),
            max_time=max_time,
            precision=precision,
            max_step=0.1,
        )
        solved1 = solver1.solve()

        if solved1 is None or (not solved1.success):
            continue  # retry with a new draw

        # Last state of the adaptive integration
        N_mid = solved1.y[:, -1]

        # Apply dilution (pulse perturbation)
        N_perturb = N_mid * dilution_coeff

        # Guard against NaN/Inf instability
        if np.any(~np.isfinite(N_perturb)):
            continue

        # ---------- Stage 2: post-perturbation integration ----------
        solver2 = ODESolverWithTimeout(
            derivative,
            [T_preturb, T_preturb + T_posturb],
            y0=N_perturb,
            args=(r_alpha, r_beta, a_alpha, a_beta, alpha, beta),
            max_time=max_time,
            precision=precision,
            max_step=0.1,
        )
        solved2 = solver2.solve()
        if (solved2 is None) or (not solved2.success):
            continue

        N_post_point = solved2.y[:, -1]

        # Treat tiny abundances as extinct
        extinct_count = np.sum(N_post_point < 1e-8)
        if extinct_count > 0:
            extinction_counter += 1
            continue

        # Compute leading eigenvalue at the terminal state
        eigvals = compute_eigenvalues(N_post_point, r_alpha, r_beta,
                                      a_alpha, a_beta, alpha, beta)
        leading_val = np.max(np.real(eigvals))
        leading_eig = leading_val
        break
    return [alpha, N_species, w, leading_eig, extinction_counter, attempt + 1]

# ==============================
# Parallel sweep + save results
# ==============================
if __name__ == "__main__":
    # Build combinations (can stay outside, but keeping it here is tidy)
    param_combinations = list(product(alpha_values, N_species_values, w_values))

    results_list = []
    with ProcessPoolExecutor(max_workers=N_workers) as executor:
        futures = [executor.submit(simulate_combination, params) for params in param_combinations]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Parameter combinations"):
            results_list.append(future.result())

    result_array = np.array(results_list, dtype=float)
    np.savetxt(
        "result_array.txt",
        result_array,
        delimiter=",",
        header="alpha,N_species,w,leading_eigenvalue,extinction_counter,attempt_counter",
        comments="",
        fmt="%.6f",
    )
