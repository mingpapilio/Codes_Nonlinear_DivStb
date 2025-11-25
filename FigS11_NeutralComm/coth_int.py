#!/usr/bin/env python3
"""
Overview (supplementary model with positive and negative interactions)
------------------------------------
This script implements the "integrated" formulation where interactions are split
by sign and routed through two channels (alpha and beta) before applying
nonlinear exponents, then assesses local stability after a pulse dilution.

Key differences from the main-text (half-normal) version:
  • We start from a full Normal base matrix A ~ Normal(mean=a_alpha_mean, sd=a_alpha_std).
  • We split by sign: negatives (competitive) go into a "positive bucket" by
    magnitude (|negatives|); nonnegative entries go into a "negative bucket". 
    Mixing across alpha/beta is via w.
  • Diagonals: the alpha_neg bucket gets diag = 1 (self-regulation); all others diag = 0.

Model form (per species i):
  dN_i/dt = N_i * [ (r_alpha[i]/alpha) * ( (1 + Σ_j A_alpha_pos[i,j] N_j)^alpha - (Σ_j A_alpha_neg[i,j] N_j)^alpha )
                  + (r_beta[i] / beta)  * ( (1 + Σ_j A_beta_pos[i,j]  N_j)^beta  - (Σ_j A_beta_neg[i,j]  N_j)^beta  ) ]

For each parameter triple (alpha, S, w), the code:
  • builds sign-split interaction matrices (A_alpha_pos, A_alpha_neg, A_beta_pos, A_beta_neg),
  • integrates the ODEs to a pre-perturbation state,
  • applies a pulse dilution (multiplicative knock-down), then integrates again,
  • rejects runs with any extinctions, and
  • computes the leading real part of the Jacobian‘s eigenvalues at the terminal
    post-perturbation state as a local-stability proxy (more negative ⇒ faster
    return to equilibrium).

Outputs:
  For each (alpha, S, w) we save a row to result_array.txt:
    [alpha, N_species, w, leading_eigenvalue, extinction_counter, attempt_counter]
  where leading_eigenvalue is from the first full-survival success (or NaN if none).

Example execution code:
  # from the script directory
  python3 coth_int.py &
"""

# =========================
# Import Python packages
# =========================
import os
# Pin BLAS/OpenMP threads inside each worker to avoid oversubscription in parallel runs.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

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
Connectivity = 1.0      # Probability that an off-diagonal a[i,j] is realised
dilution_coeff = 0.1    # Multiplicative pulse knock-down after the pre-perturbation phase
N_ini = 0.1             # Initial abundance for all species
beta = 1.0              # exponent for the beta channel (held fixed here)
a_alpha_mean = 0.05     # Mean of the base random matrix A (full Normal here; differs from half-normal in coth_parl)
a_alpha_std = 0.002     # SD of the base random matrix A

# =============================
# Parameter grids for the sweep
# =============================
alpha_values = np.linspace(-2.125, 2.125, 18)     # exponent grid (include sub- and super-linear)
N_species_values = [10, 25, 50, 100]            # Community sizes
Analysis_mode = 0                               # 0 = vary w  (current behaviour)
                                                # 1 = vary mu (mean of base Normal)
w_values = np.linspace(0.0, 1.0, 21)            # Mix between alpha and beta channels
w_fixed = 0.5
mu_values = np.linspace(-0.1, 0.1, 21)          # Change mean interaction values

# ======================
# ODE: community dynamics
# ======================

def derivative(t, N, r_alpha, r_beta, a_alpha_pos, a_alpha_neg, a_beta_pos, a_beta_neg, alpha, beta):
    """Compute dN/dt for the pos/neg-split channels.

    We compute Σ_pos = Σ_j A_pos[i,j] N_j and Σ_neg = Σ_j A_neg[i,j] N_j, then use
    (1 + Σ_pos)^exponent and (Σ_neg)^exponent for each channel, with clipping to ε>0.
    """
    N = np.array(N, dtype=float)

    # Interaction sums
    interaction_alpha_pos = np.dot(a_alpha_pos, N)
    interaction_alpha_neg = np.dot(a_alpha_neg, N)
    interaction_beta_pos  = np.dot(a_beta_pos,  N)
    interaction_beta_neg  = np.dot(a_beta_neg,  N)

    # Guard powers / gradients
    eps = 1e-12
    base_alpha_pos = 1.0 + np.clip(interaction_alpha_pos, eps, None)
    base_alpha_neg = np.clip(interaction_alpha_neg, eps, None)
    base_beta_pos  = 1.0 + np.clip(interaction_beta_pos,  eps, None)
    base_beta_neg  = np.clip(interaction_beta_neg,  eps, None)

    growth = ((r_alpha / alpha) * (np.power(base_alpha_pos, alpha) - np.power(base_alpha_neg, alpha)) +
              (r_beta  / beta)  * (np.power(base_beta_pos,  beta)  - np.power(base_beta_neg,  beta  )))

    return N * growth

# ==============================
# Random base interaction matrix
# ==============================

def generate_interaction_matrix(num_species, connectivity, mean, sd):
    """Generate a random interaction matrix with unit diagonal (overwritten downstream).

    Off-diagonals are realised with probability `connectivity`; values are Normal(mean, sd).
    The diagonal is set to 1.0 here, but downstream code overwrites it per channel policy.
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
    """Run solve_ivp in a thread and join with a wall-time limit.

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
                method='BDF',           # Stiff solver; alternatives: 'Radau', 'LSODA'
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

def compute_eigenvalues(N, r_alpha, r_beta, a_alpha_pos, a_alpha_neg, a_beta_pos, a_beta_neg, alpha, beta):
    """Approximate Jacobian eigenvalues at state N (pos/neg split, δ-term omitted).

    J[i,j] = N_i * [ r_alpha[i] ( A_alpha_pos[i,j] (1+Σ_alpha_pos)^{alpha-1} - A_alpha_neg[i,j] (Σ_alpha_neg)^{alpha-1} )
                   + r_beta[i]  ( A_beta_pos[i,j]  (1+Σ_beta_pos)^{beta-1}  - A_beta_neg[i,j]  (Σ_beta_neg)^{beta-1}  ) ].
    Consider clipping sums to ε if alpha < 1 or beta < 1 to avoid undefined powers.
    """
    S = len(N)
    if np.isscalar(r_alpha):
        r_alpha = np.full(S, r_alpha, dtype=float)
    if np.isscalar(r_beta):
        r_beta = np.full(S, r_beta, dtype=float)
    # Sum each terms
    sum_alpha_pos = np.dot(a_alpha_pos, N)
    sum_alpha_neg = np.dot(a_alpha_neg, N)
    sum_beta_pos  = np.dot(a_beta_pos,  N)
    sum_beta_neg  = np.dot(a_beta_neg,  N)
    # Safe clipping
    eps = 1e-12
    base_alpha_pos = 1.0 + np.clip(sum_alpha_pos, eps, None)
    base_alpha_neg = np.clip(sum_alpha_neg, eps, None)
    base_beta_pos  = 1.0 + np.clip(sum_beta_pos,  eps, None)
    base_beta_neg  = np.clip(sum_beta_neg,  eps, None)
    # Make the Jacobian matrix
    J = np.zeros((S, S), dtype=float)
    for i in range(S):
        for j in range(S):
            term_alpha_pos = r_alpha[i] * a_alpha_pos[i, j] * (base_alpha_pos[i] ** (alpha - 1))
            term_alpha_neg = r_alpha[i] * a_alpha_neg[i, j] * (base_alpha_neg[i] ** (alpha - 1))
            term_beta_pos  = r_beta[i]  * a_beta_pos[i,  j] * (base_beta_pos[i]  ** (beta  - 1))
            term_beta_neg  = r_beta[i]  * a_beta_neg[i,  j] * (base_beta_neg[i]  ** (beta  - 1))
            J[i, j] = N[i] * (term_alpha_pos - term_alpha_neg + term_beta_pos - term_beta_neg)

    eigenvalues = np.linalg.eigvals(J)
    return eigenvalues

# ========================================
# Single (alpha, S, w) simulation with retries
# ========================================

def simulate_combination(params):
    """Run one parameter triple with up to `max_attempts` tries to achieve zero extinctions.

    Returns a record: [alpha, N_species, w, leading_eigenvalue, extinction_counter, attempt_counter]
    leading_eigenvalue is the leading real part from the first full-survival success (or NaN if none).
    """
    # params = (alpha, N_species, third)
    alpha, N_species, third = params
    if Analysis_mode == 0:
        w = float(third)
        mu = a_alpha_mean   # keep the file's global mean fixed
    else:
        mu = float(third)
        w  = w_fixed        # keep w fixed when varying mu
    leading_eig = np.nan
    extinction_counter = 0

    for attempt in range(max_attempts):
        # Base interaction matrix (will be split by sign, then mixed across channels)
        a_raw = generate_interaction_matrix(N_species, connectivity=Connectivity, mean=mu, sd=a_alpha_std)

        # Split by sign
        neg_mask = (a_raw < 0.0)
        a_alpha_pos = np.zeros_like(a_raw)
        a_alpha_pos[neg_mask] = np.abs(a_raw[neg_mask])  # magnitudes of negatives
        a_alpha_neg = a_raw.copy()
        a_alpha_neg[neg_mask] = 0.0                      # only nonnegative entries remain

        # Mix into alpha and beta channels with weights
        a_beta_pos  = a_alpha_pos * (1.0 - w)
        a_beta_neg  = a_alpha_neg * (1.0 - w)
        a_alpha_pos = a_alpha_pos * w
        a_alpha_neg = a_alpha_neg * w

        # Channel-specific diagonals
        np.fill_diagonal(a_alpha_pos, 0.0)
        np.fill_diagonal(a_alpha_neg, 1.0)   # self-regulation via alpha_neg
        np.fill_diagonal(a_beta_pos,  0.0)
        np.fill_diagonal(a_beta_neg,  0.0)

        # Species-specific intrinsic rates (can be vectors; here constants)
        r_alpha = np.repeat(1.0,  N_species)
        r_beta  = np.repeat(0.25, N_species)

        # ---------- Stage 1: pre-perturbation integration ----------
        N_initial = np.repeat(N_ini, N_species)
        solver1 = ODESolverWithTimeout(
            derivative,
            [0.0, T_preturb],
            y0=N_initial,
            args=(r_alpha, r_beta, a_alpha_pos, a_alpha_neg, a_beta_pos, a_beta_neg, alpha, beta),
            max_time=max_time,
            precision=precision,
            max_step=0.1,
        )
        solved1 = solver1.solve()
        if (solved1 is None) or (not solved1.success):
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
            args=(r_alpha, r_beta, a_alpha_pos, a_alpha_neg, a_beta_pos, a_beta_neg, alpha, beta),
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
                                      a_alpha_pos, a_alpha_neg, a_beta_pos, a_beta_neg, alpha, beta)
        leading_val = np.max(np.real(eigvals))
        leading_eig = leading_val
        break

    return [alpha, N_species, w, leading_eig, extinction_counter, attempt + 1]

# ==============================
# Parallel sweep + save results
# ==============================
if __name__ == "__main__":
    if Analysis_mode == 0:
        param_combinations = list(product(alpha_values, N_species_values, w_values))
        third_name = "w"
    else:
        param_combinations = list(product(alpha_values, N_species_values, mu_values))
        third_name = "mu"

    results_list = []
    with ProcessPoolExecutor(max_workers=N_workers) as executor:
        futures = [executor.submit(simulate_combination, p) for p in param_combinations]
        for f in tqdm(as_completed(futures), total=len(futures), desc="Parameter combinations"):
            results_list.append(f.result())

    result_array = np.array(results_list, dtype=float)
    np.savetxt(
        "result_array.txt",
        result_array,
        delimiter=",",
        header=f"alpha,N_species,{third_name},leading_eigenvalue,extinction_counter,attempt_counter",
        comments="",
        fmt="%.6f",
    )
