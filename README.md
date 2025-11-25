# Species diversity typically reduces ecological stability  
This repository contains code and processed data used to generate the main and supplementary figures for the manuscript titled above. Each top-level folder corresponds to a figure (or supplementary figure) in the manuscript.

## Requirements

Analyses were run with:

- **R** (v4.x recommended) with typical tidy/Bayesian ecosystem packages, including:
  - `data.table`, `dplyr`, `tidyr`, `ggplot2`
  - `brms`, `future`, `future.apply`
- **Python 3** with:
  - `numpy`, `scipy`, `pandas`, `matplotlib`
- **Jupyter** (optional) for running the `.ipynb` notebooks.

Exact versions are not critical; scripts were developed with recent versions of R and Python.

## Folder structure

### Fig. 3 | Theory Heat Map

`Fig3_TheoryHeatMap/` contains numerical simulations and plotting code for the theoretical heat maps.

Contents:

- `coth_parl.py` – Python script that computes stability metrics over a grid of parameter values and writes results to `result_array.txt`.
- `result_array.txt` – Pre-computed simulation output used to make the heat maps.
- `MakeHmap.R` – R script that reads `result_array.txt` and generates the figure panels.

---

### Fig. 4 | Hatton Re-analysis

`Fig4_HattonReanalysis/` is the re-analysis of the Hatton et al. macroecological dataset and generation of the corresponding figure.

Contents:

- `metadata_Hatton_v1.2.csv` – Metadata for datasets / taxa.
- `input/fig4.csv` – Processed data file used for model fitting.
- `output/Hatton_v1.2.rds` – Saved `brms` fit (posterior draws).
- `output/mean_Hatton_v1.2.csv` – Posterior summaries used for plotting.
- `codes/brms_Hatton_v1.2.R` – Main R script to fit the θ-logistic model with `brms`.
- `codes/Make_Mean_summary_Hatton.R` – Creates summary CSV from the posterior draws.
- `codes/MakeErrbar_Hatton.R` – Builds error-bar summaries for plotting.
- `codes/Plot_posterior_Hatton.R` – Plotting script for posterior curves and main-text figure.
- `codes/K_prior_justification.R` – Helper script for prior checks.
- `codes/match_plot_Hatton_CI_v2.R` – Script to align fitted credible intervals with the original Hatton presentation.

Set the working directory to `Fig4_HattonReanalysis/` and run the scripts in `codes/` in the above logical order to reproduce fits and plots.

---

### Fig. 5 | Microbial Fitting

`Fig5_MicroFitting/` contains Bayesian fits for microbial monoculture growth curves and associated figure.

Contents:

- `metadata_microbial_v1.2.csv` – Metadata describing microbial datasets.
- `input/processed_microbial.csv` – Pre-processed microbial growth data.
- `output/microbial_v1.2.rds` – Saved `brms` fit (posterior draws).
- `output/mean_microbial_v1.2.csv` – Posterior summaries used for plotting.
- `codes/brms_microbial_v1.2.R` – Main R script fitting θ-logistic models to the microbial data.
- `codes/Make_Mean_summary.R` – Creates summary CSV from the posterior draws.
- `codes/MakeErrbar_microbial.R` – Builds error-bar summaries for the figure.
- `codes/Plot_posterior_microbial.R` – Plotting script for microbial posterior curves.
- `codes/match_plot_microbial_CI_v1.R` – Script for matching fitted curves and credible intervals to the figure layout.

Run the scripts in `codes/` (after ensuring paths to `input/processed_microbial.csv` are correct) to reproduce the fits and figure.

---

### Fig. S4 | Noisy gLV

`FigS4_NoisyGLV/` provides supporting simulations for Supplementary Figure S4: effects of observational noise on fitted exponents in generalized Lotka–Volterra (gLV) models.

Contents:

- `stationary/NoisyAlpha.ipynb` – Jupyter notebook for simulations and fitting under stationary conditions.
- `stationary/fitted_alphas.csv` – Output summary of fitted exponents from the stationary simulations.
- `growth/NoisyAlpha.ipynb` – Jupyter notebook for simulations and fitting for growth-phase data.
- `growth/fitted_alphas.csv` – Output summary of fitted exponents from the growth simulations.
- `.ipynb_checkpoints/` – Auto-generated Jupyter checkpoint files (can be ignored).

Open the notebooks in Jupyter and run all cells to regenerate the `fitted_alphas.csv` files and associated plots.

---

### Fig. S11 | Neutral Community

`FigS11_NeutralComm/` contains the simulations for Supplementary Figure S11: effects of mean interspecific interactions and neutrality/competition on stability.

Contents:

- `coth_int.py` – Python script performing simulations for neutral and competitive communities and writing results to text files.
- `neutral/result_array.txt` – Simulation output for the neutral community scenario.
- `competitive/result_array.txt` – Simulation output for the competitive community scenario.

Use `coth_int.py` to regenerate the `result_array.txt` files if desired, then visualise them with your plotting code following the manuscript / SI description.

---

### Fig. S12 | Mean Interactions

`FigS12_MeanInteractions/` provides additional simulations for Supplementary Figure S12, exploring mean interspecific interactions.

Contents:

- `result_array.txt` – Pre-computed numerical results used to build Supplementary Figure S12.


## Notes

- Run scripts from within their respective figure folders (or set your working directory accordingly) so that relative file paths resolve correctly.
- Some R scripts refer to data paths such as `"data/..."`; in this archive the corresponding files are under `input/` and `output/`. Adjust working directories or paths as needed.
- For full methodological details (model structure, priors, and data processing), please refer to the main manuscript and Supplementary Information .