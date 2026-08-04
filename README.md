# BayEZ: Bayesian Hierarchical EZ-DDM

Authors: 
- First: Anne Giacobello (Stan implementation)
- Second: Adriana F. Chávez De la Peña (JAGS counterpart)
- Senior authors: Julia Haaf and Joachim Vandekerckhove

This project implements a Bayesian hierarchical version of the EZ Drift Diffusion Model (EZ-DDM) with automatic design matrix support. The STAN implementation allows predictors to influence DDM parameters - drift rate ($\nu$), boundary separation ($\alpha$), and non-decision time ($\tau$) - without manually rewriting the model.

---

## Project Overview

**Main goals:**
1. Implement a flexible hierarchical EZ-DDM in STAN with design matrix automatization
2. Validate through simulation studies (parameter recovery)
3. Apply to real data (Grabowska et al., 2025 - Krakow dataset)
4. Develop a parallel JAGS implementation

**Paper with the data:** `Grabowska et al., 2025 Individual differences in neurophysiological correlates of post response adaptation.pdf`

---

## Folder Structure

### `/7 params ddm/`
This is an implementation of hierarchical models using the Seven-parameter Diffusion Model implementation for Stan by Henrich et al. (2024). We implemented several models and tested them on the Krakow dataset (Grabowska et al., 2025) using a tryout version (about 30 % of participants, 50% of trials, warmup = 500, sampling = 300) and a full version (all participants, all trials, warmup = 3000, sampling = 3000). 


| File | Description |
|------|-------------|
| `full_DDM_grabowska_comparison.R`| Runs three models in quicktest (super fast check whether the model has any misspecifications or errors), tryout and full mode: "zero", "free", or "sv_only". |
| `hierarchical_fullddm_no_intertrial_variability_nu_alpha.stan`| The "zero" model. The effect structure is the same for all models that we use with the Krakow data, for comparison. Meaning that nu (drift) has intercept fixed, condition/resp_type/interaction random, alpha (bound.) has intercept random (log scale), resp_type/condition fixed and tau (t0) has intercept random only. In the three models "zero", "free", or "sv_only", w (bias) is alway fixed at 0.5 and sw (variance of bias) is fixed at 0. This "zero" model also has sv (inter-trial variability in drift rate) and st0 (inter-trial variability in NDT) fixed to 0. |
| `hierarchical_fullddm_sv_only_nu_alpha.stan` | The "sv only" modoel. Again same effect structure. This model is similar to the zero but sv (inter-trial variability in drift rate) is a free, single population value and not fixed. |
| `hierarchical_fullddm_intertrial_variability_nu_alpha.stan` | The "free" model. bias and variance of bias are fixed, same effect structure as the other two models. But sv is a free, single population value and st0 is also a free, single population value. |
| `grabowska_comparison_more_models.R`| A new R script that can run all 10 models and can alternate between quicksafe, tryout and full model. | 
| `outputs` | This folder has all the estimates and figures that were saved by the script for every model. For an overview O only keep the output from the "full" models here. |




### `/working simulation + output/` (MAIN SIMULATION STUDY)
The **primary, working simulation study** with design matrix automatization.

| File | Description |
|------|-------------|
| `model2_v2_ncp.stan` | **Main STAN model** - Non-centered parameterization with: $\nu$ = condition + random slope, $\alpha$ = covariate, $\tau$ = intercept only. Uses design matrices $X_\nu$, $X_\alpha$, $X_\tau$ for flexible predictor specification. |
| `sim_and_rec_model2_ncp.R` | **Main simulation script** - Generates data, fits model, checks parameter recovery. Runs grid of $I \in \{80,160\}$ participants × $J \in \{80,160\}$ trials. Includes EZ point estimates for prior derivation, convergence checks, and recovery plots. |
| `output_model2_v2_ncp_new_params/` | Output folder with simulation results, RDS files, and recovery plots |

**Key design features:**
- Design matrices ($X_\nu[I,K,P]$, $X_\alpha[I,K,P]$, $X_\tau[I,K,P]$) allow automatic model adjustment based on predictors
- Non-centered parameterization (NCP) for random slopes: $v_i = \beta_{\nu,2} + \sigma_v \cdot z_{v,i}$
- Priors derived from pilot EZ point estimates

---

### `real data scripts and models/` (REAL DATA APPLICATION)
Application to the Krakow dataset (Grabowska et al., 2025).

| File | Description |
|------|-------------|
| `grabowska_ez_model_new.stan` | STAN model for real data - 4 cells (congruent/incongruent × post-correct/post-error), random slopes for condition, response type, and interaction effects on $\nu$; response type slope on $\alpha$ |
| `grabowska_ez_fit_krakow_new.R` | Main fitting script - Loads Krakow data, sets up design matrices with contrast coding (±1), runs STAN, produces forest plots and caterpillar plots |
| `grabowska_ez_model_ppn_intercept.stan` | Variant with person-specific random intercepts |
| `grabowska_ez_fit_krakow_randint.R` | Fitting script for random intercept model |
| `grabowska_ez_fit_krakow.R` | Earlier version of fitting script |
| `grabowska_ez_model.stan` | Earlier version of STAN model |
| `krakow_rt_checks.R` | RT descriptives and quality checks |
| `krakow_stability_check.R` | Multiple-run stability analysis |
| `fit_krakow_model_with_intercepts.R` | Additional intercept model fitting |

**4-cell design:**
- $k=1$: congruent, post-correct
- $k=2$: congruent, post-error  
- $k=3$: incongruent, post-correct
- $k=4$: incongruent, post-error

---

### `real data outputs/`
Results from real data analyses.

| File/Folder | Description |
|-------------|-------------|
| `output_krakow_new/` | Main output folder with `krakow_fit.rds`, `krakow_estimates.csv`, forest plots, caterpillar plots, cell-level violin plots |
| `output_krakow/` | Earlier analysis output |
| `stability_*.csv` | Results from stability check runs |
| `hist_rt_*.png` | RT distribution histograms |
| `rt_first_vs_rest*.png` | First trials vs rest comparison plots |

---

### `other sim and rec tryouts/`
Archive of model variants tested during development.

| File | Description |
|------|-------------|
| `model1_*.stan` | Models with condition effect on $\alpha$, covariate on $\nu$ |
| `model2_*.stan` | Models with condition effect on $\alpha$ only |
| `model3_*.stan` | Condition on $\alpha$, covariate on $\nu$, no random slope |
| `model4_*.stan` | Condition on both $\alpha$ and $\nu$, no covariate, no slope |
| `model5_*.stan` | Condition on $\alpha$, condition+covariate on $\nu$, no slope |
| `model6_*.stan` | Condition on both, slope on $\nu$ |
| `simulation_*.R` | Corresponding simulation scripts |
| `m*_recovery_*.png` | Parameter recovery plots from each model variant |
| `*.Rmd` | Documentation notebooks for various tryouts |

---

### `first tryouts Stan/`
Initial STAN model development (archived?).

| File | Description |
|------|-------------|
| `model.stan`, `old_model.stan`, `new_model.stan` | Early model versions |
| `*_w_cond_matrices*.stan` | Models experimenting with condition design matrices |
| `*_random_slope*.stan` | Models adding random slopes |
| `final_model_w_cond.stan` | Consolidated version before moving to main folder |
| `new_ez_ddm_stan_pipeline2.Rmd` | Early pipeline documentation |

---

### `JAGS/` (JAGS COUNTERPART)

Parallel JAGS implementation of the hierarchical EZ-DDM. 

**Top-level folders:**

| Folder | Description |
|--------|-------------|
| `src/` | Shared custom R functions used across reports (data generation, summary statistics, $\hat{R}$, recovery plots, simulation helpers) |
| `demos/` | Self-contained demo reports and their caches/outputs: `simulation-study/` (parameter recovery) and `applications/` (real data) |

**Shared `src/` files:**

| File | Description |
|------|-------------|
| `generate_truePars.R` | Generate true individual DDM parameters from hierarchical designs |
| `generate_trialData.R` | Simulate trial-level DDM data (random-walk emulation) |
| `generate_sumStats.R` | Sample summary statistics directly from EZ sampling distributions |
| `calculate_sumStats.R` | Compute $C$, mean RT, and RT variance from trial data |
| `calculate_Rhat.R` | $\hat{R}$ convergence diagnostic |
| `plot_recovery.R` | Base-R recovery plots (individual and population) |
| `simulation_settings.R` | Helpers for flexible design-matrix simulation output extraction |

**Three main reports:**

#### 1. `intro.Rmd` (simplest hierarchical EZ-DDM)

Introductory tutorial with background on the DDM / EZ-DDM and two short simulations (fixed $P$ and $T$): trial-level data generation vs direct summary-statistic sampling.

| Path | Role |
|------|------|
| `JAGS/intro.Rmd` | Main report (+ knitted HTML alongside it) |
| `JAGS/basic_model.bug` | Simple hierarchical EZ-DDM JAGS model (written/used by the report) |
| `JAGS/demos/simulation-study/simplest_model/` | Cached `.RData` simulation inputs/results |

#### 2. `sample_simStudy.Rmd` (design-matrix simulation study)

JAGS replication of Stan simulation (`working simulation + output/`). Uses design arrays and `inprod()` so the same model can fit different $P$ / $X$ structures. Includes two recovery studies.

| Path | Role |
|------|------|
| `JAGS/demos/simulation-study/design_matrix/sample_simStudy.Rmd` | Main report (+ knitted HTML) |
| `JAGS/demos/simulation-study/design_matrix/figures/` | Saved recovery plots (Study 1) |
| `JAGS/demos/simulation-study/design_matrix/figures_tauCond/` | Saved recovery plots (Study 2) |
| `JAGS/demos/simulation-study/design_matrix/*.RData` | Cached simulation results |

#### 3. `demo_realData.Rmd` (Krakow real-data application)

JAGS replication of Anne's Grabowska et al. (2025) Krakow analysis: 4 cells (congruence $\times$ previous accuracy), contrast coding, multiple random slopes. Behavioral predictors only (no EEG).

| Path | Role |
|------|------|
| `JAGS/demos/applications/demo_realData.Rmd` | Main report (+ knitted HTML) |
| `JAGS/demos/applications/model_krakow.bug` | Krakow JAGS model (written by the report via `writeLines()`) |
| `JAGS/demos/applications/data/krakow_data_standardized.csv` | Trial-level Krakow data ([OSF](https://osf.io/36pvg/overview)) |
| `JAGS/demos/applications/output_krakow/` | Cached JAGS fit / outputs |

---

### `paper structure/`
Paper draft materials.

| File | Description |
|------|-------------|
| `Bayes EZ DDM second draft.Rmd` | Main paper draft |
| `structure idea paper.Rmd` | Outline/structure notes |

---

### Root Files

| File | Description |
|------|-------------|
| `Outcome notes.Rmd` | **Important** - Notes comparing model results to Grabowska paper, discusses discrepancies (missing EEG predictor), stability checks, and lessons learned |
| `Grabowska et al., 2025...pdf` | Reference paper for real data application |

---

## Key Concepts

### EZ-DDM Parameters

| Parameter | Symbol | STAN code | JAGS code | Description |
|-----------|--------|-----------|-----------|-------------|
| Drift rate | $\nu$ | `nu` | `drift` | Evidence accumulation rate; higher = faster/more accurate decisions |
| Boundary separation | $\alpha$ | `alpha` | `bound` | Response caution; higher = slower but more accurate |
| Non-decision time | $\tau$ | `tau` | `nondt` | Motor/perceptual time not part of decision |

### EZ Forward Equations

$$
e = \exp(-\alpha \cdot \nu)
$$

$$
P_c = \frac{1}{1 + e}
$$

$$
\mu_{RT} = \tau + \frac{\alpha}{2\nu} \cdot \frac{1-e}{1+e}
$$

$$
\sigma^2_{RT} = \frac{\alpha}{2\nu^3} \cdot \frac{1 - 2\alpha\nu e - e^2}{(1+e)^2}
$$

### Design Matrix Approach
The model uses design matrices to flexibly specify predictors:
- $X_\nu[i,k,p]$: Predictors for drift rate (intercept, condition, covariates)
- $X_\alpha[i,k,p]$: Predictors for boundary (intercept, response type, etc.)
- $X_\tau[i,k,p]$: Predictors for non-decision time (typically just intercept)

This allows changing the model structure by modifying the design matrices rather than rewriting the model code.
