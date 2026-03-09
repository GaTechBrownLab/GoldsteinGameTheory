# GoldsteinGameTheory

Host–pathogen coevolution simulator based on the Goldstein et al. framework for punctuated evolutionary dynamics.

This model investigates how the evolution of conditional response strategies (evolved rules, ER) versus fixed trait values (evolved traits, ET) shapes the long-term dynamics of host–pathogen coevolution. The simulator implements a continuous-time Markov jump process under a strong selection, weak mutation regime, comparing four experimental conditions: **ET–ET**, **ET–ER**, **ER–ET**, and **ER–ER**.

This code is based on the original Java code written by J. Goldstein, and the original paper draft.

## Repository Structure

```
GoldsteinGameTheory/
├── simulation.py          # Core simulation engine
├── run_experiments.py     # CLI experiment runner
├── Plots.R                # Visualization and analysis (all figures)
├── figures/               # Generated figures (PDF + PNG)
├── results/               # Simulation output (generated)
│   └── {fitness_model}/
│       └── {condition}[_tags]/
│           ├── simulation.csv
│           └── config.json
└── README.md
```


## Files

### `simulation.py` — Core Engine

Contains all model architecture, classes, and functions. Designed as an importable module used by `run_experiments.py`.

- **Fitness models:** Acute infection (`fH = s/(s+d)`, `fP = v^β/(s+d)`), chronic infection (immunity-modulated virulence), minimal polynomial (`wH = c(1−c)(1−v)`, `wP = v(1−v)(1−c)`), and Taylor et al. 2006 (`H = [c/(v+c)]·[b/(m₀+c)]`, `P = vⁿ/(v+c)`)
- **`Simulation` class:** Manages evolutionary state, mutation proposals (equal-probability Gaussian quantile bins), selection via Kimura fixation probabilities (haploid or diploid), and the continuous-time substitution clock
- **Gillespie dynamics:** Each generation evaluates all mutations for both players, computes cumulative substitution rates weighted by mutation asymmetry (γ), chooses who mutates proportional to total rate, and draws exponential dwell time
- **Equilibrium solving:** Computes interior and boundary intersections of linear response strategies for ER scenarios, with 2-cycle fallback for unstable cases
- **`set_fitness_model()`:** Switches active fitness functions and adjusts trait domain (e.g., `[0,1]` for acute/minimal, `[0.001, 20]` for Taylor)

**Notation mapping (paper → code):**

| Paper | Code | Description |
|-------|------|-------------|
| ET (Evolved Trait) | EI / `evolved_strategy=False` | Fixed traits |
| ER (Evolved Response) | ES / `evolved_strategy=True` | Linear reaction norms |
| *c* (clearance) | `s` | Host immune clearance |
| *v* (virulence) | `v` | Pathogen virulence |
| *m_c*, *m_v* | `mS`, `mV` | Response slopes |
| *c₀*, *v₀* | `bS`, `bV` | Response intercepts |
| *W_H*, *W_P* | `hostFit`, `pathFit` | Host/pathogen fitness |
| γ | `prob_host_mutate` | Mutation rate asymmetry |

### `run_experiments.py` — Experiment Runner

CLI tool that configures and launches simulation runs across conditions. Sets `simulation.py` globals via `importlib` and manages output directory structure.

- **Condition control:** Maps `EThost_ETpath`, `EThost_ERpath`, `ERhost_ETpath`, `ERhost_ERpath` to the appropriate reactivity toggles
- **Fitness model selection:** `--fitness {acute, chronic, minimal, taylor}`
- **Runtime modes:** `--quick` (10K generations), custom `--gens N`, or default full (1M generations)
- **Mutation asymmetry:** `--gamma 0.01` (pathogen-fast, default), `--gamma 0.5` (equal), `--gamma 0.99` (host-fast)
- **Gamma sweeps:** `--gamma-sweep 0.01,0.1,0.5,0.9,0.99` runs all specified γ values across conditions
- **Nash dwell sweeps:** `--nash-sweep` runs parameter sweeps (σ, N, nS, nV) for weak Nash dwell analysis
- **Step-size sweeps:** `--step-sizes 0.01,0.05,0.2,0.5`
- **Diploid Kimura:** `--diploid` switches fixation probability to 4Ns denominator
- **Trait pinning:** `--fix-host 0.5` or `--fix-path 0.3` locks one player's trait
- **Population size:** `--pop-size N` overrides HOST_POP_N (PATH_POP_N scales to 100×)
- **Cost parameters:** `--nS` and `--nV` override clearance and virulence costs
- **Replicates:** `--rep N` offsets seed by N−1, tags output directory with `repN`
- **Tags:** `--tag name` appends to output directory (e.g., `--tag final_r1` for replicate tracking)
- **Mutation bins:** `--bins 11` reduces from default 51 for faster test runs
- **Reproducibility:** Saves full configuration to `config.json` per run

**Output directory naming:** `{condition}[_sigma{σ}][_diploid][_gamma{γ}][_fixH{v}][_fixP{v}][_N{n}][_rep{n}][_tag]/`

### `Plots.R` — Visualization

R script for generating all publication figures. Requires `tidyverse`, `ggplot2`, `patchwork`, `pracma`, `scales`, and `zoo`.

**Data loading:**
- `discover_experiments()` — Auto-scans `results/` for all `config.json` files, builds a cached catalog
- `load_sim()` — Loads a single experiment by model + scenario with automatic best-match selection
- `load_all_conditions()` — Loads all 4 conditions for one model into a single data frame
- `load_replicates()` — Loads replicate runs via `tag_prefix` matching (e.g., `"final_r"` finds `final_r1`, `final_r2`, `final_r3`)

**Figure functions:**

| Section | Function | Description |
|---------|----------|-------------|
| §4 | `fig_landscape` | Fitness landscape contour plots with best-response curves |
| §5 | `fig_strategy_panels` | Strategy line configurations at different slopes |
| §6 | `fig_timeseries` | 6×N grid: v, c, W_P, W_H, ω_P, ω_H across conditions |
| §7 | `fig_hex_combined` | Phase-space hexbin density on fitness landscapes |
| §8 | `fig_snapshots` | Strategy-line snapshots at different generations |
| §9 | `fig_nash_combined` | Nash violation map + slope stability distribution |
| §10 | `fig_strategy_evolution` | ER parameter (bS, mS, bV, mV) time series |
| §11 | `fig_ts_stats` | CV, spectral slope, correlation length distributions |
| §11 | `fig_acf_decay` | Autocorrelation function decay |
| §12 | `fig_step_sizes` | Mutational step-size distributions |
| §13 | `fig_neutral_drift` | Neutral mutation fraction analysis |
| §14 | `fig_dwell_times` | Dwell-time distributions in Nash/stable/boundary regions |
| §14 | `fig_boundary_occupancy` | Boundary occupancy fractions |
| §14 | `fig_trait_density` | 2D trait distributions with Nash overlay |
| §11 | `fig_gamma_*` | Gamma sweep analysis (timeseries, summary, symmetry, tempo, ACF) |
| — | `fig_replicate_timeseries` | Replicate overlay (v, c time series) |
| — | `fig_replicate_density` | Replicate trait distribution comparison |

**Replicate support:** Most figure functions accept `tag_prefix` to automatically load and overlay/pool replicates:
```r
# Single run
fig_timeseries("minimal", "ts_single", sigma = 0.01, diploid = TRUE,
               tag_filter = "final_r1")

# Overlay 3 replicates as colored semi-transparent lines
fig_timeseries("minimal", "ts_reps", sigma = 0.01, diploid = TRUE,
               tag_prefix = "final_r")

# Pool replicates for density estimation
fig_trait_density("minimal", sigma = 0.01, diploid = TRUE,
                  tag_prefix = "final_r")
```


## Quick Start

### 1. Run simulations (Python 3.8+, scipy optional)

```bash
# Quick test — all 4 conditions, acute model (~5–10 min)
python run_experiments.py --quick

# Full run — all 4 conditions, acute model
python run_experiments.py

# Single condition with a specific fitness model
python run_experiments.py --fitness minimal --condition ERhost_ERpath

# Diploid with small step size
python run_experiments.py -f minimal --diploid --step-sizes 0.01

# Gamma sweep across conditions
python run_experiments.py -f acute --gamma-sweep 0.01,0.1,0.5,0.9,0.99 --diploid

# Replicates with tagged output
python run_experiments.py -f minimal --diploid --step-sizes 0.01 --tag final_r1
python run_experiments.py -f minimal --diploid --step-sizes 0.01 --tag final_r2 --seed 3248233
python run_experiments.py -f minimal --diploid --step-sizes 0.01 --tag final_r3 --seed 3248234

# List all options
python run_experiments.py --list
```

### 2. Generate figures (R)

```r
source("Plots.R")
```

Reads simulation CSVs from `results/` and assembles multi-panel figures using `patchwork`. Call `refresh_catalog()` after running new simulations. Use `list_experiments()` to see all discovered runs.


## Experimental Conditions

| Condition | Host | Pathogen | Expected Dynamics |
|-----------|------|----------|-------------------|
| **EThost_ETpath** | Fixed trait *c* | Fixed trait *v* | Stable convergence to Nash equilibrium |
| **EThost_ERpath** | Fixed trait *c* | Linear rule *v(c) = v₀ + m_v · c* | Intermediate volatility |
| **ERhost_ETpath** | Linear rule *c(v) = c₀ + m_c · v* | Fixed trait *v* | Intermediate volatility |
| **ERhost_ERpath** | Linear rule *c(v)* | Linear rule *v(c)* | Punctuated equilibrium dynamics |


## Fitness Models

| Model | Host fitness | Pathogen fitness | Trait domain |
|-------|-------------|-----------------|--------------|
| **acute** | *s / (s + d)* | *v^β / (s + d)* | [0, 1] |
| **chronic** | *1 / d*, immunity modulates virulence | *(1−s) · v^β / d* | [0, 1] |
| **minimal** | *c(1−c)(1−v)* | *v(1−v)(1−c)* | [0, 1] |
| **taylor** | *[c/(v+c)] · [b/(m₀+c)]* | *vⁿ / (v+c)* | [0.001, 20] |


## Key Parameters

| Parameter | Default | CLI flag | Description |
|-----------|---------|----------|-------------|
| `std_dev_move` (σ) | 0.1 | `--step-sizes` | Mutation step size |
| `prob_host_mutate` (γ) | 0.01 | `--gamma` | Host mutation probability per event (~1:99 ratio) |
| `HOST_POP_N` | 10⁴ | `--pop-size` | Effective host population size |
| `PATH_POP_N` | 10⁶ | (100× host) | Effective pathogen population size |
| `nS` / `n_c` | 0.1 | `--nS` | Host cost-of-clearance scaling |
| `nV` / `n_v` | 1.0 | `--nV` | Pathogen cost-of-virulence scaling |
| `d₀` | 0.1 | — | Baseline mortality |
| `β` | 1.0 | — | Transmission–virulence exponent |
| `max_gens` | 10⁶ | `--gens` | Total generations after burn-in |
| `DIPLOID_KIMURA` | False | `--diploid` | Use 4Ns (diploid) vs 2Ns (haploid) fixation |
| `num_step_bins` | 51 | `--bins` | Mutation quantile bins (51² = 2601 ER candidates/player) |


## Output Files

Each run produces a directory at `results/{fitness_model}/{condition}[_tags]/` containing:

- **`simulation.csv`** — Per-substitution trajectory: evolutionary time, traits (*v*, *s*), strategy parameters (*bS*, *mS*, *bV*, *mV*), fitness values, omega (ω = substitution rate), mutator identity, mutation selection coefficient, Nash equilibrium status, and dwell times. Each recorded generation has a `pre` and `post` row.

- **`config.json`** — Full run configuration: condition, fitness model, reactivity toggles, population sizes, mutation parameters, random seed, timestamp, and all overrides.


## Performance Notes

- **ET/ET** runs are fast (~seconds per 1M generations): 51 mutations/player, direct fitness evaluation
- **ER conditions** are slower: 51×51 = 2,601 candidates per ER player, each requiring equilibrium solving
- **Small σ** (e.g., 0.01) with ER can take 1–3 days per condition at 1M generations
- Use `--bins 11` for ~20× faster test runs (11² = 121 ER candidates vs 2,601)
- Use `--gens 10000 --tag test` for quick parameter checks without overwriting production data


## Citation

Goldstein et al. "Punctuated instabilities of the immune system" *(in preparation)*

## Contact

[GaTechBrownLab](https://github.com/GaTechBrownLab/GoldsteinGameTheory)

Canan Karakoc (GitHub: @canankarakoc, canankarakoc@gmail.com)
