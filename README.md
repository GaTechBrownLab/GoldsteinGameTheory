# GoldsteinGameTheory

Host–pathogen coevolution simulator based on the Goldstein et al. framework for punctuated evolutionary dynamics.

This model investigates how the evolution of conditional response strategies (evolved rules, ER) versus fixed trait values (evolved traits, ET) shapes the long-term dynamics of host–pathogen coevolution. The simulator implements a continuous-time Markov jump process under a strong selection, weak mutation regime, comparing four experimental conditions: **ET–ET**, **ET–ER**, **ER–ET**, and **ER–ER**.

This code is based on the original java code written by J. Goldstein, and the original paper draft. 

## Repository Structure

```
GoldsteinGameTheory/
├── simulation.py          # Core simulation engine (classes, fitness models, equilibrium solving)
├── run_experiments.py     # CLI experiment runner (conditions, sweeps, configuration)
├── Plots.R                # Visualization and statistical analysis (all figures) saved in "figures" folder
├── results/               # Simulation output (generated)
│   └── {fitness_model}/
│       └── {condition}/
│           ├── simulation.csv
│           └── config.json
└── README.md
```


## Files

### `simulation.py` — Core Engine

Contains all model architecture, classes, and functions. Designed as an importable module used by `run_experiments.py`.

- **Fitness models:** Acute infection (`fH = s/(s+d)`, `fP = v^β/(s+d)`), chronic infection (immunity-modulated virulence), minimal polynomial (`wH = c(1−c)(1−v)`, `wP = v(1−v)(1−c)`), and Taylor et al. 2006 (`H = [c/(v+c)]·[b/(m₀+c)]`, `P = vⁿ/(v+c)`)
- **`Simulation` class:** Manages evolutionary state, mutation proposals (equal-probability Gaussian quantile bins), selection via Kimura fixation probabilities (haploid or diploid), and the continuous-time substitution clock
- **Equilibrium solving:** Computes interior and boundary intersections of linear response strategies for ER scenarios
- **Nash analysis:** Weak Nash classification for ER, strong Nash for ET
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

### `run_experiments.py` — Experiment Runner

CLI tool that configures and launches simulation runs across conditions. Sets `simulation.py` globals via `importlib` and manages output directory structure.

- **Condition control:** Maps `EThost_ETpath`, `EThost_ERpath`, `ERhost_ETpath`, `ERhost_ERpath` to the appropriate reactivity toggles
- **Fitness model selection:** `--fitness {acute, chronic, minimal, taylor}`
- **Runtime modes:** `--quick` (10K generations for testing) or default full (1M generations)
- **Step-size sweeps:** `--step-sizes 0.01,0.05,0.2,0.5` runs the same condition across multiple mutation step sizes
- **Diploid Kimura:** `--diploid` switches fixation probability to 4Ns denominator
- **Trait pinning:** `--fix-host 0.5` or `--fix-path 0.3` locks one player's trait at a constant value
- **Reproducibility:** Saves full configuration to `config.json` per run

### `Plots.R` — Visualization

R script for generating all publication figures. Requires `tidyverse`, `ggplot2`, `patchwork`, `pracma`, `scales`, and `zoo`.

### STRUCTURE:
-   §0  Setup (libraries, paths, theme)
-   §1  Fitness functions (all 4 models, defined ONCE)
-   §2  Utility functions (helpers used across figures)
-   §3  Data loading (unified loader for all experiments)
-   §4  FIGURE 1  — Fitness landscapes (any model)
-   §5  FIGURE 2  — Strategy lines & Nash equilibria
-   §6  FIGURE 3  — Time series (v, c, W_H, W_P, omega)
-   §7  FIGURE 4  — Phase-space density (hex plots on landscapes)
-   §8  FIGURE 5  — Snapshots of strategy lines (fig_snapshots)
-   §9  FIGURE 6  — Nash region & slope distribution (fig_nash_combined)
-   §10 FIGURE 7  — Strategy parameter evolution & stability
-   §11 Time series statistics (CV, spectral slope, ACF, memory)
-   §12 Step size analysis
-   §13 Neutral drift analysis
-   §14 Boundary analysis
-   §15 Cross-condition comparison figures (ts stats, steps, drift, dwell)


### 1. Run simulations (Python 3.8+, no external dependencies)

```bash
# Quick test — all 4 conditions, acute model (~5–10 min)
python run_experiments.py --quick

# Full run — all 4 conditions, acute model (~3–6 hours per condition)
python run_experiments.py

# Single condition with a specific fitness model
python run_experiments.py --fitness minimal --condition ERhost_ERpath

# Step-size sweep
python run_experiments.py --fitness acute --condition EThost_ETpath --step-sizes 0.01,0.05,0.2,0.5

# List all options
python run_experiments.py --list
```

### 2. Generate figures (R)

```r
source("Plots.R")
```

Reads simulation CSVs from `results/` and assembles multi-panel figures using `patchwork`.


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
| **chronic** | *s / (s + d)*, immunity modulates virulence | *(1−s)(v₀+v) / (s+d)* | [0, 1] |
| **minimal** | *c(1−c)(1−v)* | *v(1−v)(1−c)* | [0, 1] |
| **taylor** | *[c/(v+c)] · [b/(m₀+c)]* | *vⁿ / (v+c)* | [0.001, 20] |


## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nS` / `n_c` | 0.1 | Host cost-of-clearance scaling |
| `nV` / `n_v` | 1.0 | Pathogen cost-of-virulence scaling |
| `d₀` | 0.1 | Baseline mortality |
| `β` | 1.0 | Transmission–virulence exponent |
| `prob_host_mutate` | 0.01 | Host mutation probability per event (~1:99 ratio) |
| `HOST_POP_N` | 10⁴ | Effective host population size |
| `PATH_POP_N` | 10⁶ | Effective pathogen population size |
| `std_dev_move` | 0.1 | Mutation step size (σ) |
| `max_gens` | 10⁶ | Total generations after burn-in |


## Output Files

Each run produces a directory at `results/{fitness_model}/{condition}/` containing:

- **`simulation.csv`** — Per-substitution trajectory: evolutionary time, traits (*v*, *s*), strategy parameters (*bS*, *mS*, *bV*, *mV*), fitness values, omega (ω = substitution rate / neutral rate), mutator identity, mutation classification, Nash equilibrium status, and dwell times. Each generation records a `pre` and `post` row.

- **`config.json`** — Full run configuration: condition, fitness model, reactivity toggles, population sizes, mutation parameters, random seed, and timestamp.


## Citation

Goldstein et al. "Punctuated instabilities of the immune system" *(in preparation)*

## Contact

[GaTechBrownLab](https://github.com/GaTechBrownLab/GoldsteinGameTheory)

Canan Karakoç (GitHub: @canankarakoc, canankarakoc@gmail.com)