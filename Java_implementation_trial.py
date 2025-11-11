"""
Host–Pathogen Coevolution Simulator (EI vs ES) - ACUTE INFECTION MODEL
==========================================================================

PYTHON: 3.8+ (tested on 3.9, 3.10, 3.11)
DEPENDENCIES: None (standard library only)
  - math, random, csv, os, json: Built-in
  - functools.lru_cache: Python 3.2+
  - bisect: Built-in
  - typing: Python 3.5+

INSTALLATION:
  No installation needed - just run:
    python Java_implementation_trial.py

QUICK START:
  1. First run (test mode):
     - Set max_gens = 10_000 (line ~58)
     - Run: python Java_implementation_trial.py
     - Check: results/*.csv files created
     - Time: ~5-10 minutes
  
  2. Full run (publication mode):
     - Set max_gens = 1_000_000 (line ~58)
     - Run: python Java_implementation_trial.py
     - Time: ~3-6 hours
  
  3. Analyze in R:
     - Open Goldstein_model_notes_plots.Rmd
     - Results automatically loaded from results/*.csv

NEW FEATURES (v1.4):
  ==================
  
  1. REACTIVITY TOGGLES (ES mode only):
     - FIX_HOST_REACTIVITY = True/False (line ~202)
       * When True: Host slope (mS) fixed at 0, only intercept (bS) evolves
       * Host becomes non-reactive but can still evolve baseline clearance
       * Expected: Lower "noisy dynamics" than dual ES
     
     - FIX_PATH_REACTIVITY = True/False (line ~203)
       * When True: Pathogen slope (mV) fixed at 0, only intercept (bV) evolves  
       * Pathogen becomes non-reactive but can still evolve baseline virulence
       * Expected: Lower "noisy dynamics" than dual ES
     
     Usage example:
       # Host evolves reactively, pathogen is non-reactive
       FIX_HOST_REACTIVITY = False
       FIX_PATH_REACTIVITY = True
     
     Important: 
       - These only affect ES mode (evolved_strategy=True)
       - In EI mode, these are ignored (traits evolve directly)
       - Can fix both (mS=0, mV=0), one, or neither
  
  2. FITNESS MODEL TOGGLE:
     - FITNESS_MODEL = "acute" or "minimal" (line ~206)
       
     * "acute" (default): Current biological model
       - Host: fH = c/(c+d), where d = mortality with convex costs
       - Path: fP = v^β/(c+d)
       - Steep fitness landscape, strong selection
     
     * "minimal": Simplified multiplicative model
       - Host: wh = c(1-c)(1-v)
       - Path: wp = v(1-v)(1-c)
       - Gentler fitness landscape
     
     Usage:
       FITNESS_MODEL = "minimal"  # Switch to minimal model

REPRODUCIBILITY INFORMATION
---------------------------
Model: Goldstein et al. "Punctuated instabilities of the immune system"
Implementation: Python 3.8+ (tested on 3.9, 3.10)
Required libraries: None (uses only standard library)
Expected runtime: 
  - 10K gens: ~5-10 minutes
  - 100K gens: ~30-60 minutes  
  - 1M gens: ~3-6 hours (depends on CPU)
Memory usage: ~500MB-2GB (scales with max_gens and write frequency)

FILE STRUCTURE
--------------
Output files created:
  results/
    ei_simulation.csv       # EI mode full trajectory
    es_simulation.csv       # ES mode full trajectory
    run_config.json         # Parameter record for this run
  postproc/
    fig6_timeseries.csv     # Time series for plotting
    fig5_panel_*.csv        # Fitness landscape grids
    fig5_panel_*.meta.json  # Strategy line metadata
    fig7_region_acute.csv   # Nash equilibrium feasibility map
    fig7_occupancy.csv      # Dwell-weighted occupancy heatmap

SCIENTIFIC BACKGROUND
--------------------
This version uses the ACUTE INFECTION MODEL with biological parameters from the paper,
NOT the test Gaussian fitness functions that the Java reference uses.

Key Corrections from Java Implementation:
1. Fitness: Uses ACUTE model (s/(s+m) and v^β/(s+m)), not Gaussian test functions
2. Parameters: Biological values from paper (nS=0.1, nV=1.0, d0=0.1, β=1.0)
3. Equilibrium finding: Interior first (most common), then boundaries
4. Best intersection: Evaluates ALL intersections, picks best for mutator
5. EI factor: Applied to omega calculation ONLY, not to mutation rate selection
6. Omega: No "-1" subtraction (neutral evolution → ω≈1)
7. Gamma: 0.01 per paper (pathogens 100× faster than hosts)

Omega is a DIMENSIONLESS ratio = (substitution rate) / (neutral rate).

WHY EI HAS LOW OMEGA (~0.02-0.10):
  - EI evolves near Nash equilibrium where fitness is nearly flat
  - Most mutations are nearly neutral or slightly deleterious  
  - Only ~2-4% of mutations are advantageous (paper states this explicitly!)
  - This naturally gives omega ≈ 0.02-0.04 without any corrections

WHY ES HAS VARIABLE OMEGA (0.01-1000+):
  - During stasis: near Nash, omega low (~0.01-0.1) like EI
  - During bursts: many adaptive mutations, omega spikes (10-1000+)
  - Punctuated equilibrium pattern visible in omega time series

Biological Model (Acute Infection):
- Host fitness: f_H = c/(c+d) where d = d0 + cost(c) + cost(v)
- Path fitness: f_P = v^β/(c+d)
- Costs are convex (approach infinity as c→1 or v→1)
- c ∈ [0,1]: host immune clearance rate
- v ∈ [0,1]: pathogen virulence (affects transmission and host mortality)

EVOLUTIONARY MODES
------------------
EI (Evolved Interaction): 
  - Host and pathogen evolve fixed trait values (c̃, ṽ)
  - Each generation: mutate trait directly, evaluate fitness, accept/reject
  - Expected: Slow, stable convergence to Nash equilibrium
  
ES (Evolved Strategy):
  - Host and pathogen evolve response strategies: c(v) = b_S + m_S·v, v(c) = b_V + m_V·c
  - Each generation: mutate strategy parameters, find equilibrium, accept/reject
  - Expected: Punctuated equilibrium (bursts of rapid evolution, then stasis)

EXPECTED RESULTS (based on paper)
----------------
Good simulation should show:
  ✓ Traits explore interior (0.1 < v,s < 0.9 range for ES; narrow for EI)
  ✓ ES: Bursts visible in time series (omega spikes)
  ✓ EI: Smooth, gradual changes (flat-ish in plots)
  ✓ Path/Host mutation ratio ≈ 99:1 (from γ=0.01)
  ✓ ES interior equilibrium used 70-90% of time
  
OMEGA VALUES (key insight from paper):
  ✓ EI omegas: LOW (~0.02-0.10) - near Nash equilibrium, weak selection
  ✓ ES omegas: VARIABLE (0.01-1000+) - bursts during punctuated evolution
  ✓ ES baseline: Similar to EI (~0.01-0.1) during stasis
  ✓ ES bursts: Spikes to 10-1000+ during adaptive cascades
  
The paper states EI has "average ratio of substitution rate to neutral 
substitution rate is approximately 2-4%" - this is omega ≈ 0.02-0.04!
  
Specific expectations:
  - EI: Small v/s ranges (e.g., v: 0.50-0.55, s: 0.75-0.80)
  - EI: Few unique values (~20-100 for v, ~100-500 for s)
  - EI: omega_host ≈ 0.05-0.15, omega_path ≈ 0.02-0.10
  - ES: Full exploration (v: 0-1, s: 0-1)
  - ES: Many unique values (>10,000 for both)
  - ES: omega mean ≈ 1-10, but with huge spikes (10-1000+) during bursts

TROUBLESHOOTING
---------------
Issue: Traits stuck at 0 or 1
  → Check: cost parameters may be too small (nS, nV)
  → Try: Increase std_dev_move to 0.1
  
Issue: No evolution (zero rate streak)
  → Check: May have reached exact Nash equilibrium (rare but possible)
  → Try: Restart with different seed
  
Issue: ES takes forever
  → Check: Equilibrium finding may be slow on this trajectory
  → Try: Reduce CACHE_EQ size or disable (line ~222)
  
Issue: Path/Host ratio not 99:1
  → Check: gamma parameter (should be 0.01)
  → Check: factor multiplication in EI proposals

CITATION
--------
If using this model, please cite:
Goldstein et al. "Punctuated instabilities of the immune system"
[Add full citation when published]

CODE MAINTENANCE
----------------
Last updated: 2025
Contact: https://github.com/GaTechBrownLab/GoldsteinGameTheory
License: [Add license - MIT, GPL, etc.]

 NOTE: Java reference code also has a chronic infection model with
      'effectiveness' parameter for immune suppression. Not implemented
      here but could be added later if needed. See Simulate.java line ~30.
"""

import math
import random
import csv
import os
import sys
from bisect import bisect_left
from typing import Optional, Tuple, List
import json

__version__ = "1.4-toggles"
__date__ = "2025"
__model__ = "coevolution_with_reactivity_and_fitness_toggles"

# =========================
# FITNESS FUNCTION TOGGLE
# =========================
# CRITICAL: Set this flag to choose fitness model
#
# True:  Gaussian test functions
#        - Gentle bowl-shaped landscape centered at (0.5, 0.5)
#        - Weak selection (fitness varies 1.0 to 1.0001)
#        - Stable Nash equilibria form naturally
#        - Shows clear punctuated equilibrium
#        - Mean omega ~1-10, median ~0.1-1
#
# False: Biological acute infection model
#        - Steep fitness costs (nS, nV parameters)
#        - Strong selection throughout
#        - Unstable dynamics, no clear stasis
#        - Mean omega >100k (continuous evolution)

USE_GAUSSIAN = False  # ← CHANGE THIS: True=test, False=biology

# =========================
# Global Parameters
# =========================

# ==================== NEW TOGGLES (as of v1.4) ====================
# Reactivity controls (ES mode only)
FIX_HOST_REACTIVITY = False    # If True, force mS=0 (host non-reactive, but bS evolves)
FIX_PATH_REACTIVITY = False    # If True, force mV=0 (pathogen non-reactive, but bV evolves)

# Fitness model selection
FITNESS_MODEL = "acute"         # "acute" (current model) or "minimal" (wh=c(1-c)(1-v), wp=v(1-v)(1-c))
# ==================================================================

# Biological parameters (only used if USE_GAUSSIAN = False)
nS = 0.1   # Cost parameter for host clearance rate
nV = 1.0   # Cost parameter for pathogen clearance rate
d0 = 0.1   # Initial pathogen density
beta = 1.0 # Transmission rate
ONE_PLUS_EPSILON = 1.0 + 1.0/999.0

# Population sizes
Ne_H = 1.0e4 
Ne_P = 1.0e6

# Mutation parameters
host_vs_pathogen_mut_rate = 0.01  # γ = pathogen mutation rate / host mutation rate
num_step_bins = 51  # 51 per paper; increase to 101 for smoother but slower
std_dev_move = 0.1 # 0.01 per paper; increase to 0.1 if traits get stuck
std_dev_angle = 0.1 * math.pi # radians; increase to 0.1*pi if stuck

# Simulation parameters
burn_in_gens = 10_000    # Generations to skip before recording stats
max_gens = 1_000_000     # ← CHANGE THIS: 10K for test, 1M for full run
seed = 3248232           # ← CHANGE THIS: Any integer for RNG seed
write_every = 1          # Write every N generations (only in "full" mode)

# Runtime mode
RUNTIME_MODE = "full"   # "full" or "fast"
WRITE_EVERY_FAST = 100  # Write every N gens in "fast" mode
DWELL_MIN = 1e-12       # Minimum dwell time to record (to avoid log(0))

# Numerical safeties
ANGLE_EPS = 1e-4 # To avoid tan(±π/2)
FDH = 1e-6  # Finite difference step for numerical gradients
NEUTRAL_THRESH = 1e-8   # Selection coeff below which is neutral
MAX_FP_ITERS = 80       # Max fixed-point iterations for ES equilibrium
PROGRESS_EVERY = 10_000 
CACHE_EQ = True 

# Postproc sizes
GRID_N_FULL = 151;  BINS_FULL = 60 
GRID_N_FAST = 101;  BINS_FAST = 40 

def _runtime_write_every(gen: int) -> bool:         
    if RUNTIME_MODE == "full":
        return (gen >= 0) and (gen % write_every == 0)
    else:
        return (gen >= 0) and (gen % WRITE_EVERY_FAST == 0)

def _postproc_sizes():
    return (GRID_N_FULL, BINS_FULL) if RUNTIME_MODE == "full" else (GRID_N_FAST, BINS_FAST)

# =========================
# Utilities
# =========================

def clamp01(x: float) -> float:
    return 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)

def safe_tan(theta: float) -> float:
    half_pi = 0.5 * math.pi
    if theta >= half_pi - ANGLE_EPS:
        theta = half_pi - ANGLE_EPS
    elif theta <= -half_pi + ANGLE_EPS:
        theta = -half_pi + ANGLE_EPS
    return math.tan(theta)

def wrap_angle(theta: float) -> float:
    theta = (theta + math.pi) % (2.0 * math.pi) - math.pi
    if theta >= math.pi/2:
        theta -= math.pi
    elif theta <= -math.pi/2:
        theta += math.pi
    if abs(abs(theta) - math.pi/2) < ANGLE_EPS:
        theta = math.copysign(math.pi/2 - ANGLE_EPS, theta)
    return theta

def gaussian_step_bins(n: int) -> List[float]:
    from math import sqrt, log
    def inv_norm_cdf(p: float) -> float:
        a = [-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00]
        b = [-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01]
        c = [-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00]
        d = [7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00]
        plow = 0.02425; phigh = 1 - plow
        if p < plow:
            q = math.sqrt(-2*log(p))
            return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        if p > phigh:
            q = math.sqrt(-2*log(1-p))
            return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
        q = p - 0.5; r = q*q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    return [inv_norm_cdf((i+0.5)/n) for i in range(n)]

# =========================
# Fitness Functions - BOTH MODELS
# =========================

# -------- BIOLOGICAL ACUTE INFECTION MODEL --------
def compute_mortality(v: float, s: float) -> float:
    """
    Mortality rate with convex costs.
    Costs approach infinity as traits → 1 (vertical asymptotes).
    """
    s_term = (nS * ONE_PLUS_EPSILON * s) / (ONE_PLUS_EPSILON - s)
    v_term = (nV * ONE_PLUS_EPSILON * v) / (ONE_PLUS_EPSILON - v)
    return d0 + s_term + v_term

def host_fitness_biological(v: float, s: float) -> float:
    """
    Biological model: f_H = c/(c+d)
    Host fitness = probability of clearing infection.
    Creates steep fitness landscape with strong selection.
    """
    m = compute_mortality(v, s)
    return s / (s + m)

def path_fitness_biological(v: float, s: float) -> float:
    """
    Biological model: f_P = v^β/(c+d)
    Pathogen fitness = transmission rate × infection duration.
    Creates steep fitness landscape with strong selection.
    """
    m = compute_mortality(v, s) 
    tr = v**beta
    return tr / (s + m)


# -------- MINIMAL FITNESS MODEL --------
def host_fitness_minimal(v: float, s: float) -> float:
    """
    Minimal model: wh = c(1-c)(1-v)
    Simple multiplicative fitness without complex mortality.
    """
    return s * (1.0 - s) * (1.0 - v)

def path_fitness_minimal(v: float, s: float) -> float:
    """
    Minimal model: wp = v(1-v)(1-c)  
    Simple multiplicative fitness without complex mortality.
    """
    return v * (1.0 - v) * (1.0 - s)


# -------- GAUSSIAN TEST FUNCTIONS --------
def host_fitness_gaussian(v: float, s: float) -> float:
    """
    Gaussian test function (reproduces paper's results).
    
    Creates gentle bowl-shaped landscape centered at (0.5, 0.5).
    Fitness varies minimally (1.0 to 1.0001) → weak selection.
    This allows Nash equilibria to form and persist.
    """
    central = (v - 0.5)**2 + (s - 0.5)**2
    return 1.0 + 0.0001 * math.exp(-central / 0.5)

def path_fitness_gaussian(v: float, s: float) -> float:
    """
    Gaussian test function (reproduces paper's results).
    
    Even weaker selection than host (1.0 to 1.00001).
    Creates stable coevolutionary dynamics with clear punctuated equilibrium.
    """
    central = (v - 0.5)**2 + (s - 0.5)**2
    return 1.0 + 0.00001 * math.exp(-central / 0.5)


# -------- MASTER SWITCH --------
def host_fitness(v: float, s: float) -> float:
    """
    Master host fitness function.
    Delegates based on USE_GAUSSIAN and FITNESS_MODEL flags.
    """
    if USE_GAUSSIAN:
        return host_fitness_gaussian(v, s)
    elif FITNESS_MODEL == "minimal":
        return host_fitness_minimal(v, s)
    else:  # "acute"
        return host_fitness_biological(v, s)

def path_fitness(v: float, s: float) -> float:
    """
    Master pathogen fitness function.
    Delegates based on USE_GAUSSIAN and FITNESS_MODEL flags.
    """
    if USE_GAUSSIAN:
        return path_fitness_gaussian(v, s)
    elif FITNESS_MODEL == "minimal":
        return path_fitness_minimal(v, s)
    else:  # "acute"
        return path_fitness_biological(v, s)


# =========================
# Helper Functions
# =========================

def selection_coeff(new_fit: float, old_fit: float) -> float:
    denom = old_fit if old_fit > 1e-12 else 1e-12
    return (new_fit - old_fit) / denom

def kimura_intensity(s: float, Ne: float) -> float:
    if abs(Ne * s) < 1e-8:
        return 1.0
    if Ne * s < -25:
        return 0.0
    if Ne * s > 25:
        return 2.0 * Ne
    num = 1.0 - math.exp(-2.0 * s)
    den = 1.0 - math.exp(-2.0 * Ne * s)
    if abs(den) < 1e-16:
        return 2.0 * Ne if s > 0 else 0.0
    return 2.0 * Ne * num / den

def grad_host(v: float, s: float, h: float = FDH) -> Tuple[float,float]:
    v1, v2 = clamp01(v-h), clamp01(v+h)
    s1, s2 = clamp01(s-h), clamp01(s+h)
    dfdv = (host_fitness(v2, s) - host_fitness(v1, s)) / max((v2 - v1), 1e-12)
    dfds = (host_fitness(v, s2) - host_fitness(v, s1)) / max((s2 - s1), 1e-12)
    return dfdv, dfds

def grad_path(v: float, s: float, h: float = FDH) -> Tuple[float,float]:
    v1, v2 = clamp01(v-h), clamp01(v+h)
    s1, s2 = clamp01(s-h), clamp01(s+h)
    dfdv = (path_fitness(v2, s) - path_fitness(v1, s)) / max((v2 - v1), 1e-12)
    dfds = (path_fitness(v, s2) - path_fitness(v, s1)) / max((s2 - s1), 1e-12)
    return dfdv, dfds

# =========================
# ES Equilibrium Finding
# =========================

def apply_host_strategy(v: float, bS: float, mS: float) -> float:
    return clamp01(bS + mS * v)

def apply_path_strategy(s: float, bV: float, mV: float) -> float:
    return clamp01(bV + mV * s)

def find_all_intersections(bV: float, mV: float, bS: float, mS: float) -> Optional[List[Tuple[float, float]]]:
    """
    Find where host and pathogen strategy lines intersect.
    Returns list of (v, s) equilibrium points or None.
    """
    # Try interior analytic solution first (most common)
    denom = 1.0 - mS * mV
    if abs(denom) > 1e-12:
        s_int = (bS + mS * bV) / denom
        v_int = bV + mV * s_int
        if (0.0 < v_int < 1.0 and 0.0 < s_int < 1.0 and abs(mS * mV) < 1.0):
            return [(v_int, s_int)]
    
    # Check boundary cases
    candidates = []
    
    v_at_s0 = clamp01(bV + mV * 0.0)
    s_at_v = clamp01(bS + mS * v_at_s0)
    if abs(s_at_v) < 1e-6:
        candidates.append((v_at_s0, 0.0))
    
    v_at_s1 = clamp01(bV + mV * 1.0)
    s_at_v = clamp01(bS + mS * v_at_s1)
    if abs(s_at_v - 1.0) < 1e-6:
        candidates.append((v_at_s1, 1.0))
    
    s_at_v0 = clamp01(bS + mS * 0.0)
    v_at_s = clamp01(bV + mV * s_at_v0)
    if abs(v_at_s) < 1e-6:
        candidates.append((0.0, s_at_v0))
    
    s_at_v1 = clamp01(bS + mS * 1.0)
    v_at_s = clamp01(bV + mV * s_at_v1)
    if abs(v_at_s - 1.0) < 1e-6:
        candidates.append((1.0, s_at_v1))
    
    unique = []
    for cand in candidates:
        is_dup = any(abs(cand[0] - u[0]) + abs(cand[1] - u[1]) < 1e-6 for u in unique)
        if not is_dup:
            unique.append(cand)
    
    return unique if unique else None

# Caching wrapper
from functools import lru_cache

def _q(x: float, q: int = 6) -> float:
    return round(float(x), q)

@lru_cache(maxsize=200_000)
def _find_all_intersections_cached(bV: float, mV: float, bS: float, mS: float):
    return find_all_intersections(bV, mV, bS, mS)

def find_intersections_fast(bV, mV, bS, mS):
    if CACHE_EQ:
        return _find_all_intersections_cached(_q(bV), _q(mV), _q(bS), _q(mS))
    else:
        return find_all_intersections(bV, mV, bS, mS)

# =========================
# Nash Checks
# =========================

def es_weak_nash(v: float, s: float, bV: float, mV: float, bS: float, mS: float) -> Tuple[bool, bool, bool]:
    gh_v, gh_s = grad_host(v, s)
    gp_v, gp_s = grad_path(v, s)
    dot_host = gh_v * mV + gh_s * 1.0
    dot_path = gp_v * 1.0 + gp_s * mS
    tangency_ok = (abs(dot_host) < 1e-6) and (abs(dot_path) < 1e-6)

    ds = 1e-4
    s1, s2 = clamp01(s - ds), clamp01(s + ds)
    v1, v2 = clamp01(bV + mV * s1), clamp01(bV + mV * s2)
    Hc = host_fitness(v, s); H1 = host_fitness(v1, s1); H2 = host_fitness(v2, s2)
    host_max = (Hc >= H1 - 1e-10) and (Hc >= H2 - 1e-10)

    dv = 1e-4
    v1b, v2b = clamp01(v - dv), clamp01(v + dv)
    s1b, s2b = clamp01(bS + mS * v1b), clamp01(bS + mS * v2b)
    Pc = path_fitness(v, s); P1 = path_fitness(v1b, s1b); P2 = path_fitness(v2b, s2b)
    path_max = (Pc >= P1 - 1e-10) and (Pc >= P2 - 1e-10)

    return (tangency_ok and host_max and path_max), host_max, path_max

def ei_strong_nash(v: float, s: float) -> bool:
    ds = 1e-4; dv = 1e-4
    s1, s2 = clamp01(s - ds), clamp01(s + ds)
    v1, v2 = clamp01(v - dv), clamp01(v + dv)
    Hc = host_fitness(v, s)
    Pc = path_fitness(v, s)
    return (Hc >= host_fitness(v, s1) - 1e-10 and Hc >= host_fitness(v, s2) - 1e-10) and \
           (Pc >= path_fitness(v1, s) - 1e-10 and Pc >= path_fitness(v2, s) - 1e-10)

# =========================
# Helpers for Postproc
# =========================

ONE_PLUS_EPS_HLP = 1.0 + 1.0/999.0
nS_HLP, nV_HLP, d0_HLP, beta_HLP = 0.1, 1.0, 0.1, 1.0

def clamp01_hlp(x): 
    return 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)

def mortality_acute(v: float, s: float) -> float:
    s_term = (nS_HLP * ONE_PLUS_EPS_HLP * s) / (ONE_PLUS_EPS_HLP - s)
    v_term = (nV_HLP * ONE_PLUS_EPS_HLP * v) / (ONE_PLUS_EPS_HLP - v)
    return d0_HLP + s_term + v_term

def host_fit_acute(v: float, s: float) -> float:
    m = mortality_acute(v, s)
    return s / (s + m)

def path_fit_acute(v: float, s: float) -> float:
    m = mortality_acute(v, s)
    return (v**beta_HLP) / (s + m)

FDH_HLP = 1e-5

def grad_host_hlp(v, s, fH=host_fit_acute):
    v1, v2 = clamp01_hlp(v - FDH_HLP), clamp01_hlp(v + FDH_HLP)
    s1, s2 = clamp01_hlp(s - FDH_HLP), clamp01_hlp(s + FDH_HLP)
    dfdv = (fH(v2, s) - fH(v1, s)) / max(v2 - v1, 1e-12)
    dfds = (fH(v, s2) - fH(v, s1)) / max(s2 - s1, 1e-12)
    return dfdv, dfds

def grad_path_hlp(v, s, fP=path_fit_acute):
    v1, v2 = clamp01_hlp(v - FDH_HLP), clamp01_hlp(v + FDH_HLP)
    s1, s2 = clamp01_hlp(s - FDH_HLP), clamp01_hlp(s + FDH_HLP)
    dfdv = (fP(v2, s) - fP(v1, s)) / max(v2 - v1, 1e-12)
    dfds = (fP(v, s2) - fP(v, s1)) / max(s2 - s1, 1e-12)
    return dfdv, dfds

def read_sim_csv(path: str):
    with open(path, newline="") as f:
        r = csv.DictReader(f);  rows = list(r)
    return rows

# =========================
# Simulation Class
# =========================

class Simulation:
    def __init__(self, evolved_strategy: bool, rng: random.Random):
        self.rng = rng
        self.es = evolved_strategy
        self.zero_rate_streak = 0
        
        self.interior_count_total = 0
        self.boundary_count_total = 0

        self.v = self.rng.random()
        self.s = self.rng.random()

        if self.es:
            # Initialize angles (will be set to 0 if reactivity is fixed)
            if FIX_PATH_REACTIVITY:
                self.v_angle = 0.0
                self.mV = 0.0
            else:
                self.v_angle = wrap_angle(self.rng.uniform(-math.pi/2 + 0.2, math.pi/2 - 0.2))
                self.mV = safe_tan(self.v_angle)
            
            if FIX_HOST_REACTIVITY:
                self.s_angle = 0.0
                self.mS = 0.0
            else:
                self.s_angle = wrap_angle(self.rng.uniform(-math.pi/2 + 0.2, math.pi/2 - 0.2))
                self.mS = safe_tan(self.s_angle)
            
            # Set intercepts based on current (v,s) and slopes
            self.bV = self.v - self.mV * self.s
            self.bS = self.s - self.mS * self.v
        else:
            self.v_angle = 0.0; self.s_angle = 0.0
            self.mV = 0.0; self.mS = 0.0
            self.bV = self.v;  self.bS = self.s

        self.host_fit = host_fitness(self.v, self.s)
        self.path_fit = path_fitness(self.v, self.s)

        self.step_bins = gaussian_step_bins(num_step_bins)

    def propose_host_ES_fixed(self):
        """
        Host ES proposals (original - iterates over both position and angle).
        Use this when host reactivity is NOT fixed.
        """
        proposals: List[Tuple[float, dict]] = []
        cum = 0.0
        neutral_count = 0
        interior_count = 0
        
        # Determine if we should mutate angles
        if FIX_HOST_REACTIVITY:
            # Only mutate position (bS), keep mS=0
            angle_steps = [0.0]  # No angle change
        else:
            # Mutate both position and angle
            angle_steps = self.step_bins
        
        for dz in self.step_bins:
            new_s_guess = clamp01(self.s + std_dev_move * dz)
            for da in angle_steps:
                if FIX_HOST_REACTIVITY:
                    # Keep slope at 0
                    new_angle = 0.0
                    new_mS = 0.0
                else:
                    # Mutate slope
                    new_angle = wrap_angle(self.s_angle + std_dev_angle * da)
                    new_mS = safe_tan(new_angle)
                
                new_bS = new_s_guess - new_mS * self.v
                
                intersections = find_intersections_fast(self.bV, self.mV, new_bS, new_mS)
                
                if intersections is None:
                    continue
                
                neutral_count += 1
                if len(intersections) == 1:
                    interior_count += 1
                
                best_v, best_s = None, None
                best_fit = -1e9
                
                for v2, s2 in intersections:
                    fit = host_fitness(v2, s2)
                    if fit > best_fit:
                        best_fit = fit
                        best_v, best_s = v2, s2
                
                s_coef = selection_coeff(best_fit, self.host_fit)
                w = kimura_intensity(s_coef, Ne_H)
                cum += w
                proposals.append((cum, {
                    "bS": new_bS, "mS": new_mS, "s_angle": new_angle,
                    "v": best_v, "s": best_s,
                    "host_fit": best_fit,
                    "path_fit": path_fitness(best_v, best_s),
                    "sel_coef_mut": s_coef
                }))
        
        self.interior_count_total += interior_count
        self.boundary_count_total += (neutral_count - interior_count)
        return proposals, cum, neutral_count

    def propose_path_ES_fixed(self):
        """
        Pathogen ES proposals (original - iterates over both position and angle).
        Use this when pathogen reactivity is NOT fixed.
        """
        proposals: List[Tuple[float, dict]] = []
        cum = 0.0
        neutral_count = 0
        interior_count = 0
        
        # Determine if we should mutate angles
        if FIX_PATH_REACTIVITY:
            # Only mutate position (bV), keep mV=0
            angle_steps = [0.0]  # No angle change
        else:
            # Mutate both position and angle
            angle_steps = self.step_bins
        
        for dz in self.step_bins:
            new_v_guess = clamp01(self.v + std_dev_move * dz)
            for da in angle_steps:
                if FIX_PATH_REACTIVITY:
                    # Keep slope at 0
                    new_angle = 0.0
                    new_mV = 0.0
                else:
                    # Mutate slope
                    new_angle = wrap_angle(self.v_angle + std_dev_angle * da)
                    new_mV = safe_tan(new_angle)
                
                new_bV = new_v_guess - new_mV * self.s
                
                intersections = find_intersections_fast(new_bV, new_mV, self.bS, self.mS)
                
                if intersections is None:
                    continue
                
                neutral_count += 1
                if len(intersections) == 1:
                    interior_count += 1
                
                best_v, best_s = None, None
                best_fit = -1e9
                
                for v2, s2 in intersections:
                    fit = path_fitness(v2, s2)
                    if fit > best_fit:
                        best_fit = fit
                        best_v, best_s = v2, s2
                
                s_coef = selection_coeff(best_fit, self.path_fit)
                w = kimura_intensity(s_coef, Ne_P)
                cum += w
                proposals.append((cum, {
                    "bV": new_bV, "mV": new_mV, "v_angle": new_angle,
                    "v": best_v, "s": best_s,
                    "path_fit": best_fit,
                    "host_fit": host_fitness(best_v, best_s),
                    "sel_coef_mut": s_coef
                }))
        
        self.interior_count_total += interior_count
        self.boundary_count_total += (neutral_count - interior_count)
        return proposals, cum, neutral_count

    def propose_host_EI(self):
        proposals: List[Tuple[float, dict]] = []
        cum = 0.0
        neutral_count = 0
        
        for dz in self.step_bins:
            v2, s2 = self.v, clamp01(self.s + std_dev_move * dz)
            neutral_count += 1
            new_fit = host_fitness(v2, s2)
            s_coef = selection_coeff(new_fit, self.host_fit)
            w = kimura_intensity(s_coef, Ne_H)
            cum += w
            proposals.append((cum, {
                "v": v2, "s": s2,
                "host_fit": new_fit,
                "path_fit": path_fitness(v2, s2),
                "sel_coef_mut": s_coef
            }))
        
        return proposals, cum, neutral_count

    def propose_path_EI(self):
        proposals: List[Tuple[float, dict]] = []
        cum = 0.0
        neutral_count = 0
        
        for dz in self.step_bins:
            v2, s2 = clamp01(self.v + std_dev_move * dz), self.s
            neutral_count += 1
            new_fit = path_fitness(v2, s2)
            s_coef = selection_coeff(new_fit, self.path_fit)
            w = kimura_intensity(s_coef, Ne_P)
            cum += w
            proposals.append((cum, {
                "v": v2, "s": s2,
                "path_fit": new_fit,
                "host_fit": host_fitness(v2, s2),
                "sel_coef_mut": s_coef
            }))
        
        return proposals, cum, neutral_count

    def step_generation(self):
        if self.es:
            host_props, cum_host, host_neut = self.propose_host_ES_fixed()
            path_props, cum_path, path_neut = self.propose_path_ES_fixed()
        else:
            host_props, cum_host, host_neut = self.propose_host_EI()
            path_props, cum_path, path_neut = self.propose_path_EI()

        omega_host = (cum_host / host_neut) if host_neut > 0 else 0.0
        omega_path = (cum_path / path_neut) if path_neut > 0 else 0.0

        total_rate = host_vs_pathogen_mut_rate * cum_host + (1.0 - host_vs_pathogen_mut_rate) * cum_path
        self.zero_rate_streak = self.zero_rate_streak + 1 if total_rate <= 0.0 else 0
        
        if self.zero_rate_streak > 100:
            print(f"WARNING: Stalled for {self.zero_rate_streak} gens at v={self.v:.3f}, s={self.s:.3f}")

        mutator = "None";  chosen = None
        if total_rate > 0.0:
            r_top = self.rng.random() * total_rate
            if r_top < host_vs_pathogen_mut_rate * cum_host and cum_host > 0.0:
                r = self.rng.random() * cum_host
                idx = bisect_left([x[0] for x in host_props], r);  idx = min(idx, len(host_props)-1)
                state = host_props[idx][1]
                if self.es:
                    self.bS = state["bS"]; self.mS = state["mS"]; self.s_angle = wrap_angle(math.atan(self.mS))
                self.v, self.s = state["v"], state["s"]
                old_host, old_path = self.host_fit, self.path_fit
                self.host_fit, self.path_fit = state["host_fit"], state["path_fit"]
                mutator = "Host"
                s_coef_mut = state["sel_coef_mut"]
                mut_effect = "neutral" if abs(s_coef_mut) < NEUTRAL_THRESH else ("beneficial" if s_coef_mut > 0 else "deleterious")
                chosen = (mutator, s_coef_mut, mut_effect, old_host, old_path)
            else:
                if cum_path > 0.0:
                    r = self.rng.random() * cum_path
                    idx = bisect_left([x[0] for x in path_props], r);  idx = min(idx, len(path_props)-1)
                    state = path_props[idx][1]
                    if self.es:
                        self.bV = state["bV"]; self.mV = state["mV"]; self.v_angle = wrap_angle(math.atan(self.mV))
                    self.v, self.s = state["v"], state["s"]
                    old_host, old_path = self.host_fit, self.path_fit
                    self.path_fit, self.host_fit = state["path_fit"], state["host_fit"]
                    mutator = "Path"
                    s_coef_mut = state["sel_coef_mut"]
                    mut_effect = "neutral" if abs(s_coef_mut) < NEUTRAL_THRESH else ("beneficial" if s_coef_mut > 0 else "deleterious")
                    chosen = (mutator, s_coef_mut, mut_effect, old_host, old_path)

        dwell = self.rng.expovariate(total_rate) if total_rate > 0 else 0.0

        nash_type = "none";  host_max_flag = "";  path_max_flag = "";  stable_fp = ""
        if self.es:
            stable_fp = "stable" if abs(self.mS * self.mV) < 1.0 else "unstable"
            is_weak, hmax, pmax = es_weak_nash(self.v, self.s, self.bV, self.mV, self.bS, self.mS)
            if is_weak: nash_type = "weak"
            host_max_flag = "Hmax" if hmax else "Hnot"
            path_max_flag = "Pmax" if pmax else "Pnot"
        else:
            nash_type = "strong" if ei_strong_nash(self.v, self.s) else "none"

        return {
            "mutator": mutator,
            "dwell": dwell,
            "omega_host": omega_host,
            "omega_path": omega_path,
            "nash": nash_type,
            "host_max_flag": host_max_flag,
            "path_max_flag": path_max_flag,
            "stable_fp": stable_fp,
            "chosen": chosen, 
            "total_rate": total_rate
        }

# =========================
# Runtime Wrapper
# =========================

def run_with_runtime_modes(sim: Simulation, out_csv: str):
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "event","gen","time","mutator",
            "v","bV","vAngle","mV",
            "s","bS","sAngle","mS",
            "pathFit","hostFit",
            "omegaPath","omegaHost",
            "nash","stableFP","hostLineMax","pathLineMax",
            "mutSelCoeff","mutClass","dwell"
        ])

        t = 0.0
        for gen in range(-burn_in_gens, max_gens):
            if gen == 0:  # First generation after burn-in
                t = 0.0   # Reset evolutionary time

            if gen % PROGRESS_EVERY == 0:
                mode = 'ES' if sim.es else 'EI'
                if sim.es:
                    total = sim.interior_count_total + sim.boundary_count_total
                    pct = 100 * sim.interior_count_total / total if total > 0 else 0
                    print(f"[{mode}] gen={gen} t={t:.3e} v={sim.v:.3f} s={sim.s:.3f} interior={pct:.1f}% zr={sim.zero_rate_streak}")
                else:
                    print(f"[{mode}] gen={gen} t={t:.3e} v={sim.v:.3f} s={sim.s:.3f} zr={sim.zero_rate_streak}")

            record_now = _runtime_write_every(gen)

            if record_now:
                w.writerow([
                    "pre", gen, f"{t:.12e}", "NA",
                    f"{sim.v:.6f}", f"{(getattr(sim,'bV',sim.v)):.6f}", f"{sim.v_angle:.6f}", f"{sim.mV:.6f}",
                    f"{sim.s:.6f}", f"{(getattr(sim,'bS',sim.s)):.6f}", f"{sim.s_angle:.6f}", f"{sim.mS:.6f}",
                    f"{sim.path_fit:.6f}", f"{sim.host_fit:.6f}",
                    "", "",
                    "", "", "", "",
                    "", "",
                    ""
                ])

            result = sim.step_generation()
            dwell_out = max(result["dwell"], DWELL_MIN)
            t += dwell_out

            if record_now:
                mutator = result["mutator"]
                mutSelCoeff = "";  mutClass = ""
                if result["chosen"] is not None:
                    _, scoef, cls, _, _ = result["chosen"]
                    mutSelCoeff = f"{scoef:.3e}";  mutClass = cls

                w.writerow([
                    "post", gen, f"{t:.12e}", mutator,
                    f"{sim.v:.6f}", f"{(getattr(sim,'bV',sim.v)):.6f}", f"{sim.v_angle:.6f}", f"{sim.mV:.6f}",
                    f"{sim.s:.6f}", f"{(getattr(sim,'bS',sim.s)):.6f}", f"{sim.s_angle:.6f}", f"{sim.mS:.6f}",
                    f"{sim.path_fit:.6f}", f"{sim.host_fit:.6f}",
                    f"{result['omega_path']:.6f}", f"{result['omega_host']:.6f}",
                    result["nash"], result["stable_fp"], result["host_max_flag"], result["path_max_flag"],
                    mutSelCoeff, mutClass,
                    f"{dwell_out:.6e}"
                ])

# =========================
# Post-processing
# =========================

def write_fig6_timeseries(sim_csv: str, out_csv: str):
    rows = read_sim_csv(sim_csv)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time","v","s","mutator","mutClass","nash"])
        for row in rows:
            if row["event"] != "post":
                continue
            w.writerow([
                f"{float(row['time']):.8f}",
                f"{float(row['v']):.6f}",
                f"{float(row['s']):.6f}",
                row.get("mutator","NA"),
                row.get("mutClass",""),
                row.get("nash","")
            ])

def write_fig5_panel(snapshot: dict, grid_n: int, out_prefix: str,
                     fH=host_fit_acute, fP=path_fit_acute):
    bS = float(snapshot["bS"]); mS = float(snapshot["mS"])
    bV = float(snapshot["bV"]); mV = float(snapshot["mV"])
    v0 = float(snapshot["v"]);  s0 = float(snapshot["s"])

    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    out_csv = f"{out_prefix}.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["v","s","f_H","f_P"])
        for i in range(grid_n):
            v = i/(grid_n-1)
            for j in range(grid_n):
                s = j/(grid_n-1)
                w.writerow([f"{v:.6f}", f"{s:.6f}", f"{fH(v,s):.6f}", f"{fP(v,s):.6f}"])
    meta = {"bS": bS, "mS": mS, "bV": bV, "mV": mV, "intersection": {"v": v0, "s": s0}}
    with open(f"{out_prefix}.meta.json","w") as g:
        json.dump(meta, g, indent=2)

def weak_nash_possible(v: float, s: float,
                       fH=host_fit_acute, fP=path_fit_acute) -> bool:
    ghv, ghs = grad_host_hlp(v,s,fH);  gpv, gps = grad_path_hlp(v,s,fP)
    if abs(ghv)+abs(ghs) < 1e-12 or abs(gpv)+abs(gps) < 1e-12:
        return False
    if abs(gps) < 1e-12 or abs(ghv) < 1e-12:
        return False
    mS_star = -gpv / gps
    mV_star = -ghs / ghv
    if abs(mS_star) > 1e3 or abs(mV_star) > 1e3:
        return False
    bS_star = s - mS_star * v
    bV_star = v - mV_star * s
    def H(vx): return clamp01(bS_star + mS_star * vx)
    def P(sx): return clamp01(bV_star + mV_star * sx)
    s1 = H(v); v1 = P(s1)
    if abs(v1 - v) + abs(s1 - s) > 1e-3:
        return False
    dv = 1e-4; ds = 1e-4
    vL, vR = clamp01(v - dv), clamp01(v + dv)
    sL, sR = H(vL), H(vR)
    P0, PL, PR = fP(v,s), fP(vL,sL), fP(vR,sR)
    if not (P0 >= PL - 1e-10 and P0 >= PR - 1e-10):
        return False
    sL2, sR2 = clamp01(s - ds), clamp01(s + ds)
    vL2, vR2 = P(sL2), P(sR2)
    H0, HL, HR = fH(v,s), fH(vL2,sL2), fH(vR2,sR2)
    return (H0 >= HL - 1e-10 and H0 >= HR - 1e-10)

def write_fig7_region(grid_n: int, out_csv: str,
                      fH=host_fit_acute, fP=path_fit_acute):
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f);  w.writerow(["v","s","weak_nash_possible"])
        for i in range(grid_n):
            v = i/(grid_n-1)
            for j in range(grid_n):
                s = j/(grid_n-1)
                w.writerow([f"{v:.6f}", f"{s:.6f}", int(weak_nash_possible(v,s,fH,fP))])

def write_fig7_occupancy(sim_csv: str, out_csv: str,
                         dwell_threshold: float = 0.01,
                         bins: int = 50):
    rows = read_sim_csv(sim_csv)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    H = [[0.0 for _ in range(bins)] for _ in range(bins)]
    last_t = None
    for row in rows:
        if row["event"] != "post":
            continue
        t = float(row["time"])
        v = float(row["v"]); s = float(row["s"])
        dwell = row.get("dwell","")
        if dwell == "" or dwell.lower()=="nan":
            dwell_val = t - (last_t if last_t is not None else t)
        else:
            dwell_val = max(float(dwell), DWELL_MIN)
        last_t = t
        if dwell_val < dwell_threshold:
            continue
        i = min(bins-1, max(0, int(v * bins)))
        j = min(bins-1, max(0, int(s * bins)))
        H[i][j] += dwell_val
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f);  w.writerow(["v_bin_mid","s_bin_mid","occupancy"])
        for i in range(bins):
            for j in range(bins):
                v_mid = (i + 0.5)/bins
                s_mid = (j + 0.5)/bins
                w.writerow([f"{v_mid:.6f}", f"{s_mid:.6f}", f"{H[i][j]:.6e}"])

def pick_snapshots_for_fig5(sim_csv: str) -> List[dict]:
    rows = read_sim_csv(sim_csv)
    weak = [r for r in rows if r["event"]=="post" and r.get("nash","")=="weak" and r.get("stableFP","")=="stable"]
    host_neu = [r for r in rows if r["event"]=="post" and r.get("mutator","")=="Host" and r.get("mutClass","")=="neutral"]
    path_neu = [r for r in rows if r["event"]=="post" and r.get("mutator","")=="Path" and r.get("mutClass","")=="neutral"]
    def to_snap(r):
        return {"bS": r.get("bS","0"), "mS": r.get("mS","0"),
                "bV": r.get("bV","0"), "mV": r.get("mV","0"),
                "v": r.get("v","0"), "s": r.get("s","0")}
    picks = []
    if weak: picks.append(to_snap(weak[0]))
    if host_neu: picks.append(to_snap(host_neu[0]))
    if path_neu: picks.append(to_snap(path_neu[0]))
    return picks

def validate_simulation_results(csv_path: str, mode: str):
    rows = read_sim_csv(csv_path)
    post_rows = [r for r in rows if r["event"] == "post"]
    
    if len(post_rows) == 0:
        return {"error": "No data found"}
    
    v_vals = [float(r["v"]) for r in post_rows]
    s_vals = [float(r["s"]) for r in post_rows]
    times = [float(r["time"]) for r in post_rows]
    
    v_min, v_max = min(v_vals), max(v_vals)
    s_min, s_max = min(s_vals), max(s_vals)
    
    v_unique = len(set([round(v, 3) for v in v_vals]))
    s_unique = len(set([round(s, 3) for s in s_vals]))
    
    mutations = [r.get("mutator", "None") for r in post_rows]
    n_host = mutations.count("Host")
    n_path = mutations.count("Path")
    
    report = {
        "mode": mode,
        "n_generations": len(post_rows),
        "final_time": times[-1] if times else 0,
        "v_range": (v_min, v_max),
        "s_range": (s_min, s_max),
        "v_diversity": v_unique,
        "s_diversity": s_unique,
        "host_mutations": n_host,
        "path_mutations": n_path,
        "mutation_ratio": n_path / n_host if n_host > 0 else float('inf')
    }
    
    warnings = []
    if v_unique < 10 or s_unique < 10:
        warnings.append("LOW DIVERSITY: Traits may be stuck at boundaries")
    if report["mutation_ratio"] < 50 or report["mutation_ratio"] > 150:
        warnings.append(f"UNEXPECTED RATIO: Path/Host mutations = {report['mutation_ratio']:.1f} (expect ~99)")
    if v_min < 0.01 or v_max > 0.99 or s_min < 0.01 or s_max > 0.99:
        warnings.append("BOUNDARY HUGGING: Traits spending time at edges")
    
    report["warnings"] = warnings
    return report

# =========================
# Main
# =========================

def main(run_postproc=True, diag=False):
    print("=" * 70)
    if USE_GAUSSIAN:
        print("GAUSSIAN FITNESS MODE - Reproducing Paper Results")
        print("  Gentle bowl-shaped landscape, weak selection")
        print("  Expected: Clear stasis + bursts, omega ~1-10")
    else:
        print("BIOLOGICAL FITNESS MODE - Acute Infection Model")
        print(f"  Steep costs: nS={nS}, nV={nV}, d0={d0}, β={beta}")
        print("  Expected: Strong selection, omega >1000, no clear stasis")
    print("=" * 70)
    print(f"Populations: Ne_H={Ne_H:.0e}, Ne_P={Ne_P:.0e}")
    print(f"Mutations: γ={host_vs_pathogen_mut_rate} (host:pathogen = 1:99)")
    print(f"Generations: {max_gens:,} (burn-in: {burn_in_gens:,})")
    print(f"Runtime mode: {RUNTIME_MODE}")
    print(f"Random seed: {seed}")
    print("=" * 70)
    
    if diag:
        print("\n[Diagnostic Mode]")
        print(f"Host fitness at (0.5, 0.5): {host_fitness(0.5, 0.5):.6e}")
        print(f"Path fitness at (0.5, 0.5): {path_fitness(0.5, 0.5):.6e}")
        bv, mv, bs, ms = 0.2, 0.5, 0.3, -0.4
        sol = find_intersections_fast(bv, mv, bs, ms)
        print(f"Equilibrium test: {sol}")
        return

    os.makedirs("results", exist_ok=True)
    rng = random.Random(seed)
    
    import datetime
    config = {
        "version": __version__,
        "model": "gaussian" if USE_GAUSSIAN else "acute_infection",
        "fitness_mode": "Gaussian bowl" if USE_GAUSSIAN else "Biological acute",
        "USE_GAUSSIAN": USE_GAUSSIAN,
        "timestamp": datetime.datetime.now().isoformat(),
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "parameters": {
            "Ne_H": int(Ne_H), "Ne_P": int(Ne_P),
            "gamma": host_vs_pathogen_mut_rate,
            "std_dev_move": std_dev_move,
            "std_dev_angle": std_dev_angle / math.pi,
        },
        "simulation": {
            "burn_in_gens": burn_in_gens,
            "max_gens": max_gens,
            "seed": seed,
            "runtime_mode": RUNTIME_MODE,
            "num_step_bins": num_step_bins
        }
    }
    
    if not USE_GAUSSIAN:
        config["parameters"].update({
            "nS": nS, "nV": nV, "d0": d0, "beta": beta
        })
    
    with open("results/run_config.json", "w") as f:
        json.dump(config, f, indent=2)
    print(f"Configuration saved to results/run_config.json\n")

    print("Running ES (Evolved Strategy) simulation...")
    es = Simulation(evolved_strategy=True, rng=rng)
    run_with_runtime_modes(es, "results/es_simulation.csv")
    if hasattr(es, 'interior_count_total'):
        print(f"ES Stats: Interior={es.interior_count_total}, Boundary={es.boundary_count_total}")

    print("Running EI (Evolved Interaction) simulation...")
    ei = Simulation(evolved_strategy=False, rng=random.Random(seed+1))
    run_with_runtime_modes(ei, "results/ei_simulation.csv")

    print("\nSimulations complete!")
    
    print("\n" + "="*70)
    print("VALIDATION")
    print("="*70)
    
    for csv_file, mode in [("results/es_simulation.csv", "ES"), 
                           ("results/ei_simulation.csv", "EI")]:
        if os.path.exists(csv_file):
            report = validate_simulation_results(csv_file, mode)
            print(f"\n{mode} Mode Summary:")
            print(f"  Generations: {report['n_generations']:,}")
            print(f"  Final time: {report['final_time']:.2e}")
            print(f"  v range: [{report['v_range'][0]:.3f}, {report['v_range'][1]:.3f}]")
            print(f"  s range: [{report['s_range'][0]:.3f}, {report['s_range'][1]:.3f}]")
            print(f"  Unique v values: {report['v_diversity']}")
            print(f"  Unique s values: {report['s_diversity']}")
            print(f"  Mutations - Host: {report['host_mutations']}, Path: {report['path_mutations']}")
            print(f"  Path/Host ratio: {report['mutation_ratio']:.1f}× (expect ~99×)")
            
            if report['warnings']:
                print(f"WARNINGS:")
                for w in report['warnings']:
                    print(f"     - {w}")
            else:
                print(f"All checks passed")
    
    print("="*70 + "\n")

    if run_postproc:
        print("Running post-processing...")
        os.makedirs("postproc", exist_ok=True)
        write_fig6_timeseries("results/es_simulation.csv", "postproc/fig6_timeseries.csv")
        grid_n, bins = _postproc_sizes()
        snaps = pick_snapshots_for_fig5("results/es_simulation.csv")
        for k, snap in enumerate(snaps, start=1):
            write_fig5_panel(snap, grid_n=grid_n, out_prefix=f"postproc/fig5_panel_{k}")
        write_fig7_region(grid_n=grid_n, out_csv="postproc/fig7_region_acute.csv")
        write_fig7_occupancy("results/es_simulation.csv", "postproc/fig7_occupancy.csv",
                             dwell_threshold=0.01, bins=bins)
        print("Post-processing done.")

if __name__ == "__main__":
    main(run_postproc=True, diag=False)