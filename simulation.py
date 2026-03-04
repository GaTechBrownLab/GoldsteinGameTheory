# =========================
# simulation.py
# Unified host-pathogen coevolution simulation
# Supports fitness models: acute, minimal, taylor, chronic
# Compatible with run_experiments.py
# Author: Canan Karakoc, Karan Gosrani
# Origin: Jan Goldstein's Java code 2020
# =========================
#
# Notation mapping (paper -> code):
#   ET (Evolved Trait)       -> EI / evolved_strategy=False
#   ER (Evolved Response)    -> ES / evolved_strategy=True
#   c  (clearance)           -> s
#   v  (virulence)           -> v
#   m_c (host slope)         -> mS
#   m_v (pathogen slope)     -> mV
#   c_0 (host intercept)     -> bS
#   v_0 (pathogen intercept) -> bV
#   W_H                      -> hostFit
#   W_P                      -> pathFit

from __future__ import annotations

import csv
import math
import os
import random
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

try:
    from scipy.special import erfinv as _erfinv
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ============================================================
# Runner-overridden globals (run_experiments.py sets these)
# ============================================================

burn_in_gens = 10_000
max_gens = 1_000_000
seed = 3248232

RUNTIME_MODE = "fast"
WRITE_EVERY_FAST = 100
write_every = 1
PROGRESS_EVERY = 50_000
DWELL_MIN = 1e-12

FIX_HOST_REACTIVITY = False
FIX_PATH_REACTIVITY = False

# Fix a player's trait at a constant value (None = evolves normally).
# When set, that player never mutates; their trait stays at this value.
FIX_HOST_TRAIT: Optional[float] = None   # e.g. 0.5 to pin s=0.5
FIX_PATH_TRAIT: Optional[float] = None   # e.g. 0.5 to pin v=0.5

USE_BOUNDED_TRAITS = True
USE_GAUSSIAN = False
DIPLOID_KIMURA = False  # True: 4Ns denominator (diploid semi-dominant); False: 2Ns (haploid)

# Model-specific trait domain.  set_fitness_model() adjusts these.
TRAIT_MIN = 0.0
TRAIT_MAX = 1.0

# ============================================================
# Evolution / mutation controls
# ============================================================

num_step_bins = 51
std_dev_move = 0.1
std_dev_angle = 0.1 * math.pi

# Mutation asymmetry: probability of host mutation per event.
# Default 0.01 reflects ~100:1 pathogen-to-host mutation ratio
# (Goldstein 2020). gamma=1.0 gives equal rates.
prob_host_mutate = 0.01

NEUTRAL_THRESH = 1e-12

# Separate effective population sizes (haploid)
HOST_POP_N = 10_000
PATH_POP_N = 1_000_000

ANGLE_EPS = 1e-4
TINY = 1e-12

# ============================================================
# Model-specific parameters
# ============================================================

# --- Acute model ---
# mortality_acute(v,s) = d0 + nS*(1+eps)*s/(1+eps-s) + nV*(1+eps)*v/(1+eps-v)
d0_HLP = 0.1
nS_HLP = 0.1
nV_HLP = 1
eps_HLP = 1e-3
ONE_PLUS_EPS_HLP = 1.0 + eps_HLP
beta_HLP = 1.0  # exponent on v in pathogen transmission-like term

# --- Chronic model (Goldstein 2020, §Methods) ---
# Mortality includes immunity-modulated virulence: (1-s) dampens v damage
# mortality_chronic = d0 + nS*(1+eps)*s/(1+eps-s) + (1-s)*nV*(1+eps)*v/(1+eps-v)
# W_H = 1 / mortality              (expected lifetime)
# W_P = (1-s)*v^beta / mortality    (transmission × lifetime)
# "Unless stated otherwise, the parameters are the same as used for acute infections."


# --- Taylor model ---
# Host: H(v,c) = [c/(v+c)] * [b/(m0+c)]
# Path: P(v,c) = v^n / (v+c)
# Taylor et al. 2006, Eq. 2.2: n=3/4, m0=b=1
b_paper = 1.0
m0_paper = 1.0     #  Taylor paper uses m0=b=1
n_paper = 0.75     #  paper: "n between 0 and 1", uses 3/4

# ============================================================
# Helpers
# ============================================================

def clamp_trait(x: float) -> float:
    """Clamp x to [TRAIT_MIN, TRAIT_MAX]."""
    if x <= TRAIT_MIN:
        return TRAIT_MIN
    if x >= TRAIT_MAX:
        return TRAIT_MAX
    return x

# Backward-compatible alias (used internally everywhere)
def clamp01(x: float) -> float:
    return clamp_trait(x)

def clamp01_or_unbounded(x: float) -> float:
    return clamp_trait(x) if USE_BOUNDED_TRAITS else x

def angle_to_slope(theta: float) -> float:
    """Convert angle to slope. Angles live on [0, pi) so slopes span all reals."""
    th = max(min(theta, math.pi - ANGLE_EPS), ANGLE_EPS)
    return math.tan(th)

def wrap_angle(theta: float) -> float:
    """Wrap angle to [0, pi) — Java-style wrapping where slopes cycle through +/-inf."""
    return (theta + 2.0 * math.pi) % math.pi


def make_equal_prob_steps(n_bins: int) -> List[float]:
    """
    Reproduce Java makeSteps(): compute expected z-scores within each of
    n_bins equal-probability Gaussian quantile bins.

    Each step represents the expected value of Z within an equal-probability
    slice of the standard normal, so all steps are equally likely.
    Falls back to linearly-spaced z-scores if scipy is unavailable.
    """
    if not HAS_SCIPY or n_bins < 3:
        # Fallback: linearly spaced z-scores
        zmax = 3.0
        return [(-zmax + 2 * zmax * i / (n_bins - 1)) for i in range(n_bins)]

    sqrt2pi = math.sqrt(2.0 * math.pi)
    sqrt2 = math.sqrt(2.0)

    # Compute quantile boundaries (divides[0..n_bins])
    divides = [0.0] * (n_bins + 1)
    for i in range(1, n_bins):
        divides[i] = sqrt2 * _erfinv((2.0 * i) / n_bins - 1.0)
    # Extend tails
    divides[0] = divides[1] * 10.0
    divides[n_bins] = divides[n_bins - 1] * 10.0

    # Expected z within each bin: (phi(lo) - phi(hi)) / (1/n_bins * sqrt(2pi))
    steps = []
    for i in range(n_bins):
        lo, hi = divides[i], divides[i + 1]
        val = (math.exp(-lo * lo / 2.0) - math.exp(-hi * hi / 2.0)) / (sqrt2pi / n_bins)
        steps.append(val)

    return steps


def kimura_fixation_prob(scoef: float, N: int) -> float:
    """
    Kimura fixation probability.
    Haploid:           (1 - exp(-2s)) / (1 - exp(-2Ns))
    Diploid semi-dom:  (1 - exp(-2s)) / (1 - exp(-4Ns))
    Near-neutral (|s| < threshold) returns 1/N (neutral drift).
    """
    if abs(scoef) <= NEUTRAL_THRESH:
        return 1.0 / N   # Neutral drift: fixation prob = 1/N
    factor = 4.0 if DIPLOID_KIMURA else 2.0
    x = -factor * N * scoef
    if x > 700:
        return 0.0        # Strongly deleterious
    if x < -700:
        return 2.0 * scoef  # Strongly beneficial: approx 2s
    denom = 1.0 - math.exp(x)
    if abs(denom) < TINY:
        return 1.0 / N
    return (1.0 - math.exp(-2.0 * scoef)) / denom


def kimura_rate(scoef: float, N: int) -> float:
    """
    Substitution rate = N * P_fix.
    For equal-probability step bins, each step has equal mutation probability
    (1/n_bins), so the rate weight is N * P_fix.
    Neutral mutations contribute rate = 1 (N * 1/N = 1).
    """
    return N * kimura_fixation_prob(scoef, N)


# ============================================================
# Fitness functions (keyed by model name)
# ============================================================
# All fitness functions take (v, s) where s = clearance
# (called 'c' in the Taylor paper notation).

def _mortality_acute(v: float, s: float) -> float:
    dv = max(ONE_PLUS_EPS_HLP - v, 1e-12)
    ds = max(ONE_PLUS_EPS_HLP - s, 1e-12)
    s_term = (nS_HLP * ONE_PLUS_EPS_HLP * s) / ds
    v_term = (nV_HLP * ONE_PLUS_EPS_HLP * v) / dv
    return d0_HLP + s_term + v_term

def _host_acute(v: float, s: float) -> float:
    m = _mortality_acute(v, s)
    return s / (s + m) if (s + m) > 1e-12 else 0.0

def _path_acute(v: float, s: float) -> float:
    m = _mortality_acute(v, s)
    denom = s + m
    if denom <= 1e-12:
        return 0.0
    return (v ** beta_HLP) / denom

def _host_minimal(v: float, s: float) -> float:
    # wh = c(1-c)(1-v)
    return s * (1.0 - s) * (1.0 - v)

def _path_minimal(v: float, s: float) -> float:
    # wp = v(1-v)(1-c)
    return v * (1.0 - v) * (1.0 - s)

def _mortality_chronic(v: float, s: float) -> float:
    """Chronic model: immunity modulates virulence damage via (1-s) factor.
    Uses same parameters as acute (d0, nS, nV, eps from _HLP)."""
    dv = max(ONE_PLUS_EPS_HLP - v, 1e-12)
    ds = max(ONE_PLUS_EPS_HLP - s, 1e-12)
    s_term = (nS_HLP * ONE_PLUS_EPS_HLP * s) / ds
    v_term = (1.0 - s) * (nV_HLP * ONE_PLUS_EPS_HLP * v) / dv
    return d0_HLP + s_term + v_term

def _host_chronic(v: float, s: float) -> float:
    """W_H = 1/m (expected lifetime). Chronic: no clearance, just survival."""
    m = _mortality_chronic(v, s)
    return 1.0 / m if m > 1e-12 else 1e12

def _path_chronic(v: float, s: float) -> float:
    """W_P = (1-s)*v^beta / m (transmission × expected lifetime)."""
    m = _mortality_chronic(v, s)
    if m <= 1e-12:
        return 0.0
    return (1.0 - s) * (v ** beta_HLP) / m

def _host_taylor(v: float, s: float) -> float:
    """s = clearance (called 'c' in Taylor et al.)"""
    denom1 = v + s
    denom2 = m0_paper + s
    if denom1 <= TINY or denom2 <= TINY:
        return 0.0
    return (s / denom1) * (b_paper / denom2)

def _path_taylor(v: float, s: float) -> float:
    """s = clearance (called 'c' in Taylor et al.)"""
    denom = v + s
    if denom <= TINY:
        return 0.0
    return (v ** n_paper) / denom


# Registry: model name -> (host_fitness, path_fitness)
FITNESS_FUNCS: Dict[str, Tuple[Callable, Callable]] = {
    "acute":   (_host_acute,   _path_acute),
    "chronic": (_host_chronic, _path_chronic),
    "minimal": (_host_minimal, _path_minimal),
    "taylor":  (_host_taylor,  _path_taylor),
}

# Module-level active fitness functions.
host_fitness: Callable[[float, float], float] = _host_acute
path_fitness: Callable[[float, float], float] = _path_acute
FITNESS_MODEL = "acute"

def set_fitness_model(model: str) -> None:
    """Switch the active fitness functions and trait domain."""
    global host_fitness, path_fitness, FITNESS_MODEL, TRAIT_MIN, TRAIT_MAX
    if model not in FITNESS_FUNCS:
        raise ValueError(f"Unknown fitness model '{model}'. Choose from: {list(FITNESS_FUNCS)}")
    host_fitness, path_fitness = FITNESS_FUNCS[model]
    FITNESS_MODEL = model
    # Taylor traits are rates (unbounded above 1); all others are proportions in [0,1]
    if model == "taylor":
        TRAIT_MIN = 0.001   # avoid division by zero at v=0,s=0
        TRAIT_MAX = 20.0    # Nash ≈ (v*=9, c*=3) for n=3/4, m0=1
    else:
        TRAIT_MIN = 0.0
        TRAIT_MAX = 1.0

# ============================================================
# ER equilibrium solving with clamping
# ============================================================

@dataclass(frozen=True)
class Equilibrium:
    v: float
    s: float
    interior: bool
    stable: bool
    host_max: bool
    path_max: bool
    nash: bool

def _best_equilibrium_for_player(eqs: List[Equilibrium], player: str) -> Optional[Equilibrium]:
    if not eqs:
        return None
    if player == "host":
        key = lambda e: host_fitness(e.v, e.s)
    else:
        key = lambda e: path_fitness(e.v, e.s)
    return max(eqs, key=key)

def _solve_interior(bS: float, mS: float, bV: float, mV: float) -> Optional[Tuple[float, float]]:
    denom = (1.0 - mS*mV)
    if abs(denom) < 1e-10:
        return None
    s = (bS + mS*bV) / denom
    v = bV + mV*s
    return (v, s)

def _clamped_line_value(b: float, m: float, x: float) -> float:
    return clamp01(b + m*x)

def find_all_equilibria(bS: float, mS: float, bV: float, mV: float) -> List[Equilibrium]:
    eqs: List[Equilibrium] = []

    interior = _solve_interior(bS, mS, bV, mV)

    if not USE_BOUNDED_TRAITS:
        # Unbounded: always use interior intersection (it's the unique crossing)
        if interior is not None:
            v0, s0 = interior
            stable = (abs(mS * mV) < 1.0)
            eqs.append(Equilibrium(v=v0, s=s0, interior=True, stable=stable,
                                   host_max=False, path_max=False, nash=True))
        return eqs

    # Bounded mode: always use interior if in bounds.
    # The simultaneous solution s=clamp(bS+mS*v), v=clamp(bV+mV*s) is exact
    # when (v0,s0) is in [0,1] — clamping has no effect on interior points.
    # Stability label tracks sensitivity but doesn't affect validity.
    if interior is not None:
        v0, s0 = interior
        if TRAIT_MIN < v0 < TRAIT_MAX and TRAIT_MIN < s0 < TRAIT_MAX:
            stable = (abs(mS*mV) < 1.0)
            eqs.append(Equilibrium(v=v0, s=s0, interior=True, stable=stable,
                                   host_max=False, path_max=False, nash=True))

    candidates: List[Tuple[float, float]] = []
    lo, hi = TRAIT_MIN, TRAIT_MAX
    candidates.extend([(lo,lo),(lo,hi),(hi,lo),(hi,hi)])

    for vfix in [lo, hi]:
        s = _clamped_line_value(bS, mS, vfix)
        v = _clamped_line_value(bV, mV, s)
        candidates.append((v, s))

    for sfix in [lo, hi]:
        v = _clamped_line_value(bV, mV, sfix)
        s = _clamped_line_value(bS, mS, v)
        candidates.append((v, s))

    uniq: List[Tuple[float, float]] = []
    for v, s in candidates:
        if all(abs(v-v2) > 1e-6 or abs(s-s2) > 1e-6 for v2, s2 in uniq):
            uniq.append((v, s))

    for v, s in uniq:
        s_check = _clamped_line_value(bS, mS, v)
        v_check = _clamped_line_value(bV, mV, s)
        if abs(s-s_check) < 1e-6 and abs(v-v_check) < 1e-6:
            host_max = (s <= TRAIT_MIN + 1e-9 or s >= TRAIT_MAX - 1e-9)
            path_max = (v <= TRAIT_MIN + 1e-9 or v >= TRAIT_MAX - 1e-9)
            eqs.append(Equilibrium(v=v, s=s, interior=False, stable=True,
                                   host_max=host_max, path_max=path_max, nash=True))

    # Fallback: when no equilibria found (unstable interior + no boundary
    # candidates pass), iterate clamped maps to find 2-cycle states.
    # These are the biologically relevant outcomes — boundary-clamped states
    # where the system actually rests (stasis).
    if len(eqs) == 0:
        eqs.extend(_solve_clamped_cycle(bS, mS, bV, mV))

    return eqs


def _solve_clamped_cycle(bS: float, mS: float, bV: float, mV: float,
                          max_iter: int = 200, tol: float = 1e-9) -> List[Equilibrium]:
    """Iterate clamped response maps to find fixed point or 2-cycle states.

    Tries multiple starting points (Brouwer guarantees a fixed point exists).
    When |mS*mV| > 1, the clamped map may produce a 2-cycle instead of
    converging.  Returns BOTH cycle states (labeled nash=False, since they
    only satisfy one player's best-response) so the selecting player can
    choose the one maximizing their fitness.
    """
    import warnings
    lo, hi = TRAIT_MIN, TRAIT_MAX
    starts = [0.5 * (lo + hi), lo + 0.01, hi - 0.01, 0.25, 0.75]

    for v0 in starts:
        v = v0
        for _ in range(max_iter):
            s_new = clamp01(bS + mS * v)
            v_new = clamp01(bV + mV * s_new)
            if abs(v_new - v) < tol:
                # Converged to a true fixed point
                s_final = clamp01(bS + mS * v_new)
                host_max = (s_final <= TRAIT_MIN + 1e-9 or s_final >= TRAIT_MAX - 1e-9)
                path_max = (v_new <= TRAIT_MIN + 1e-9 or v_new >= TRAIT_MAX - 1e-9)
                return [Equilibrium(v=v_new, s=s_final, interior=False,
                                    stable=True, host_max=host_max,
                                    path_max=path_max, nash=True)]
            v = v_new

    # No starting point converged -> extract 2-cycle states.
    # These are NOT Nash equilibria (each only satisfies one player's
    # best-response), so nash=False.
    warnings.warn(f"Clamped iteration 2-cycle fallback: bS={bS:.4f}, mS={mS:.4f}, "
                  f"bV={bV:.4f}, mV={mV:.4f}", stacklevel=3)
    v_a = v  # last iterate from final starting point
    s_a = clamp01(bS + mS * v_a)
    v_b = clamp01(bV + mV * s_a)
    s_b = clamp01(bS + mS * v_b)

    results = []
    for (vc, sc) in [(v_a, s_a), (v_b, s_b)]:
        host_max = (sc <= TRAIT_MIN + 1e-9 or sc >= TRAIT_MAX - 1e-9)
        path_max = (vc <= TRAIT_MIN + 1e-9 or vc >= TRAIT_MAX - 1e-9)
        results.append(Equilibrium(v=vc, s=sc, interior=False, stable=False,
                                   host_max=host_max, path_max=path_max, nash=False))

    # Deduplicate if both states are the same (degenerate 2-cycle = fixed point)
    if abs(results[0].v - results[1].v) < 1e-6 and abs(results[0].s - results[1].s) < 1e-6:
        results[0] = Equilibrium(v=results[0].v, s=results[0].s, interior=False,
                                 stable=True, host_max=results[0].host_max,
                                 path_max=results[0].path_max, nash=True)
        return [results[0]]

    return results

# ============================================================
# Simulation (Gillespie architecture)
# ============================================================

class Simulation:
    """
    EI (es=False): traits evolve directly (v, s)       -- ET in paper
    ES (es=True):  rules evolve (bS,mS,bV,mV),         -- ER in paper
                   realized phenotype is equilibrium (v,s)

    Gillespie dynamics: each generation evaluates ALL mutations for BOTH
    players, computes cumulative substitution rates, chooses who mutates
    proportional to total rate, and draws exponential dwell time.
    """

    def __init__(self, evolved_strategy: bool, rng: Optional[random.Random] = None,
                 model: Optional[str] = None):
        self.es = evolved_strategy
        self.rng = rng or random.Random(seed)

        if model is not None:
            set_fitness_model(model)

        # Pre-compute step bins (equal-probability Gaussian quantiles)
        self._trait_steps = make_equal_prob_steps(num_step_bins)
        self._angle_steps = make_equal_prob_steps(num_step_bins)

        # Initialize traits randomly within [TRAIT_MIN, TRAIT_MAX]
        if FIX_PATH_TRAIT is not None:
            self.v = FIX_PATH_TRAIT
        else:
            self.v = TRAIT_MIN + self.rng.random() * (TRAIT_MAX - TRAIT_MIN)

        if FIX_HOST_TRAIT is not None:
            self.s = FIX_HOST_TRAIT
        else:
            self.s = TRAIT_MIN + self.rng.random() * (TRAIT_MAX - TRAIT_MIN)

        # Initialize ER parameters
        self.s_angle = 0.0
        self.v_angle = 0.0
        if self.es:
            if not FIX_HOST_REACTIVITY:
                self.s_angle = self.rng.random() * math.pi
            if not FIX_PATH_REACTIVITY:
                self.v_angle = self.rng.random() * math.pi
        self.mS = angle_to_slope(self.s_angle)
        self.mV = angle_to_slope(self.v_angle)
        # Intercepts derived so response lines pass through current (v, s)
        self.bS = self.s - self.mS * self.v
        self.bV = self.v - self.mV * self.s

        self.path_fit = path_fitness(self.v, self.s)
        self.host_fit = host_fitness(self.v, self.s)
        self.zero_rate_streak = 0

        self.interior_count_total = 0
        self.boundary_count_total = 0

        if self.es:
            self._refresh_equilibrium(selector="host", track=True)

    def _refresh_equilibrium(self, selector: str, track: bool = False):
        eqs = find_all_equilibria(self.bS, self.mS, self.bV, self.mV)
        best = _best_equilibrium_for_player(eqs, selector)
        if best is None:
            # Should be very rare now with 2-cycle fallback
            import warnings
            warnings.warn(f"No equilibrium found: bS={self.bS:.4f}, mS={self.mS:.4f}, "
                          f"bV={self.bV:.4f}, mV={self.mV:.4f}. Using intercepts.",
                          stacklevel=2)
            self.v, self.s = clamp01(self.bV), clamp01(self.bS)
            interior = False
            stable = False
            host_max = False
            path_max = False
            nash = False
        else:
            self.v, self.s = best.v, best.s
            interior = best.interior
            stable = best.stable
            host_max = best.host_max
            path_max = best.path_max
            nash = best.nash

        self.path_fit = path_fitness(self.v, self.s)
        self.host_fit = host_fitness(self.v, self.s)

        if track:
            if interior:
                self.interior_count_total += 1
            else:
                self.boundary_count_total += 1

        return {"nash": str(nash), "stable_fp": str(stable),
                "host_max_flag": str(host_max), "path_max_flag": str(path_max)}

    # ----------------------------------------------------------
    # ET mutation proposals
    # ----------------------------------------------------------

    def _propose_ET_mutants_host(self) -> List[Tuple[float, float]]:
        """Host ET: mutate s by step * std_dev_move, keep v fixed."""
        muts = []
        for step in self._trait_steps:
            new_s = clamp01_or_unbounded(self.s + std_dev_move * step)
            muts.append((self.v, new_s))
        return muts

    def _propose_ET_mutants_path(self) -> List[Tuple[float, float]]:
        """Pathogen ET: mutate v by step * std_dev_move, keep s fixed."""
        muts = []
        for step in self._trait_steps:
            new_v = clamp01_or_unbounded(self.v + std_dev_move * step)
            muts.append((new_v, self.s))
        return muts

    # ----------------------------------------------------------
    # ER mutation proposals (Java-style: mutate trait, derive intercept)
    # ----------------------------------------------------------

    def _propose_ER_mutants_host(self) -> List[Tuple]:
        """
        Host ER mutations following Java/paper approach:
          - new_s = s + std_dev_move * trait_step    (mutate trait value)
          - new_angle = wrap(s_angle + std_dev_angle * angle_step)
          - new_mS = tan(new_angle)
          - new_bS = new_s - new_mS * v   (intercept derived: line passes through (v, new_s))
        Returns list of (new_bS, new_s_angle, new_mS).
        If FIX_HOST_REACTIVITY, only trait mutates (angle fixed).
        """
        muts = []
        for t_step in self._trait_steps:
            new_s = self.s + std_dev_move * t_step
            if FIX_HOST_REACTIVITY:
                new_angle = self.s_angle
                new_mS = self.mS
                new_bS = new_s - new_mS * self.v
                muts.append((new_bS, new_angle, new_mS))
            else:
                for a_step in self._angle_steps:
                    new_angle = wrap_angle(self.s_angle + std_dev_angle * a_step)
                    new_mS = angle_to_slope(new_angle)
                    new_bS = new_s - new_mS * self.v
                    muts.append((new_bS, new_angle, new_mS))
        return muts

    def _propose_ER_mutants_path(self) -> List[Tuple]:
        """
        Pathogen ER mutations (symmetric to host):
          - new_v = v + std_dev_move * trait_step
          - new_angle = wrap(v_angle + std_dev_angle * angle_step)
          - new_mV = tan(new_angle)
          - new_bV = new_v - new_mV * s   (line passes through (s, new_v))
        """
        muts = []
        for t_step in self._trait_steps:
            new_v = self.v + std_dev_move * t_step
            if FIX_PATH_REACTIVITY:
                new_angle = self.v_angle
                new_mV = self.mV
                new_bV = new_v - new_mV * self.s
                muts.append((new_bV, new_angle, new_mV))
            else:
                for a_step in self._angle_steps:
                    new_angle = wrap_angle(self.v_angle + std_dev_angle * a_step)
                    new_mV = angle_to_slope(new_angle)
                    new_bV = new_v - new_mV * self.s
                    muts.append((new_bV, new_angle, new_mV))
        return muts

    # ----------------------------------------------------------
    # Evaluate all mutations for one player, return candidates + cumulative rate
    # ----------------------------------------------------------

    def _evaluate_host_mutations(self) -> Tuple[List[Tuple], float]:
        """Returns (candidate_list, cumulative_rate) for host."""
        candidates = []
        cum_rate = 0.0
        current_fit = self.host_fit

        if not self.es:
            for (v2, s2) in self._propose_ET_mutants_host():
                f2 = host_fitness(v2, s2)
                scoef = (f2 - current_fit) / max(current_fit, TINY) if current_fit > TINY else 0.0
                rate = kimura_rate(scoef, HOST_POP_N)
                if rate > 1e-4:
                    candidates.append(( ("ET", v2, s2, scoef, {}), rate ))
                    cum_rate += rate
        else:
            for (bS2, s_ang2, mS2) in self._propose_ER_mutants_host():
                bS_old, s_old_ang, mS_old = self.bS, self.s_angle, self.mS
                self.bS, self.s_angle, self.mS = bS2, s_ang2, mS2
                diag = self._refresh_equilibrium(selector="host")
                f2 = self.host_fit

                self.bS, self.s_angle, self.mS = bS_old, s_old_ang, mS_old
                self._refresh_equilibrium(selector="host")

                scoef = (f2 - current_fit) / max(current_fit, TINY) if current_fit > TINY else 0.0
                rate = kimura_rate(scoef, HOST_POP_N)
                if rate > 1e-4:
                    candidates.append(( ("ER", bS2, s_ang2, mS2, scoef, diag), rate ))
                    cum_rate += rate

        return candidates, cum_rate

    def _evaluate_path_mutations(self) -> Tuple[List[Tuple], float]:
        """Returns (candidate_list, cumulative_rate) for pathogen."""
        candidates = []
        cum_rate = 0.0
        current_fit = self.path_fit

        if not self.es:
            for (v2, s2) in self._propose_ET_mutants_path():
                f2 = path_fitness(v2, s2)
                scoef = (f2 - current_fit) / max(current_fit, TINY) if current_fit > TINY else 0.0
                rate = kimura_rate(scoef, PATH_POP_N)
                if rate > 1e-4:
                    candidates.append(( ("ET", v2, s2, scoef, {}), rate ))
                    cum_rate += rate
        else:
            for (bV2, v_ang2, mV2) in self._propose_ER_mutants_path():
                bV_old, v_old_ang, mV_old = self.bV, self.v_angle, self.mV
                self.bV, self.v_angle, self.mV = bV2, v_ang2, mV2
                diag = self._refresh_equilibrium(selector="path")
                f2 = self.path_fit

                self.bV, self.v_angle, self.mV = bV_old, v_old_ang, mV_old
                self._refresh_equilibrium(selector="path")

                scoef = (f2 - current_fit) / max(current_fit, TINY) if current_fit > TINY else 0.0
                rate = kimura_rate(scoef, PATH_POP_N)
                if rate > 1e-4:
                    candidates.append(( ("ER", bV2, v_ang2, mV2, scoef, diag), rate ))
                    cum_rate += rate

        return candidates, cum_rate

    # ----------------------------------------------------------
    # Gillespie step: evaluate both players, choose proportionally
    # ----------------------------------------------------------

    def step_generation(self) -> Dict:
        """
        Full Gillespie step (Goldstein 2020, Section 5D):
        1. Enumerate all viable host AND pathogen mutations
        2. Compute substitution rates (N * P_fix for each)
        3. Weight by mutation probability (gamma for host, 1-gamma for path)
        4. Choose who mutates proportional to weighted cumulative rate
        5. Choose specific mutation proportional to rate within selected player
        6. Return exponential dwell time = 1 / total_rate
        """

        host_pinned = FIX_HOST_TRAIT is not None
        path_pinned = FIX_PATH_TRAIT is not None

        # Evaluate all mutations for both players
        if host_pinned:
            host_candidates, cum_host_rate = [], 0.0
        else:
            host_candidates, cum_host_rate = self._evaluate_host_mutations()

        if path_pinned:
            path_candidates, cum_path_rate = [], 0.0
        else:
            path_candidates, cum_path_rate = self._evaluate_path_mutations()

        # Weight by mutation probability asymmetry
        weighted_host = prob_host_mutate * cum_host_rate
        weighted_path = (1.0 - prob_host_mutate) * cum_path_rate
        total_rate = weighted_host + weighted_path

        if total_rate <= 0.0:
            self.zero_rate_streak += 1
            return {"mutator": "none", "chosen": None,
                    "omega_host": self.host_fit, "omega_path": self.path_fit,
                    "nash": "", "stable_fp": "", "host_max_flag": "", "path_max_flag": "",
                    "dwell": 1.0,
                    "cum_host_rate": 0.0, "cum_path_rate": 0.0}

        # Exponential dwell time (Gillespie)
        dwell = self.rng.expovariate(total_rate)

        # Choose host or pathogen proportional to weighted rate
        mutate_host = self.rng.random() < (weighted_host / total_rate)
        self.zero_rate_streak = 0

        if mutate_host:
            mutator = "host"
            candidates = host_candidates
            cum_rate = cum_host_rate
        else:
            mutator = "path"
            candidates = path_candidates
            cum_rate = cum_path_rate

        # Choose specific mutation proportional to rate within selected player
        r = self.rng.random() * cum_rate
        acc = 0.0
        chosen_state = candidates[-1][0]  # fallback
        for (state, rate) in candidates:
            acc += rate
            if r <= acc:
                chosen_state = state
                break

        # --- Apply the chosen mutation ---
        if not self.es:
            # ET mode: chosen_state = ("ET", v2, s2, scoef, {})
            _, v_new, s_new, scoef, _ = chosen_state
            self.v, self.s = v_new, s_new
            self.host_fit = host_fitness(self.v, self.s)
            self.path_fit = path_fitness(self.v, self.s)
            return {"mutator": mutator,
                    "chosen": ("ET", scoef, "ET", v_new, s_new),
                    "omega_host": self.host_fit, "omega_path": self.path_fit,
                    "nash": "", "stable_fp": "", "host_max_flag": "", "path_max_flag": "",
                    "dwell": dwell,
                    "cum_host_rate": cum_host_rate, "cum_path_rate": cum_path_rate}
        else:
            # ER mode
            if mutator == "host":
                _, bS2, s_ang2, mS2, scoef, _ = chosen_state
                self.bS, self.s_angle, self.mS = bS2, s_ang2, mS2
                diag = self._refresh_equilibrium(selector="host", track=True)
            else:
                _, bV2, v_ang2, mV2, scoef, _ = chosen_state
                self.bV, self.v_angle, self.mV = bV2, v_ang2, mV2
                diag = self._refresh_equilibrium(selector="path", track=True)

            return {"mutator": mutator,
                    "chosen": ("ER", scoef, "ER", self.v, self.s),
                    "omega_host": self.host_fit, "omega_path": self.path_fit,
                    "nash": diag["nash"], "stable_fp": diag["stable_fp"],
                    "host_max_flag": diag["host_max_flag"], "path_max_flag": diag["path_max_flag"],
                    "dwell": dwell,
                    "cum_host_rate": cum_host_rate, "cum_path_rate": cum_path_rate}


# ============================================================
# CSV writer
# ============================================================

def _runtime_write_every(gen: int) -> bool:
    if gen < 0:
        return False
    if RUNTIME_MODE == "full":
        return (gen % max(1, int(write_every)) == 0)
    return (gen % max(1, int(WRITE_EVERY_FAST)) == 0)

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
            if gen == 0:
                t = 0.0

            if gen % PROGRESS_EVERY == 0:
                mode = "ES" if sim.es else "EI"
                print(f"[{mode}] gen={gen} t={t:.3e} v={sim.v:.3f} s={sim.s:.3f} zr={sim.zero_rate_streak}")

            record = _runtime_write_every(gen)

            if record:
                w.writerow([
                    "pre", gen, f"{t:.12e}", "NA",
                    f"{sim.v:.6f}", f"{sim.bV:.6f}", f"{sim.v_angle:.6f}", f"{sim.mV:.6f}",
                    f"{sim.s:.6f}", f"{sim.bS:.6f}", f"{sim.s_angle:.6f}", f"{sim.mS:.6f}",
                    f"{sim.path_fit:.6f}", f"{sim.host_fit:.6f}",
                    "", "", "", "", "", "",
                    "", "", ""
                ])
                f.flush()

            result = sim.step_generation()
            dwell_out = max(float(result.get("dwell", 1.0)), DWELL_MIN)
            t += dwell_out

            if record:
                mutSelCoeff = ""
                mutClass = ""
                if result.get("chosen") is not None:
                    _, scoef, cls, _, _ = result["chosen"]
                    mutSelCoeff = f"{scoef:.3e}"
                    mutClass = cls

                w.writerow([
                    "post", gen, f"{t:.12e}", result.get("mutator","NA"),
                    f"{sim.v:.6f}", f"{sim.bV:.6f}", f"{sim.v_angle:.6f}", f"{sim.mV:.6f}",
                    f"{sim.s:.6f}", f"{sim.bS:.6f}", f"{sim.s_angle:.6f}", f"{sim.mS:.6f}",
                    f"{sim.path_fit:.6f}", f"{sim.host_fit:.6f}",
                    f"{result.get('omega_path','')}", f"{result.get('omega_host','')}",
                    result.get("nash",""), result.get("stable_fp",""),
                    result.get("host_max_flag",""), result.get("path_max_flag",""),
                    mutSelCoeff, mutClass,
                    f"{dwell_out:.6e}"
                ])
