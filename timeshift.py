"""
timeshift.py — time-shift (reciprocal cross-inoculation) assay on stored
GoldsteinGameTheory trajectories.

Design follows Gaba & Ebert (2009) TREE 24:226-232: challenge one antagonist
from time t against the other from time t + delta, for delta spanning past,
contemporary and future, and contrast same-lineage (sympatric) against
different-lineage (allopatric) pairings.

No re-simulation is needed: fitness is a deterministic function of the realised
(v, c), so the assay is a re-evaluation pass over trajectories already on disk.

TWO PROTOCOLS
-------------
"phenotype"  Cross the stored realised traits directly: W(v_path(t1), c_host(t2)).
             Empirical analogue: an assay in which induction is blocked or the
             phenotype is otherwise fixed (constitutive mutant, heat-killed
             challenge, standardised non-inducing conditions).

"rule"       Cross the reaction-norm parameters (bS, mS) x (bV, mV), re-settle
             the mutual fixed point, and evaluate fitness at the NEWLY realised
             (v*, c*). Empirical analogue: standard live reciprocal cross-
             infection, in which a conditional strategy re-expresses against
             whichever partner it now faces.

For EThost_ETpath the two protocols coincide identically (both slopes are zero,
so the rule IS the phenotype). For ERhost_ETpath and EThost_ERpath they differ
only through the one reactive player. For ERhost_ERpath they differ most, and
the gap between them is the plasticity contribution.

DATA CAVEATS THAT THIS MODULE HANDLES (established by inspecting results/)
-------------------------------------------------------------------------
1. In EThost_ETpath runs the simulation sets evolved_strategy=False and never
   calls _refresh_equilibrium, so the bV/bS/mV/mS columns keep their
   initialisation values for the whole run and do NOT describe the realised
   traits (bV is literally constant, and clamp(bV + mV*s) misses v by up to
   0.35). Rules for those runs are therefore reconstructed from the realised
   traits as (b = trait, m = 0), which is exact for a non-reactive player.
   In every evolved_strategy=True run the stored columns are consistent with
   the realised traits to <= 6e-4 and are used as-is. `validate_rules()`
   reports the residual per run.

2. Elapsed continuous time is not comparable across conditions: in the
   minimal model ER/ER runs cover t ~ 0.16 while ET/ET runs cover t ~ 1500
   over the same 100k substitutions, because ER dwell times are ~1e4x shorter. Delta on the
   continuous-time axis must therefore be interpreted within a condition, not
   across conditions. Both axes are supported ("sub" and "time").

3. The main runs record one row per 10 substitutions. They support the
   Delta=0 sympatric-vs-allopatric contrast (which needs no lag resolution)
   and the long-lag tail, but the shape of the profile at short lag needs the
   write_every=1 zoom runs (results_zoom/). `autocorr()` is provided so this is checked rather than assumed.

Depends on simulation.py for the fitness registry and find_all_equilibria.
Pure standard library, matching the rest of the repo.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import glob
import json
import math
import os
import random
import sys
from array import array
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulation as sim


# ==========================================================================
# Loading
# ==========================================================================

NUMERIC_COLS = ("gen", "time", "v", "bV", "mV", "s", "bS", "mS",
                "pathFit", "hostFit")


@dataclass
class Run:
    """One stored trajectory plus the metadata needed to interpret it."""
    path: str
    cfg: dict
    gen: List[float]
    time: List[float]
    v: List[float]          # realised virulence
    c: List[float]          # realised clearance (stored as 's')
    bV: List[float]         # effective pathogen rule intercept
    mV: List[float]         # effective pathogen rule slope
    bS: List[float]         # effective host rule intercept
    mS: List[float]         # effective host rule slope
    mutator: List[str]      # which player substituted into this state
    rule_residual: float    # max |clamp(b + m*partner) - trait| before repair
    repaired: bool          # True if rules were reconstructed from traits

    @property
    def condition(self) -> str:
        return self.cfg.get("condition", "")

    @property
    def fitness_model(self) -> str:
        return self.cfg.get("fitness_model", "")

    @property
    def rep(self) -> int:
        return int(self.cfg.get("rep") or 0)

    @property
    def tracking_k(self) -> float:
        return float(self.cfg.get("TRACKING_K") or 0.0)

    def __len__(self) -> int:
        return len(self.gen)


def load_run(run_dir: str, min_gen: float = 0.0) -> Run:
    """Load the post-substitution rows of one run directory.

    Only 'post' rows are used: a 'pre' row is the state before the candidate
    substitution and carries no fitness or dwell information.
    """
    csv_path = os.path.join(run_dir, "simulation.csv")
    with open(os.path.join(run_dir, "config.json")) as fh:
        cfg = json.load(fh)

    cols: Dict[str, List[float]] = {k: [] for k in NUMERIC_COLS}
    mutator: List[str] = []
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            if row.get("event") != "post":
                continue
            if float(row["gen"]) < min_gen:
                continue
            for k in NUMERIC_COLS:
                raw = row.get(k, "")
                cols[k].append(float(raw) if raw not in ("", "NA") else float("nan"))
            mutator.append(row.get("mutator", ""))

    v, c = cols["v"], cols["s"]
    bV, mV, bS, mS = cols["bV"], cols["mV"], cols["bS"], cols["mS"]

    # How well do the stored rules reproduce the realised traits?
    resid = 0.0
    for i in range(len(v)):
        resid = max(resid,
                    abs(sim.clamp01(bV[i] + mV[i] * c[i]) - v[i]),
                    abs(sim.clamp01(bS[i] + mS[i] * v[i]) - c[i]))

    # evolved_strategy=False runs never refresh the rule columns, so they are
    # stale initialisation values. Both players are non-reactive there, so the
    # effective rule is exactly (b = trait, m = 0).
    repaired = not bool(cfg.get("evolved_strategy", True))
    if repaired:
        bV, mV = list(v), [0.0] * len(v)
        bS, mS = list(c), [0.0] * len(c)

    return Run(path=run_dir, cfg=cfg, gen=cols["gen"], time=cols["time"],
               v=v, c=c, bV=bV, mV=mV, bS=bS, mS=mS, mutator=mutator,
               rule_residual=resid, repaired=repaired)


def discover_runs(roots: Sequence[str] = ("results",),
                  fitness_model: Optional[str] = None,
                  condition: Optional[str] = None) -> List[str]:
    """Find run directories (those containing both config.json and simulation.csv)."""
    found: List[str] = []
    for root in roots:
        for cfg_path in glob.glob(os.path.join(root, "**", "config.json"),
                                  recursive=True):
            d = os.path.dirname(cfg_path)
            if not os.path.exists(os.path.join(d, "simulation.csv")):
                continue
            try:
                with open(cfg_path) as fh:
                    cfg = json.load(fh)
            except (OSError, json.JSONDecodeError):
                continue
            if fitness_model and cfg.get("fitness_model") != fitness_model:
                continue
            if condition and cfg.get("condition") != condition:
                continue
            found.append(d)
    return sorted(found)


def activate_model(run: Run) -> None:
    """Point simulation.py's module-level fitness functions at this run's model."""
    sim.set_fitness_model(run.fitness_model)
    if run.fitness_model == "tracking":
        sim.TRACKING_K = run.tracking_k


# ==========================================================================
# Settling a cross-pairing
# ==========================================================================

@dataclass(frozen=True)
class Settled:
    v: float
    c: float
    status: str          # "fixed" | "cycle2" | "noconv"
    iters: int
    cycle_amp: float     # 0 for a fixed point; peak-to-peak in c for a 2-cycle
    slope_product: float


def settle_pair(bS: float, mS: float, bV: float, mV: float,
                c_start: Optional[float] = None,
                tol: float = 1e-13, max_iter: int = 5000) -> Settled:
    """Let a cross-paired host rule and pathogen rule settle against each other.

    The composed map is  c -> clamp(bS + mS * clamp(bV + mV * c)),  with
    derivative mS*mV away from the clamps, so it contracts iff |mS*mV| < 1.

    Why iterate rather than call _best_equilibrium_for_player: in the
    simulation, the equilibrium is chosen as the best fixed point for whichever
    player just mutated. In a cross-inoculation nobody has just mutated, so
    that convention does not transfer. Iterating instead reproduces what the
    assay physically does — put the two together and let induction settle.

    Default start is the naive (uninduced) host, c = clamp(bS), i.e. the
    constitutive level expressed before any pathogen is encountered. This
    matters: in ERhost_ERpath, 9.6% of cross-pairings settle to a different
    v* depending on the starting point, and where multiple fixed points exist
    they typically span the whole trait range. Pass c_start explicitly to test
    sensitivity to this convention.

    When |mS*mV| < -1 the iteration can land on a 2-cycle rather than a fixed
    point (17% of ERhost_ERpath cross-pairings). That is a substantive outcome
    — induction that never settles — not a numerical failure. We report the
    cycle mean, which is what a real assay integrating over the oscillation
    would measure, and flag it via `status` and `cycle_amp` so those pairings
    can be filtered or analysed separately.
    """
    clamp = sim.clamp01
    c = clamp(bS if c_start is None else c_start)
    prev: Optional[float] = None

    for i in range(max_iter):
        c_next = clamp(bS + mS * clamp(bV + mV * c))
        if abs(c_next - c) < tol:
            return Settled(v=clamp(bV + mV * c_next), c=c_next, status="fixed",
                           iters=i + 1, cycle_amp=0.0, slope_product=mS * mV)
        if prev is not None and abs(c_next - prev) < tol:
            # 2-cycle between c and c_next: report the mean of the two states
            c_a, c_b = c, c_next
            v_a, v_b = clamp(bV + mV * c_a), clamp(bV + mV * c_b)
            return Settled(v=0.5 * (v_a + v_b), c=0.5 * (c_a + c_b),
                           status="cycle2", iters=i + 1,
                           cycle_amp=abs(c_b - c_a), slope_product=mS * mV)
        prev = c
        c = c_next

    return Settled(v=clamp(bV + mV * c), c=c, status="noconv",
                   iters=max_iter, cycle_amp=float("nan"),
                   slope_product=mS * mV)


def solve_pair(bS: float, mS: float, bV: float, mV: float,
               selector: str) -> Settled:
    """Resolve a cross-pairing the way the simulation does: find every fixed
    point of the composed clamped map and take the one best for `selector`.

    This is `_refresh_equilibrium`'s convention. It is offered alongside
    `settle_pair` because the two disagree substantially in ERhost_ERpath, and
    the disagreement is a modelling choice rather than a numerical detail:

      - 32.5% of the *realised* ERhost_ERpath states in results/minimal have
        |mS*mV| > 1, i.e. the simulation is sitting on a fixed point that
        repels. Such a point is a genuine root of the algebraic fixed-point
        equation but is not reachable by any dynamical induction process.
      - Re-settling a host against its OWN contemporaneous pathogen (delta=0,
        sympatric) by iteration changes the realised trait by more than 0.1 in
        20% of records, and by up to the full trait range.

    Under an algebraic/ESS reading of the reaction norms the solver convention
    is the right one, and using it makes the delta=0 point of the assay exactly
    the state the simulation recorded. Under a dynamical reading — plasticity
    settling within a host during an infection — `settle_pair` is right and
    repelling roots should not be expressible. Both are provided; the choice
    should be stated in the methods rather than left implicit.

    All of this is confined to ERhost_ERpath. In the other three conditions at
    least one slope is ~0, so |mS*mV| ~ 0, the root is unique, and every mode
    agrees to solver tolerance (<= 1e-4).
    """
    eqs = sim.find_all_equilibria(bS, mS, bV, mV)
    best = sim._best_equilibrium_for_player(eqs, selector)
    if best is None:
        return Settled(v=sim.clamp01(bV), c=sim.clamp01(bS), status="noroot",
                       iters=0, cycle_amp=float("nan"), slope_product=mS * mV)
    return Settled(v=best.v, c=best.s,
                   status="fixed" if best.stable else "repelling",
                   iters=0, cycle_amp=0.0, slope_product=mS * mV)


def resolve(bS: float, mS: float, bV: float, mV: float,
            mode: str, c_contemporary: Optional[float] = None) -> Settled:
    """Dispatch to the requested equilibrium-resolution convention."""
    if mode == "naive":
        return settle_pair(bS, mS, bV, mV, c_start=None)
    if mode == "contemporary":
        return settle_pair(bS, mS, bV, mV, c_start=c_contemporary)
    if mode == "solver-host":
        return solve_pair(bS, mS, bV, mV, "host")
    if mode == "solver-path":
        return solve_pair(bS, mS, bV, mV, "path")
    raise ValueError(f"unknown settling mode {mode!r}")


SETTLE_MODES = ("naive", "contemporary", "solver-host", "solver-path")


def n_equilibria(bS: float, mS: float, bV: float, mV: float) -> int:
    """How many fixed points this pairing admits (diagnostic; expensive)."""
    try:
        return len(sim.find_all_equilibria(bS, mS, bV, mV))
    except Exception:
        return -1


# ==========================================================================
# Indexing on the two time axes
# ==========================================================================

def locf_index(times: Sequence[float], targets: Sequence[float]) -> List[int]:
    """Last-observation-carried-forward indices for a continuous-time grid.

    Traits are piecewise constant between substitutions, so interpolation onto
    a time grid is a step lookup, not a linear one. Recorded rows are every
    write_every-th substitution while time accumulates continuously through
    variable dwell times, so recorded points are not evenly spaced in t.
    """
    out = []
    n = len(times)
    for t in targets:
        i = bisect.bisect_right(times, t) - 1
        out.append(min(max(i, 0), n - 1))
    return out


# ==========================================================================
# The assay
# ==========================================================================

OUT_FIELDS = [
    "fitness_model", "condition", "tracking_k", "protocol", "axis", "settle",
    "rep_path", "rep_host", "sympatric",
    "t_focal", "delta", "t_host",
    "v_stored", "c_stored", "v_realized", "c_realized",
    "W_P", "W_H", "on_boundary", "status", "cycle_amp", "slope_product", "n_roots",
]

BOUNDARY_EPS = 1e-6


def assay(path_run: Run, host_run: Run,
          deltas: Sequence[float],
          focal: Sequence[float],
          protocol: str,
          axis: str = "sub",
          count_roots: bool = False,
          settle: str = "naive") -> List[dict]:
    """Cross one pathogen lineage against one host lineage over a delta grid.

    axis="sub"   focal and deltas are indices into the recorded series
                 (1 unit = write_every substitutions).
    axis="time"  focal and deltas are in continuous time units.

    rep_path == rep_host is the sympatric series; unequal reps give the
    allopatric control that separates lineage-specific matching from shared
    directional trends in trait level.
    """
    if protocol not in ("phenotype", "rule"):
        raise ValueError(f"unknown protocol {protocol!r}")
    if axis not in ("sub", "time"):
        raise ValueError(f"unknown axis {axis!r}")

    if axis == "sub":
        p_idx = [int(t) for t in focal]
    else:
        p_idx = locf_index(path_run.time, focal)

    rows: List[dict] = []
    n_host = len(host_run)

    for d in deltas:
        if axis == "sub":
            h_idx = [min(max(int(t + d), 0), n_host - 1) for t in focal]
            t_host = [host_run.gen[i] for i in h_idx]
        else:
            targets = [t + d for t in focal]
            h_idx = locf_index(host_run.time, targets)
            t_host = targets

        for k in range(len(focal)):
            pi, hi = p_idx[k], h_idx[k]
            v_stored, c_stored = path_run.v[pi], host_run.c[hi]

            if protocol == "phenotype":
                v, c = v_stored, c_stored
                status, amp, sp = "fixed", 0.0, float("nan")
            else:
                st = resolve(host_run.bS[hi], host_run.mS[hi],
                             path_run.bV[pi], path_run.mV[pi],
                             mode=settle, c_contemporary=c_stored)
                v, c = st.v, st.c
                status, amp, sp = st.status, st.cycle_amp, st.slope_product

            nr = -1
            if count_roots and protocol == "rule":
                nr = n_equilibria(host_run.bS[hi], host_run.mS[hi],
                                  path_run.bV[pi], path_run.mV[pi])

            rows.append({
                "fitness_model": path_run.fitness_model,
                "condition": path_run.condition,
                "tracking_k": path_run.tracking_k,
                "protocol": protocol,
                "axis": axis,
                "settle": settle if protocol == "rule" else "",
                "rep_path": path_run.rep,
                "rep_host": host_run.rep,
                "sympatric": int(path_run.rep == host_run.rep),
                "t_focal": focal[k],
                "delta": d,
                "t_host": t_host[k],
                "v_stored": v_stored,
                "c_stored": c_stored,
                "v_realized": v,
                "c_realized": c,
                "W_P": sim.path_fitness(v, c),
                "W_H": sim.host_fitness(v, c),
                # Mismatched rules settle onto a clamped (boundary) equilibrium
                # far more often than coevolved ones, and W_P = v(1-v)(1-c) is
                # exactly 0 at v in {0,1}. That alone moves the allopatric mean,
                # so boundary hits must be separable from any matching claim.
                "on_boundary": int(v <= BOUNDARY_EPS or v >= 1.0 - BOUNDARY_EPS or
                                   c <= BOUNDARY_EPS or c >= 1.0 - BOUNDARY_EPS),
                "status": status,
                "cycle_amp": amp,
                "slope_product": sp,
                "n_roots": nr,
            })
    return rows


def make_focal_grid(runs: Sequence[Run], deltas: Sequence[float],
                    n_focal: int, axis: str) -> List[float]:
    """Focal points valid for every delta, so the profile is not confounded.

    If the focal set changed with delta, differences across delta would mix the
    shift with a change in which moments of the run were sampled. Restricting
    to focal points where t+delta is in range for all deltas keeps the sampled
    set identical across the whole profile.
    """
    d_lo, d_hi = min(deltas), max(deltas)
    if axis == "sub":
        n = min(len(r) for r in runs)
        lo = max(0, int(math.ceil(-d_lo)))
        hi = min(n - 1, n - 1 - int(math.ceil(d_hi)))
    else:
        t_end = min(r.time[-1] for r in runs)
        lo = max(0.0, -d_lo)
        hi = t_end - max(0.0, d_hi)
    if hi <= lo:
        raise ValueError(f"delta range [{d_lo}, {d_hi}] leaves no valid focal "
                         f"window on axis={axis!r}")
    if n_focal == 1:
        return [0.5 * (lo + hi)]
    step = (hi - lo) / (n_focal - 1)
    return [lo + step * i for i in range(n_focal)]


# ==========================================================================
# Diagnostics
# ==========================================================================

def autocorr(x: Sequence[float], lags: Sequence[int]) -> List[float]:
    """Sample autocorrelation of a series at the given lags."""
    n = len(x)
    mu = sum(x) / n
    d = [xi - mu for xi in x]
    den = sum(di * di for di in d)
    if den <= 0:
        return [float("nan")] * len(lags)
    return [sum(d[i] * d[i + k] for i in range(n - k)) / den if k < n else float("nan")
            for k in lags]


def decorrelation_lag(x: Sequence[float], max_lag: int = 500) -> int:
    """First lag at which the autocorrelation drops below 1/e.

    Sets both the meaningful delta resolution and the block size for a block
    bootstrap over time blocks (time points within a run are autocorrelated,
    so effective n is far below the raw pair count).
    """
    a = autocorr(x, range(min(max_lag, len(x) - 1) + 1))
    for k, val in enumerate(a):
        if val < 1.0 / math.e:
            return k
    return len(a)


def validate_rules(run: Run) -> dict:
    """Per-run report of rule/trait consistency and mixing timescale."""
    return {
        "path": run.path,
        "condition": run.condition,
        "fitness_model": run.fitness_model,
        "rep": run.rep,
        "n_records": len(run),
        "evolved_strategy": bool(run.cfg.get("evolved_strategy", True)),
        "rules_repaired": int(run.repaired),
        "rule_residual": run.rule_residual,
        "t_max": run.time[-1] if run.time else float("nan"),
        "tau_v": decorrelation_lag(run.v),
        "tau_c": decorrelation_lag(run.c),
    }


SUMMARY_FIELDS = [
    "fitness_model", "condition", "tracking_k", "protocol", "settle", "axis",
    "sympatric", "delta", "n", "n_pairings",
    "W_P_mean", "W_P_se", "W_H_mean", "W_H_se",
    "frac_boundary", "frac_cycle2",
    "W_P_mean_interior", "E_f", "E_g", "cov_fg",
    "n_lineages", "W_P_se_lineage", "W_H_se_lineage",
]

PAIR_FIELDS = [
    "fitness_model", "condition", "tracking_k", "protocol", "settle", "axis",
    "sympatric", "delta", "rep_path", "rep_host", "n", "W_P_mean", "W_H_mean",
]


class CellAccumulator:
    """Streaming aggregation of assay rows into one record per assay cell.

    A cell is one (model, condition, k, protocol, settle, axis, sympatric,
    delta) combination. Rows are folded in as they are produced rather than
    held until a group finishes: a row dict costs ~1 KB in CPython, and the
    allopatric zoom design (64 pairings x 1500 focal x 81 deltas x 2 protocols)
    is ~15M rows per group, i.e. ~15 GB. Here a row costs 16 bytes (W_P and W_H
    as doubles, kept in arrival order for the block bootstrap) plus O(1)
    running sums, and the summary comes out the same.

    Two standard errors are reported per cell:

      *_se          moving block bootstrap over focal-time blocks. Treats
                    autocorrelated time points as the replication unit, so it
                    measures how precisely this design pins the mean, not how
                    well the mean generalises to new lineages. `n_block` is in
                    focal-grid steps and should be at least the decorrelation lag.
      *_se_lineage  the pathogen lineage (rep_path) as the unit: SD of
                    per-lineage cell means / sqrt(n_lineages). A sympatric cell
                    has one pairing per lineage; an allopatric lineage mean
                    averages over host partners that are shared across lineages,
                    so this is a mild underestimate. Blank below two lineages.

    Neither is the right SE for a sympatric-minus-allopatric contrast: a
    lineage's sympatric and allopatric means are positively correlated, so
    combining the two cell SEs in quadrature overstates the contrast SE. Take
    pair_rows() and difference within lineage instead.
    """

    def __init__(self) -> None:
        self.cells: Dict[Tuple, dict] = {}

    def add(self, rows: Iterable[dict]) -> None:
        for r in rows:
            key = (r["fitness_model"], r["condition"], r["tracking_k"],
                   r["protocol"], r["settle"], r["axis"],
                   r["sympatric"], r["delta"])
            c = self.cells.get(key)
            if c is None:
                c = self.cells[key] = {
                    "wp": array("d"), "wh": array("d"),
                    "n_bdry": 0, "n_cyc": 0, "int_sum": 0.0, "int_n": 0,
                    "sf": 0.0, "sg": 0.0, "sfg": 0.0, "pairs": {},
                }
            wp, wh = r["W_P"], r["W_H"]
            c["wp"].append(wp)
            c["wh"].append(wh)
            if r["on_boundary"]:
                c["n_bdry"] += 1
            else:
                c["int_sum"] += wp
                c["int_n"] += 1
            if r["status"] == "cycle2":
                c["n_cyc"] += 1
            if key[0] == "minimal":
                # decompose_minimal() as running sums: f(v) = v(1-v), g(c) = 1-c
                v = r["v_realized"]
                f, g = v * (1.0 - v), 1.0 - r["c_realized"]
                c["sf"] += f
                c["sg"] += g
                c["sfg"] += f * g
            p = c["pairs"].setdefault((r["rep_path"], r["rep_host"]), [0.0, 0.0, 0])
            p[0] += wp
            p[1] += wh
            p[2] += 1

    def summary(self, n_block: int = 20) -> List[dict]:
        out: List[dict] = []
        for key, c in sorted(self.cells.items(), key=lambda kv: str(kv[0])):
            n = len(c["wp"])
            n_lin, wp_se_lin = _lineage_se(c["pairs"], 0)
            _, wh_se_lin = _lineage_se(c["pairs"], 1)
            rec = {
                "fitness_model": key[0], "condition": key[1], "tracking_k": key[2],
                "protocol": key[3], "settle": key[4], "axis": key[5],
                "sympatric": key[6], "delta": key[7],
                "n": n,
                "n_pairings": len(c["pairs"]),
                "W_P_mean": sum(c["wp"]) / n,
                "W_P_se": _block_se(c["wp"], n_block),
                "W_H_mean": sum(c["wh"]) / n,
                "W_H_se": _block_se(c["wh"], n_block),
                "frac_boundary": c["n_bdry"] / n,
                "frac_cycle2": c["n_cyc"] / n,
                "W_P_mean_interior": c["int_sum"] / c["int_n"] if c["int_n"] else "",
                "E_f": "", "E_g": "", "cov_fg": "",
                "n_lineages": n_lin,
                "W_P_se_lineage": wp_se_lin,
                "W_H_se_lineage": wh_se_lin,
            }
            if key[0] == "minimal":
                ef, eg = c["sf"] / n, c["sg"] / n
                rec["E_f"], rec["E_g"] = ef, eg
                rec["cov_fg"] = c["sfg"] / n - ef * eg
            out.append(rec)
        return out

    def pair_rows(self) -> List[dict]:
        """Per-(cell, pathogen lineage, host lineage) means, for paired contrasts."""
        out: List[dict] = []
        for key, c in sorted(self.cells.items(), key=lambda kv: str(kv[0])):
            for (rp, rh), (swp, swh, n) in sorted(c["pairs"].items()):
                out.append({
                    "fitness_model": key[0], "condition": key[1], "tracking_k": key[2],
                    "protocol": key[3], "settle": key[4], "axis": key[5],
                    "sympatric": key[6], "delta": key[7],
                    "rep_path": rp, "rep_host": rh, "n": n,
                    "W_P_mean": swp / n, "W_H_mean": swh / n,
                })
        return out


def _lineage_se(pairs: Dict[Tuple[int, int], list], col: int) -> Tuple[int, object]:
    """SE of a cell mean with the pathogen lineage as the unit of replication.

    `col` picks the running sum (0 = W_P, 1 = W_H). Returns (n_lineages, SE).
    The SE is written blank rather than NaN below two lineages so the CSV
    column still parses as numeric in R.
    """
    by_path: Dict[int, List[float]] = {}
    for (rp, _rh), sums in pairs.items():
        acc = by_path.setdefault(rp, [0.0, 0])
        acc[0] += sums[col]
        acc[1] += sums[2]
    means = [s / k for s, k in by_path.values() if k]
    m = len(means)
    if m < 2:
        return m, ""
    mu = sum(means) / m
    sd = math.sqrt(sum((x - mu) ** 2 for x in means) / (m - 1))
    return m, sd / math.sqrt(m)


def summarize(rows: Sequence[dict], n_block: int = 20) -> List[dict]:
    """Collapse pairing-level rows already in memory to one record per cell.

    Thin wrapper over CellAccumulator, which the CLI uses directly so that rows
    never accumulate. See CellAccumulator for the columns and the two SEs.
    """
    acc = CellAccumulator()
    acc.add(rows)
    return acc.summary(n_block=n_block)


def _block_se(x: Sequence[float], n_block: int, n_boot: int = 200,
              seed: int = 12345) -> float:
    """Moving block bootstrap standard error of the mean."""
    n = len(x)
    if n < 2:
        return float("nan")
    b = max(1, min(n_block, n))
    n_draw = max(1, n // b)
    rng = random.Random(seed)
    means = []
    for _ in range(n_boot):
        tot = 0.0
        for _ in range(n_draw):
            start = rng.randrange(0, n - b + 1)
            tot += sum(x[start:start + b])
        means.append(tot / (n_draw * b))
    mu = sum(means) / n_boot
    return math.sqrt(sum((m - mu) ** 2 for m in means) / (n_boot - 1))


def decompose_minimal(rows: Iterable[dict]) -> dict:
    """Split mean pathogen fitness into trait-level and matching components.

    Only valid for the minimal model, where W_P = v(1-v) * (1-c) is a product
    of a pathogen-only term f(v) and a host-only term g(c). Then exactly

        E[W_P] = E[f] * E[g] + Cov(f, g)

    E[f] is where the pathogen's trait sits, E[g] is where the host's trait
    sits, and Cov(f, g) is the ONLY term that reflects the two players being
    matched to each other. Comparing this decomposition between sympatric and
    allopatric sets (or across delta) separates "reciprocal adaptation in trait
    levels" from genotype-specific local adaptation.

    This matters because the raw contrast is misleading here: in results/minimal
    the sympatric-minus-allopatric advantage in W_P is +0.018 to +0.026, which
    looks like textbook local adaptation, but the Cov term contributes
    -0.003 to +0.0001 of it. The rest is trait positioning. In ERhost_ETpath
    the pathogen is non-reactive, so f(v) is constant and Cov is identically
    zero, yet the contrast is still +0.026 — which is the cleanest available
    demonstration that the contrast does not require matching at all.
    """
    F, G, W = [], [], []
    for r in rows:
        v, c = float(r["v_realized"]), float(r["c_realized"])
        F.append(v * (1.0 - v))
        G.append(1.0 - c)
        W.append(float(r["W_P"]))
    n = len(F)
    if n == 0:
        return {}
    ef = sum(F) / n
    eg = sum(G) / n
    cov = sum((f - ef) * (g - eg) for f, g in zip(F, G)) / n
    return {"n": n, "E_f": ef, "E_g": eg, "cov": cov,
            "E_WP_reconstructed": ef * eg + cov,
            "E_WP_direct": sum(W) / n}


def reproduce_recorded_state(run: Run) -> dict:
    """End-to-end pipeline check: re-derive each recorded state from its own rules.

    Applies the simulation's exact convention — all fixed points of the composed
    clamped map, keeping the one best for the player named in the `mutator`
    column — to the (bS, mS, bV, mV) we extracted, and compares the result with
    the recorded (v, c). Residuals at solver tolerance mean the rule extraction,
    the ET-run repair and the indexing are all correct. This is a check on the
    pipeline, not a result: any large residual here is a bug, not a finding.

    Returns max/median absolute deviation and the fraction exceeding 1e-6.
    """
    devs = []
    for i in range(len(run)):
        who = run.mutator[i] if run.mutator[i] in ("host", "path") else "path"
        st = solve_pair(run.bS[i], run.mS[i], run.bV[i], run.mV[i], who)
        devs.append(max(abs(st.v - run.v[i]), abs(st.c - run.c[i])))
    devs.sort()
    n = len(devs)
    return {
        "path": run.path,
        "condition": run.condition,
        "rep": run.rep,
        "n": n,
        "median_dev": devs[n // 2] if n else float("nan"),
        "max_dev": devs[-1] if n else float("nan"),
        "frac_gt_1e6": sum(1 for d in devs if d > 1e-6) / n if n else float("nan"),
    }


def check_phenotype_profile(path_run: Run, host_run: Run,
                            deltas: Sequence[float], focal: Sequence[float],
                            axis: str = "sub") -> List[Tuple[float, float]]:
    """Closed-form check on the phenotype protocol under the minimal model.

    W_P = v(1-v)(1-c) factorises as f(v)*g(c), so the whole phenotype-shift
    profile must reduce to the lagged cross-moment E[f(v_t) g(c_{t+delta})].
    Disagreement with the assay output means a bug in the indexing or pairing,
    not a result. Meaningless for tracking (k>0), which is non-separable, and
    for acute/chronic/taylor, which are not products of a v-term and a c-term.
    """
    if path_run.fitness_model != "minimal":
        raise ValueError("closed form applies to the minimal model only")
    f = lambda v: v * (1.0 - v)
    g = lambda c: 1.0 - c
    p_idx = ([int(t) for t in focal] if axis == "sub"
             else locf_index(path_run.time, focal))
    fv = [f(path_run.v[i]) for i in p_idx]
    n_host = len(host_run)
    out = []
    for d in deltas:
        if axis == "sub":
            h_idx = [min(max(int(t + d), 0), n_host - 1) for t in focal]
        else:
            h_idx = locf_index(host_run.time, [t + d for t in focal])
        gc = [g(host_run.c[i]) for i in h_idx]
        out.append((d, sum(a * b for a, b in zip(fv, gc)) / len(fv)))
    return out


# ==========================================================================
# CLI
# ==========================================================================

def _parse_deltas(spec: str) -> List[float]:
    """'-10:10:1' -> range; '0,1,5' -> explicit list."""
    if ":" in spec:
        lo, hi, step = (float(x) for x in spec.split(":"))
        n = int(round((hi - lo) / step))
        return [lo + step * i for i in range(n + 1)]
    return [float(x) for x in spec.split(",")]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Time-shift (reciprocal cross-inoculation) assay on stored trajectories.")
    ap.add_argument("--roots", nargs="+", default=["results"],
                    help="Directories to scan for runs (default: results).")
    ap.add_argument("-f", "--fitness", default=None, help="Filter by fitness model.")
    ap.add_argument("-c", "--condition", default=None, help="Filter by condition.")
    ap.add_argument("--deltas", default="-20:20:1",
                    help="Delta grid, 'lo:hi:step' or comma list. Units follow --axis.")
    ap.add_argument("--axis", choices=["sub", "time"], default="sub",
                    help="'sub': delta in recorded steps. 'time': delta in continuous time.")
    ap.add_argument("--n-focal", type=int, default=200,
                    help="Number of focal time points per pairing.")
    ap.add_argument("--protocols", nargs="+", default=["rule", "phenotype"],
                    choices=["rule", "phenotype"])
    ap.add_argument("--allopatric", action="store_true",
                    help="Include cross-lineage pairings (all rep pairs, not just matched).")
    ap.add_argument("--settle", nargs="+", default=["naive"], choices=list(SETTLE_MODES),
                    help="Equilibrium-resolution convention(s) for the rule protocol. "
                         "'naive'/'contemporary' iterate induction to a settled state, "
                         "differing only in the starting clearance; 'solver-host'/"
                         "'solver-path' reproduce the simulation's own convention and "
                         "so make delta=0 exactly the recorded state. Pass several to "
                         "emit them side by side. They agree except in ERhost_ERpath.")
    ap.add_argument("--count-roots", action="store_true",
                    help="Also count fixed points per pairing (slow; diagnostic).")
    ap.add_argument("--raw", default=None,
                    help="Also write the un-aggregated pairing-level rows here. "
                         "Large: the full main-run design is ~1 GB.")
    ap.add_argument("--n-block", type=int, default=20,
                    help="Block length (in focal-grid steps) for the bootstrap SE.")
    ap.add_argument("--pairs", default=None,
                    help="Also write per-(cell, rep_path, rep_host) means here. Small; "
                         "needed for lineage-paired sym - allo contrasts and their SEs.")
    ap.add_argument("--min-gen", type=float, default=0.0,
                    help="Drop records below this generation (extra burn-in).")
    ap.add_argument("-o", "--out", default="timeshift.csv")
    ap.add_argument("--diagnostics", default=None,
                    help="Optional path for the per-run validation table.")
    args = ap.parse_args(argv)

    dirs = discover_runs(args.roots, args.fitness, args.condition)
    if not dirs:
        print("No runs matched.", file=sys.stderr)
        return 1

    runs = [load_run(d, min_gen=args.min_gen) for d in dirs]

    # Group by everything that must match for a pairing to be meaningful.
    groups: Dict[Tuple, List[Run]] = {}
    for r in runs:
        groups.setdefault((r.fitness_model, r.condition, r.tracking_k), []).append(r)

    if args.diagnostics:
        with open(args.diagnostics, "w", newline="") as fh:
            w = None
            for r in runs:
                d = validate_rules(r)
                if w is None:
                    w = csv.DictWriter(fh, fieldnames=list(d))
                    w.writeheader()
                w.writerow(d)
        print(f"diagnostics -> {args.diagnostics}")

    deltas = _parse_deltas(args.deltas)
    n_written = 0

    raw_fh = raw_writer = None
    if args.raw:
        raw_fh = open(args.raw, "w", newline="")
        raw_writer = csv.DictWriter(raw_fh, fieldnames=OUT_FIELDS)
        raw_writer.writeheader()

    pairs_fh = pairs_writer = None
    if args.pairs:
        pairs_fh = open(args.pairs, "w", newline="")
        pairs_writer = csv.DictWriter(pairs_fh, fieldnames=PAIR_FIELDS)
        pairs_writer.writeheader()

    n_raw = n_pair_rows = 0
    try:
        with open(args.out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
            writer.writeheader()

            for key, grp in sorted(groups.items(), key=lambda kv: str(kv[0])):
                grp.sort(key=lambda r: r.rep)
                activate_model(grp[0])
                try:
                    focal = make_focal_grid(grp, deltas, args.n_focal, args.axis)
                except ValueError as exc:
                    print(f"  skip {key}: {exc}", file=sys.stderr)
                    continue

                pairs = [(a, b) for a in grp for b in grp
                         if args.allopatric or a.rep == b.rep]
                # Rows are folded into the accumulator as each pairing finishes
                # and then dropped, so memory scales with the number of cells,
                # not with pairings x focal points.
                acc = CellAccumulator()
                for protocol in args.protocols:
                    # The phenotype protocol never settles anything, so the
                    # settling convention does not apply and is not replicated.
                    modes = args.settle if protocol == "rule" else [""]
                    for mode in modes:
                        for a, b in pairs:
                            rows = assay(a, b, deltas, focal, protocol,
                                         axis=args.axis,
                                         count_roots=args.count_roots,
                                         settle=mode)
                            acc.add(rows)
                            if raw_writer is not None:
                                raw_writer.writerows(rows)
                                n_raw += len(rows)

                summary = acc.summary(n_block=args.n_block)
                writer.writerows(summary)
                n_written += len(summary)
                if pairs_writer is not None:
                    prow = acc.pair_rows()
                    pairs_writer.writerows(prow)
                    n_pair_rows += len(prow)
                print(f"  {key[0]}/{key[1]}"
                      + (f"/k={key[2]}" if key[0] == "tracking" else "")
                      + f": {len(grp)} reps, {len(pairs)} pairings, "
                        f"{len(focal)} focal x {len(deltas)} deltas "
                        f"-> {len(summary)} summary rows")
    finally:
        if raw_fh is not None:
            raw_fh.close()
        if pairs_fh is not None:
            pairs_fh.close()

    print(f"{n_written} summary rows -> {args.out}")
    if args.pairs:
        print(f"{n_pair_rows} pair rows -> {args.pairs}")
    if args.raw:
        print(f"{n_raw} raw rows -> {args.raw}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
