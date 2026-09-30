"""
diag_equilibrium.py -- which fixed point does the ER solver realise, and would
behavioural dynamics have got there?

Runs an ordinary ER/ER simulation (the trajectory is NOT altered: the solver
and selection rule in simulation.py are used as-is) and, at every accepted
substitution, shadows the realised phenotype with three alternative readings
of the reaction norms, all started from the RESIDENT phenotype (the state the
interaction was in just before the mutant arose):

  ct[r]  continuous-time behavioural dynamics
             dc/dt = r   * (clamp(bS + mS*v) - c)
             dv/dt = 1.0 * (clamp(bV + mV*c) - v)
         for host:pathogen rate ratios r in RATE_RATIOS. Stability of an
         interior point is rate-independent (mS*mV < 1), but in the bistable
         regime (mS*mV > 1) the basin boundary -- and so which corner is
         reached -- depends on r, so we report several.
  disc   discrete lagged, simultaneous full-step update
             (c, v) <- (clamp(bS + mS*v), clamp(bV + mV*c))
         which is the reading behind the |mS*mV| < 1 condition.

Usage (defaults = production minimal-model settings in simulation.py):
    python diag_equilibrium.py --gens 20000 --burn 10000 --model minimal
    python diag_equilibrium.py --csv results/.../simulation.csv   # sign split only

Outputs a per-substitution CSV and a printed summary.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import simulation as sim

RATE_RATIOS = (0.1, 1.0, 10.0)
AGREE_TOL = 1e-3


# --------------------------------------------------------------------------
# Behavioural dynamics
# --------------------------------------------------------------------------

def _rhs(c, v, bS, mS, bV, mV, r):
    cl = sim.clamp01
    return r * (cl(bS + mS * v) - c), (cl(bV + mV * c) - v)


def flow_ct(c, v, bS, mS, bV, mV, r, tol=1e-11, max_steps=400_000):
    """RK4 integration of the clamped linear response dynamics from (c, v).
    Returns (c, v, converged). Negative divergence everywhere => no cycles,
    so the flow always ends on a fixed point."""
    h = 0.05 / max(r, 1.0)
    for _ in range(max_steps):
        k1 = _rhs(c, v, bS, mS, bV, mV, r)
        k2 = _rhs(c + .5*h*k1[0], v + .5*h*k1[1], bS, mS, bV, mV, r)
        k3 = _rhs(c + .5*h*k2[0], v + .5*h*k2[1], bS, mS, bV, mV, r)
        k4 = _rhs(c + h*k3[0], v + h*k3[1], bS, mS, bV, mV, r)
        dc = h * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6
        dv = h * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6
        c, v = sim.clamp01(c + dc), sim.clamp01(v + dv)
        if abs(dc) < tol and abs(dv) < tol:
            return c, v, True
    return c, v, False


def iterate_discrete(c, v, bS, mS, bV, mV, tol=1e-12, max_iter=5000):
    """Simultaneous lagged full-step update. Returns (c, v, status, amp)."""
    cl = sim.clamp01
    hist = [(c, v)]
    for _ in range(max_iter):
        c, v = cl(bS + mS * v), cl(bV + mV * c)
        for lag in (1, 2, 4):
            if len(hist) >= lag:
                c0, v0 = hist[-lag]
                if abs(c - c0) < tol and abs(v - v0) < tol:
                    if lag == 1:
                        return c, v, "fixed", 0.0
                    cyc = hist[-lag:]
                    amp = max(p[0] for p in cyc) - min(p[0] for p in cyc)
                    return (sum(p[0] for p in cyc) / lag,
                            sum(p[1] for p in cyc) / lag, f"cycle{lag}", amp)
        hist.append((c, v))
        if len(hist) > 8:
            hist.pop(0)
    return c, v, "noconv", float("nan")


def pclass(p):
    if p < -1:
        return "p<-1"
    if p > 1:
        return "p>+1"
    return "|p|<1"


# --------------------------------------------------------------------------
# Shadowing the simulation
# --------------------------------------------------------------------------

def run(args):
    sim.set_fitness_model(args.model)
    if args.model == "tracking":
        sim.TRACKING_K = args.k
    sim.FIX_HOST_REACTIVITY = False
    sim.FIX_PATH_REACTIVITY = False
    sim.USE_SIMPLE_PROPOSALS = False
    # Match the production runs; the module defaults are sigma 0.1 and haploid
    sim.std_dev_move = args.sigma
    sim.DIPLOID_KIMURA = not args.haploid
    if args.gamma is not None:
        sim.prob_host_mutate = args.gamma
    print(f"  model={args.model} sigma={sim.std_dev_move} "
          f"diploid={sim.DIPLOID_KIMURA} gamma={sim.prob_host_mutate} "
          f"seed={args.seed}", flush=True)

    records = []
    state = {"gen": None}
    orig = sim.Simulation._refresh_equilibrium

    def patched(self, selector, track=False):
        if not track or state["gen"] is None or state["gen"] < 0:
            return orig(self, selector, track)
        c_res, v_res = self.s, self.v           # resident phenotype
        bS, mS, bV, mV = self.bS, self.mS, self.bV, self.mV
        eqs = sim.find_all_equilibria(bS, mS, bV, mV)
        out = orig(self, selector, track)       # simulation proceeds unchanged
        c_sol, v_sol = self.s, self.v
        p = mS * mV
        unc = sim._solve_interior(bS, mS, bV, mV)
        rec = {
            "gen": state["gen"], "mutator": selector,
            "mS": mS, "mV": mV, "p": p, "pclass": pclass(p),
            "n_roots": len(eqs),
            "unclamped_in_box": (unc is not None and 0 < unc[0] < 1 and 0 < unc[1] < 1),
            "c_res": c_res, "v_res": v_res, "c_sol": c_sol, "v_sol": v_sol,
            "sol_interior": 0 < c_sol < 1 and 0 < v_sol < 1,
            "sol_is_saddle": (p > 1 and 0 < c_sol < 1 and 0 < v_sol < 1
                              and len(eqs) > 1),
            "W_H_sol": sim.host_fitness(v_sol, c_sol),
            "W_P_sol": sim.path_fitness(v_sol, c_sol),
        }
        for r in RATE_RATIOS:
            if len(eqs) == 1:     # unique root is globally attracting in ct
                cc, vv, ok = c_sol, v_sol, True
            else:
                cc, vv, ok = flow_ct(c_res, v_res, bS, mS, bV, mV, r)
            tag = f"ct{r:g}"
            rec[f"c_{tag}"], rec[f"v_{tag}"], rec[f"ok_{tag}"] = cc, vv, ok
            rec[f"agree_{tag}"] = (abs(cc - c_sol) < AGREE_TOL
                                   and abs(vv - v_sol) < AGREE_TOL)
            rec[f"W_H_{tag}"] = sim.host_fitness(vv, cc)
            rec[f"W_P_{tag}"] = sim.path_fitness(vv, cc)
        cd, vd, st, amp = iterate_discrete(c_res, v_res, bS, mS, bV, mV)
        rec.update({"c_disc": cd, "v_disc": vd, "disc_status": st,
                    "disc_amp": amp,
                    "agree_disc": abs(cd - c_sol) < AGREE_TOL
                                  and abs(vd - v_sol) < AGREE_TOL})
        records.append(rec)
        return out

    sim.Simulation._refresh_equilibrium = patched

    rng = random.Random(args.seed)
    s = sim.Simulation(evolved_strategy=True, rng=rng)
    t0 = time.time()
    for gen in range(-args.burn, args.gens):
        state["gen"] = gen
        s.step_generation()
        if gen % max(1, (args.gens + args.burn) // 20) == 0:
            print(f"  gen {gen:>8}  v={s.v:.3f} c={s.s:.3f}  "
                  f"mS*mV={s.mS*s.mV:+.2f}  [{time.time()-t0:.0f}s]", flush=True)
    sim.Simulation._refresh_equilibrium = orig
    return records


# --------------------------------------------------------------------------
# Summaries
# --------------------------------------------------------------------------

def pct(a, b):
    return f"{100*a/b:5.1f}%" if b else "   --"


def summarise(recs):
    n = len(recs)
    print(f"\n=== {n} accepted substitutions after burn-in ===\n")
    by = Counter(r["pclass"] for r in recs)
    print("Slope-product class of realised states (after each substitution):")
    for k in ("p<-1", "|p|<1", "p>+1"):
        print(f"  {k:6s} {by[k]:>8}  {pct(by[k], n)}")
    print(f"  => '|mc mv| > 1' as used in text/plots: "
          f"{pct(by['p<-1'] + by['p>+1'], n)};  signed 'mc mv > 1': {pct(by['p>+1'], n)}")

    print("\nRoot structure and what the solver realised, by class:")
    hdr = f"  {'class':6s} {'n':>7} {'multi-root':>11} {'sol=saddle':>11} " \
          f"{'sol@bound':>10} {'unclamped out of box':>21}"
    print(hdr)
    for k in ("p<-1", "|p|<1", "p>+1"):
        g = [r for r in recs if r["pclass"] == k]
        m = len(g)
        print(f"  {k:6s} {m:>7} {pct(sum(r['n_roots'] > 1 for r in g), m):>11} "
              f"{pct(sum(r['sol_is_saddle'] for r in g), m):>11} "
              f"{pct(sum(not r['sol_interior'] for r in g), m):>10} "
              f"{pct(sum(not r['unclamped_in_box'] for r in g), m):>21}")

    print("\nAgreement of solver phenotype with behavioural dynamics from the resident:")
    cols = [f"ct{r:g}" for r in RATE_RATIOS] + ["disc"]
    print(f"  {'class':6s} " + " ".join(f"{c:>9}" for c in cols))
    for k in ("p<-1", "|p|<1", "p>+1", "all"):
        g = recs if k == "all" else [r for r in recs if r["pclass"] == k]
        print(f"  {k:6s} " + " ".join(
            f"{pct(sum(r['agree_' + c] for r in g), len(g)):>9}" for c in cols))

    st = Counter(r["disc_status"] for r in recs)
    print("\nDiscrete-update outcome:", dict(st))
    neg = [r for r in recs if r["pclass"] == "p<-1"]
    if neg:
        cyc = sum(r["disc_status"].startswith("cycle") for r in neg)
        print(f"  p<-1 states that cycle under discrete update: {pct(cyc, len(neg))}"
              f" (all converge in continuous time)")

    print("\nMean fitness: solver vs continuous-time (r=1) phenotype")
    for k in ("p>+1", "all"):
        g = recs if k == "all" else [r for r in recs if r["pclass"] == k]
        if not g:
            continue
        m = len(g)
        print(f"  {k:5s} W_H {sum(r['W_H_sol'] for r in g)/m:.4f} -> "
              f"{sum(r['W_H_ct1'] for r in g)/m:.4f}   "
              f"W_P {sum(r['W_P_sol'] for r in g)/m:.4f} -> "
              f"{sum(r['W_P_ct1'] for r in g)/m:.4f}")

    # The flow only approaches a clamped boundary asymptotically, so it stops a
    # hair inside it. A strict 0 < x < 1 test counts those as interior and badly
    # understates boundary occupancy under the behavioural dynamics.
    def on_bound(x, tol=1e-6):
        return x <= tol or x >= 1 - tol
    bnd_sol = sum(on_bound(r["c_sol"]) or on_bound(r["v_sol"]) for r in recs)
    bnd_ct = sum(on_bound(r["c_ct1"]) or on_bound(r["v_ct1"]) for r in recs)
    print(f"\nBoundary occupancy (either trait clamped): solver {pct(bnd_sol, n)}, "
          f"continuous-time r=1 {pct(bnd_ct, n)}")
    nc = sum(not r[f"ok_ct{x:g}"] for r in recs for x in RATE_RATIOS)
    if nc:
        print(f"WARNING: {nc} continuous-time integrations did not converge")


def sign_split_csv(path):
    """Sign split of an existing simulation.csv (post rows)."""
    cnt, n = Counter(), 0
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if row.get("event") != "post":
                continue
            try:
                p = float(row["mS"]) * float(row["mV"])
            except (ValueError, KeyError):
                continue
            cnt[pclass(p)] += 1
            n += 1
    print(f"{path}: {n} post rows")
    for k in ("p<-1", "|p|<1", "p>+1"):
        print(f"  {k:6s} {cnt[k]:>8}  {pct(cnt[k], n)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="minimal")
    ap.add_argument("--k", type=float, default=0.0)
    ap.add_argument("--gens", type=int, default=20_000)
    ap.add_argument("--burn", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=3248232)
    ap.add_argument("--sigma", type=float, default=0.01,
                    help="Mutation step size (production: 0.01)")
    ap.add_argument("--haploid", action="store_true",
                    help="Haploid Kimura; default is diploid, as in production")
    ap.add_argument("--gamma", type=float, default=None)
    ap.add_argument("--out", default="diag_equilibrium.csv")
    ap.add_argument("--csv", nargs="+", default=None,
                    help="Only sign-split existing simulation.csv file(s).")
    ap.add_argument("--resummarise", nargs="+", default=None,
                    help="Re-print the summary from diagnostic CSV(s) already written.")
    args = ap.parse_args()

    if args.csv:
        for p in args.csv:
            sign_split_csv(p)
        return

    if args.resummarise:
        recs = []
        for p in args.resummarise:
            with open(p) as fh:
                for row in csv.DictReader(fh):
                    rec = {}
                    for k, v in row.items():
                        if v in ("True", "False"):
                            rec[k] = v == "True"
                        elif k in ("mutator", "pclass", "disc_status"):
                            rec[k] = v
                        else:
                            rec[k] = float(v)
                    recs.append(rec)
        print(f"({len(args.resummarise)} file(s) pooled)")
        summarise(recs)
        return

    recs = run(args)
    if recs:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]))
            w.writeheader()
            w.writerows(recs)
        print(f"\nPer-substitution records -> {args.out}")
    summarise(recs)


if __name__ == "__main__":
    main()
