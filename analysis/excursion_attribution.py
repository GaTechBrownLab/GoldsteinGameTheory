"""
excursion_attribution.py — what carries the punctuated excursions?

Two candidate drivers for a large excursion in the realised phenotype:
  escalation   the slope product crosses +1, so the rules admit several
               phenotypes and the pair can run to an extreme
  steepness    one slope is large, so a small shift in the partner's rule is
               amplified: dc = m_c * delta / (1 - m_c m_v), regardless of sign

Tabulates excursions by slope-product class and by amplification, for runs made
under either equilibrium-selection rule.

  python3 excursion_attribution.py results_anchor/minimal/*/simulation.csv
"""
import csv, os, sys
from collections import defaultdict

W_P_LOW = 0.10      # pathogen fitness collapses
D_TRAIT = 0.25      # one substitution moves a trait this far
TOL     = 1e-6

# The two kinds of event are not the same and must not be pooled: a jump is a
# large single-substitution move of the realised phenotype, a dip is low
# pathogen fitness, which is common in ordinary states.
DEFS = {
    "jump: a trait moves > %.2f in one substitution" % D_TRAIT:
        lambda r: r["dc"] > D_TRAIT or r["dv"] > D_TRAIT,
    "big jump: > 0.50":
        lambda r: r["dc"] > 0.5 or r["dv"] > 0.5,
    "dip: W_P < %.2f" % W_P_LOW:
        lambda r: r["W_P"] < W_P_LOW,
}

def pclass(p):
    return "p<-1" if p < -1 else ("p>+1" if p > 1 else "|p|<1")

def abin(a):
    return "<1" if a < 1 else ("1-3" if a < 3 else ("3-10" if a < 10 else ">10"))

def scan(path):
    rows = []
    pre = None
    for row in csv.DictReader(open(path)):
        if row["event"] == "pre":
            pre = row; continue
        if pre is None: continue
        f = lambda r, k: float(r[k])
        mS, mV = f(row, "mS"), f(row, "mV")
        bS, bV = f(row, "bS"), f(row, "bV")
        c, v = f(row, "s"), f(row, "v")
        dc, dv = c - f(pre, "s"), v - f(pre, "v")
        pre = None
        p = mS * mV
        den = abs(1.0 - p)
        # amplification of a partner shift into each trait
        amp = max(abs(mS), abs(mV)) / den if den > 1e-12 else float("inf")
        saddle = False
        if p > 1 and den > 1e-12:
            cu, vu = (bS + mS * bV) / (1 - p), (bV + mV * bS) / (1 - p)
            if 0 < cu < 1 and 0 < vu < 1:
                saddle = abs(c - cu) < TOL and abs(v - vu) < TOL
        rows.append(dict(p=p, pclass=pclass(p), steep=max(abs(mS), abs(mV)),
                         amp=amp, W_P=f(row, "pathFit"), W_H=f(row, "hostFit"),
                         dc=abs(dc), dv=abs(dv), saddle=saddle,
                         boundary=(c <= 1e-9 or c >= 1-1e-9 or v <= 1e-9 or v >= 1-1e-9)))
    return rows

def pct(a, b): return f"{100*a/b:5.1f}%" if b else "   --"

def report(label, rows):
    n = len(rows)
    print(f"\n=== {label}: {n} substitutions ===")
    print(f"\n  {'event':50s} {'rate':>7} {'p>+1':>7} {'|p|<1':>7} {'p<-1':>7} {'amp>3':>7}")
    for name, f in DEFS.items():
        e = [r for r in rows if f(r)]
        if not e: continue
        sh = lambda k: 100 * sum(r["pclass"] == k for r in e) / len(e)
        print(f"  {name:50s} {100*len(e)/n:6.1f}% {sh('p>+1'):6.1f}% {sh('|p|<1'):6.1f}% "
              f"{sh('p<-1'):6.1f}% {100*sum(r['amp']>3 for r in e)/len(e):6.1f}%")
    sh = lambda k: 100 * sum(r["pclass"] == k for r in rows) / n
    print(f"  {'(all states, for comparison)':50s} {100.0:6.1f}% {sh('p>+1'):6.1f}% {sh('|p|<1'):6.1f}% "
          f"{sh('p<-1'):6.1f}% {100*sum(r['amp']>3 for r in rows)/n:6.1f}%")

    # enrichment tables use the jump definition: those are the punctuated events
    exc = [r for r in rows if r["dc"] > D_TRAIT or r["dv"] > D_TRAIT]
    print(f"\n  breakdown of jumps ({len(exc)} events)")
    print(f"\n  {'slope product':14s} {'% of states':>12} {'% of excursions':>16} {'enrichment':>11}")
    for k in ("p<-1", "|p|<1", "p>+1"):
        a = sum(r["pclass"] == k for r in rows)
        b = sum(r["pclass"] == k for r in exc)
        enr = (b / len(exc)) / (a / n) if a and exc else float("nan")
        print(f"  {k:14s} {pct(a, n):>12} {pct(b, len(exc)):>16} {enr:10.2f}x")
    print(f"\n  {'amplification':14s} {'% of states':>12} {'% of excursions':>16} {'enrichment':>11}")
    for k in ("<1", "1-3", "3-10", ">10"):
        a = sum(abin(r["amp"]) == k for r in rows)
        b = sum(abin(r["amp"]) == k for r in exc)
        enr = (b / len(exc)) / (a / n) if a and exc else float("nan")
        print(f"  {k:14s} {pct(a, n):>12} {pct(b, len(exc)):>16} {enr:10.2f}x")
    # the decisive cell for the steepness hypothesis
    key = [r for r in exc if r["pclass"] != "p>+1" and r["steep"] >= 3]
    print(f"\n  excursions with NO escalation (m_c m_v < 1) but a slope >= 3: "
          f"{pct(len(key), len(exc))} of excursions")
    print(f"  excursions that are at a trait boundary: {pct(sum(r['boundary'] for r in exc), len(exc))}")
    sad = sum(r["saddle"] for r in rows)
    print(f"  realised states that are the interior repelling root: {pct(sad, n)}  ({sad} states)")

def main(paths):
    groups = defaultdict(list)
    for p in paths:
        run = os.path.basename(os.path.dirname(p))   # not the parent tree's name
        rule = "anchored" if "_anchor" in run else "best for the mutant"
        groups[rule] += scan(p)
    for rule in ("best for the mutant", "anchored"):
        if groups[rule]:
            report(rule, groups[rule])

if __name__ == "__main__":
    main(sys.argv[1:])
