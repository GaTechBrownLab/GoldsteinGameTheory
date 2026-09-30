"""
test_root_finder.py — acceptance test for the exact fixed-point enumeration.

Two levels, because a dense scan over a million parameter sets is not feasible:
  1. residual test at full scale: every returned root must satisfy |g(c)| < 1e-12
     and the count must be odd-ish/sane; run over ~1e6 random parameter sets
     drawn from the regimes the simulation actually visits, plus edge cases.
  2. dense-scan comparison on a hard-case subsample: any sign change the dense
     scan brackets must correspond to a root the enumerator returned.

  python3 test_root_finder.py [n_random] [n_dense]
"""
import math, os, random, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import simulation as sim

def residual(c, bS, mS, bV, mV):
    return sim.clamp01(bS + mS * sim.clamp01(bV + mV * c)) - c

def dense_roots(bS, mS, bV, mV, n):
    lo, hi = sim.TRAIT_MIN, sim.TRAIT_MAX
    step = (hi - lo) / n
    tol = 1e-10 * max(1.0, abs(mS * mV))      # g amplifies rounding by the slope product
    out, prev_c, prev_g = [], lo, residual(lo, bS, mS, bV, mV)
    if abs(prev_g) < tol: out.append(lo)
    for i in range(1, n + 1):
        c = lo + i * step
        g = residual(c, bS, mS, bV, mV)
        if (prev_g < 0) != (g < 0):           # a strict sign change only
            a, b = prev_c, c                      # bisect the bracket
            for _ in range(60):
                m = 0.5 * (a + b)
                if (residual(a, bS, mS, bV, mV) < 0) != (residual(m, bS, mS, bV, mV) < 0):
                    b = m
                else:
                    a = m
            m = 0.5 * (a + b)
            if abs(residual(m, bS, mS, bV, mV)) < tol:   # keep it only if it really is a root
                out.append(m)
        prev_c, prev_g = c, g
    if abs(prev_g) < tol: out.append(hi)
    ded = []
    for c in sorted(out):
        if all(abs(c - d) > 1e-6 for d in ded): ded.append(c)
    return ded

def cases(n, rng):
    """Parameter sets spanning the regimes plus the awkward ones."""
    big = 1.0 / sim.ANGLE_EPS      # the steepest slope angle_to_slope can produce
    for i in range(n):
        k = i % 8
        if k == 0:                                   # everyday
            mS, mV = rng.uniform(-3, 3), rng.uniform(-3, 3)
        elif k == 1:                                 # bistable regime
            mS = rng.uniform(1, 4); mV = rng.uniform(1, 4)
            if rng.random() < 0.5: mS, mV = -mS, -mV
        elif k == 2:                                 # near-vertical rules
            mS = rng.choice([1, -1]) * rng.uniform(1e2, big)
            mV = rng.uniform(-3, 3)
        elif k == 3:                                 # both near-vertical
            mS = rng.choice([1, -1]) * rng.uniform(1e2, big)
            mV = rng.choice([1, -1]) * rng.uniform(1e2, big)
        elif k == 4:                                 # a zero slope
            mS, mV = (0.0, rng.uniform(-3, 3)) if rng.random() < 0.5 else (rng.uniform(-3, 3), 0.0)
        elif k == 5:                                 # product within 1e-6 of 1
            mS = rng.uniform(0.5, 3)
            mV = (1.0 + rng.uniform(-1e-6, 1e-6)) / mS
        elif k == 6:                                 # near-fold: roots almost touching
            mS = rng.uniform(1.01, 2.0); mV = (1.0 + rng.uniform(1e-6, 1e-3)) / mS
        else:                                        # intersection sitting on a clamp
            mS, mV = rng.uniform(-3, 3), rng.uniform(-3, 3)
        bS = rng.choice([rng.uniform(-1, 2), 0.0, sim.TRAIT_MAX, -mS * rng.uniform(0, 1)])
        bV = rng.choice([rng.uniform(-1, 2), 0.0, sim.TRAIT_MAX, -mV * rng.uniform(0, 1)])
        yield bS, mS, bV, mV

def main():
    n_random = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    n_dense = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    for model in ("minimal", "taylor"):
        sim.set_fitness_model(model)
        rng = random.Random(7)
        lo, hi = sim.TRAIT_MIN, sim.TRAIT_MAX
        bad_res = bad_dom = n = 0
        worst = 0.0
        t0 = time.time()
        for bS, mS, bV, mV in cases(n_random, rng):
            n += 1
            tol = 1e-10 * max(1.0, abs(mS * mV))
            for v, c in sim._find_all_roots(bS, mS, bV, mV):
                r = abs(residual(c, bS, mS, bV, mV))
                worst = max(worst, r / tol)
                if r > tol: bad_res += 1
                if not (lo - 1e-12 <= c <= hi + 1e-12 and lo - 1e-12 <= v <= hi + 1e-12):
                    bad_dom += 1
        print(f"[{model}] domain [{lo}, {hi}]  {n} parameter sets, {time.time()-t0:.0f}s")
        print(f"   roots failing |g| < 1e-10*max(1,|m_c m_v|): {bad_res}"
              f"   worst residual, in units of that tolerance: {worst:.2e}")
        print(f"   roots outside the model's trait domain: {bad_dom}")

        rng = random.Random(11)
        missed = extra = m = 0
        for bS, mS, bV, mV in cases(n_dense, rng):
            m += 1
            ex = sorted(c for _, c in sim._find_all_roots(bS, mS, bV, mV))
            dn = dense_roots(bS, mS, bV, mV, 20000)
            for c in dn:
                if all(abs(c - e) > 1e-5 for e in ex): missed += 1
            for c in ex:
                if all(abs(c - d) > 1e-5 for d in dn): extra += 1
        print(f"   dense-scan comparison on {m} hard cases (20,000 points each):")
        print(f"     roots the dense scan found and the enumerator missed: {missed}")
        print(f"     roots only the enumerator reports (dense scan cannot bracket): {extra}")

if __name__ == "__main__":
    main()
