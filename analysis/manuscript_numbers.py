#!/usr/bin/env python3
"""
manuscript_numbers.py -- every simulation-derived number quoted in the
manuscript, computed from the finished runs (nothing is re-simulated).

    python3 manuscript_numbers.py > output/manuscript_numbers.txt

Definitions follow Plots.R, so the numbers match the figures:
  * rows come in pre/post pairs, one pair per recorded substitution; the post
    row's dwell is the evolutionary time spent in the pre state
  * time-weighted = pre state weighted by that dwell (Fig 3F, Fig 4D-E)
  * substitution-weighted = post rows (Fig 3D, Fig 4A-B, Fig S3, S8)
  * boundary = within 0.02 of 0 or 1; near Nash = within 0.1 of 0.5
  * slope classes: 'stable' |mc mv| < 1 (Fig 3D); escalating mc mv > 1
"""
import glob, math, os, re, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

SCEN = [("EThost_ETpath", "ET/ET"), ("EThost_ERpath", "ETh/ERp"),
        ("ERhost_ETpath", "ERh/ETp"), ("ERhost_ERpath", "ER/ER")]
LAB = dict(SCEN)
BND, NEAR, SIG = 0.02, 0.1, 0.01

def hdr(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)
def pc(x): return f"{100 * x:6.2f}%"

def run_dirs(pattern):
    return [d for d in sorted(glob.glob(pattern))
            if os.path.exists(f"{d}/simulation.csv") and os.path.getsize(f"{d}/simulation.csv") > 100]

def read_runs(pattern):
    out = []
    for d in run_dirs(pattern):
        df = pd.read_csv(f"{d}/simulation.csv", low_memory=False)
        df["rep"] = int(re.search(r"_rep(\d+)", os.path.basename(d)).group(1))
        for c in ("omegaPath", "omegaHost", "mutSelCoeff", "dwell"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        out.append(df)
    return pd.concat(out, ignore_index=True) if out else None

def pairs(df):
    pre = df[df.event == "pre"].set_index(["rep", "gen"])
    post = df[df.event == "post"].set_index(["rep", "gen"])
    return pre.join(post, lsuffix="0", rsuffix="1", how="inner").reset_index()

def fH(v, c): return c * (1 - c) * (1 - v)
def fP(v, c): return v * (1 - v) * (1 - c)

def zones(x):
    b = (x < BND) | (x > 1 - BND)
    n = (~b) & (np.abs(x - 0.5) < NEAR)
    return b, n, ~(b | n)

def wmean(x, w): return float(np.sum(np.asarray(x) * w) / np.sum(w))

# ---------------------------------------------------------------------------
def scenario_block(model, root, cond, fH_=fH, fP_=fP, nashH=0.125, nashP=0.125):
    df = read_runs(f"{root}/{model}/{cond}_sigma0.01_diploid_rep*")
    if df is None: return None
    post, P = df[df.event == "post"], pairs(df)
    w = P.dwell1.values
    r = dict(n_reps=df.rep.nunique(), n_post=len(post))
    r["v_mean"], r["c_mean"] = post.v.mean(), post.s.mean()
    r["v_sd"], r["c_sd"] = post.v.std(), post.s.std()
    r["v_sd_rep"] = post.groupby("rep").v.std().mean()
    r["c_sd_rep"] = post.groupby("rep").s.std().mean()
    r["WH_mean"], r["WP_mean"] = post.hostFit.mean(), post.pathFit.mean()
    r["WH_mean_time"] = wmean(fH_(P.v0.values, P.s0.values), w)
    r["WP_mean_time"] = wmean(fP_(P.v0.values, P.s0.values), w)
    r["v_mean_time"], r["c_mean_time"] = wmean(P.v0, w), wmean(P.s0, w)
    r["WH_rel_nash"], r["WP_rel_nash"] = r["WH_mean"] / nashH, r["WP_mean"] / nashP
    r["frac_WH_below_nash"] = (post.hostFit < nashH).mean()
    r["frac_WP_below_nash"] = (post.pathFit < nashP).mean()
    r["v_min"], r["v_max"], r["c_min"], r["c_max"] = post.v.min(), post.v.max(), post.s.min(), post.s.max()
    r["v_q01"], r["v_q99"] = post.v.quantile(.01), post.v.quantile(.99)
    r["c_q01"], r["c_q99"] = post.s.quantile(.01), post.s.quantile(.99)
    r["c_absdev_max"], r["v_absdev_max"] = (post.s - .5).abs().max(), (post.v - .5).abs().max()
    for t, col in (("c", "s0"), ("v", "v0")):
        b, n, i = zones(P[col].values)
        r[f"{t}_time_boundary"], r[f"{t}_time_nearNash"], r[f"{t}_time_interior"] = wmean(b, w), wmean(n, w), wmean(i, w)
        r[f"{t}_time_low"] = wmean(P[col].values < BND, w)
        r[f"{t}_time_high"] = wmean(P[col].values > 1 - BND, w)
        x = post["s" if t == "c" else "v"].values
        r[f"{t}_sub_boundary"] = ((x < BND) | (x > 1 - BND)).mean()
    r["any_boundary_exact_sub"] = ((post.s <= 1e-9) | (post.s >= 1 - 1e-9) | (post.v <= 1e-9) | (post.v >= 1 - 1e-9)).mean()
    r["any_boundary_002_sub"] = ((post.s < BND) | (post.s > 1 - BND) | (post.v < BND) | (post.v > 1 - BND)).mean()
    # substitution rates
    for who, col, N in (("H", "omegaHost", "host"), ("P", "omegaPath", "path")):
        o = post[col]; o = o[np.isfinite(o) & (o > 0)]
        r[f"omega{who}_median"] = o.median()
        r[f"omega{who}_frac_below1"] = (o < 1).mean()
        per = post.groupby("rep")[col].apply(lambda s: s[np.isfinite(s) & (s > 0)].median())
        r[f"omega{who}_rep_medians"] = ", ".join(f"{x:.3g}" for x in per)
    r["frac_path_subs"] = (post.mutator == "path").mean()
    # per-substitution steps
    P = P.assign(dc=(P.s1 - P.s0).abs(), dv=(P.v1 - P.v0).abs())
    host_er, path_er = cond.startswith("ERhost"), cond.endswith("ERpath")
    steps = {}
    for t, col, owner, er in (("c", "dc", "host", host_er), ("v", "dv", "path", path_er)):
        own = P[P.mutator1 == owner][col]
        steps[(t, "own")] = own
        if er:
            steps[(t, "opp")] = P[P.mutator1 != owner][col]
    for (t, src), x in steps.items():
        nz = x[x > 0]
        r[f"step_{t}_{src}_median"] = nz.median()
        r[f"step_{t}_{src}_q90"] = nz.quantile(.9)
        r[f"step_{t}_{src}_q99"] = nz.quantile(.99)
        r[f"step_{t}_{src}_frac_gt_sigma"] = (x > SIG).mean()
        r[f"step_{t}_{src}_frac_gt_0.1"] = (x > .1).mean()
        r[f"step_{t}_{src}_frac_zero"] = (x == 0).mean()
    # neutral substitutions and genotype/phenotype decoupling (Fig 3A-B)
    rms = lambda db, dm: np.sqrt(np.maximum(db ** 2 + db * dm + dm ** 2 / 3, 0))
    for er, who, b, m, tr, N in ((host_er, "host", "bS", "mS", "s", 1e4), (path_er, "path", "bV", "mV", "v", 1e6)):
        if not er: continue
        E = P[P.mutator1 == who]
        geno = rms(E[b + "1"] - E[b + "0"], E[m + "1"] - E[m + "0"])
        pheno = (E[tr + "1"] - E[tr + "0"]).abs()
        ratio = geno / np.maximum(pheno, 1e-12)
        neutral = (geno > SIG) & (ratio > 10)
        nn = (E.mutSelCoeff1.abs() * N) <= 0.01
        r[f"neutral_frac_{who}"] = neutral.mean()
        keep = (~nn) & (pheno > 1e-9)
        r[f"geno_pheno_median_{who}"] = ratio[keep].median()
        r[f"near_neutral_frac_{who}"] = nn.mean()
    return r, df, post, P

def show(r, keys=None):
    for k, v in r.items():
        if keys and not any(k.startswith(x) for x in keys): continue
        print(f"  {k:32s} {v:.4g}" if isinstance(v, (float, np.floating)) else f"  {k:32s} {v}")

# ---------------------------------------------------------------------------
hdr("1. MINIMAL MODEL, main runs (data/main/minimal)")
MAIN = {}
for cond, lab in SCEN:
    res = scenario_block("minimal", "results", cond)
    if res is None: continue
    MAIN[cond] = res
    print(f"\n--- {lab} ({cond})"); show(res[0])

# ---------------------------------------------------------------------------
hdr("2. ER/ER slope products, quadrants, Nash-violation map (minimal main)")
r, df, post, P = MAIN["ERhost_ERpath"]
p = post.mS * post.mV
p0 = (P.mS0 * P.mV0).values; w = P.dwell1.values
print(f"  |mc mv| < 1 (post rows, Fig 3D 'stable')   {pc((p.abs() < 1).mean())}")
print(f"  mc mv > 1  escalating (post rows)          {pc((p > 1).mean())}")
print(f"  mc mv < -1 (post rows)                     {pc((p < -1).mean())}")
print(f"  |mc mv| < 1 time-weighted                  {pc(wmean(np.abs(p0) < 1, w))}")
print(f"  mc mv > 1 time-weighted                    {pc(wmean(p0 > 1, w))}")
print(f"  mc mv < -1 time-weighted                   {pc(wmean(p0 < -1, w))}")
q = {"policing (mc>0,mv<0)": (post.mS > 0) & (post.mV < 0), "arms race (+,+)": (post.mS > 0) & (post.mV > 0),
     "mutual retreat (-,-)": (post.mS < 0) & (post.mV < 0), "punishment (mc<0,mv>0)": (post.mS < 0) & (post.mV > 0)}
for k, m in q.items():
    print(f"  quadrant {k:28s} {pc(m.mean())}   of which mc mv>1: {pc((p[m] > 1).mean())}")
# best-response slope product at the realised state (Fig 3E)
def br_prod(v, c):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (1 - 2 * c) * (1 - 2 * v) / (4 * (1 - v) * (1 - c))
bp = br_prod(post.v.values, post.s.values); bp0 = br_prod(P.v0.values, P.s0.values)
print(f"  realised state in BR product <= -1: subs {pc(np.mean(bp <= -1))}, time {pc(wmean(bp0 <= -1, w))}")
print(f"  realised state in BR product >= +1: subs {pc(np.mean(bp >= 1))}, time {pc(wmean(bp0 >= 1, w))}")
lr = (bp <= -1) & (post.v.values > .5); ul = (bp <= -1) & (post.v.values <= .5)
print(f"     of <= -1: high-v/low-c corner {pc(lr.mean())}, low-v/high-c corner {pc(ul.mean())}")
# corners of trait space
V, C = post.v.values, post.s.values
V0, C0 = P.v0.values, P.s0.values
corner = {"c>0.5 & v>0.5": lambda c, v: (c > .5) & (v > .5), "c>0.75 & v>0.75": lambda c, v: (c > .75) & (v > .75),
          "c<0.25 & v>0.75": lambda c, v: (c < .25) & (v > .75), "c>0.75 & v<0.25": lambda c, v: (c > .75) & (v < .25),
          "c<0.5": lambda c, v: c < .5, "c<0.25": lambda c, v: c < .25, "c<0.1": lambda c, v: c < .1, "v>0.75": lambda c, v: v > .75}
for nm, f in corner.items():
    print(f"  realised {nm:18s} subs {pc(f(C, V).mean())}  time {pc(wmean(f(C0, V0), w))}")
# realised state classes by slope-product class
Pp = P.assign(p1=P.mS1 * P.mV1)
atb = lambda d: ((d.s1 <= 1e-9) | (d.s1 >= 1 - 1e-9) | (d.v1 <= 1e-9) | (d.v1 >= 1 - 1e-9))
for nm, m in (("mc mv > 1", Pp.p1 > 1), ("|mc mv| < 1", Pp.p1.abs() < 1), ("mc mv < -1", Pp.p1 < -1)):
    print(f"  realised phenotype on a trait boundary when {nm:12s}: {pc(atb(Pp[m]).mean())}  (n={m.sum()})")

# ---------------------------------------------------------------------------
hdr("3. Jumps and escalation (ER/ER)")
def jump_table(P, lab):
    d = P.assign(p1=P.mS1 * P.mV1, dc=(P.s1 - P.s0).abs(), dv=(P.v1 - P.v0).abs())
    j = (d.dc > .25) | (d.dv > .25); bj = (d.dc > .5) | (d.dv > .5)
    print(f"  [{lab}] n={len(d)}  jumps {pc(j.mean())}; of jumps at mc mv>1 {pc((d.p1[j] > 1).mean())} "
          f"(base {pc((d.p1 > 1).mean())}); big jumps {pc(bj.mean())}, at mc mv>1 {pc((d.p1[bj] > 1).mean())}; "
          f"jumps at |p|<1 {pc((d.p1[j].abs() < 1).mean())}, at p<-1 {pc((d.p1[j] < -1).mean())}")
    # entries into escalation within one substitution
    p0 = d.mS0 * d.mV0
    ent = (p0 <= 1) & (d.p1 > 1)
    sign = (np.sign(d.mS0) != np.sign(d.mS1)) | (np.sign(d.mV0) != np.sign(d.mV1))
    print(f"      entries into mc mv>1: {ent.sum()}; from mc mv<-1 {pc((p0[ent] < -1).mean())}; with a slope sign change {pc(sign[ent].mean())}")
jump_table(MAIN["ERhost_ERpath"][3], "main 3x100K, every 10th substitution")
Z = read_runs("data/zoom/minimal/ERhost_ERpath_sigma0.01_diploid_rep*_zoom")
ZP = pairs(Z)
jump_table(ZP, "zoom 8x10K, every substitution")
# episodes: maximal runs of consecutive substitutions with mc mv > 1 (zoom runs)
eps = []
for rep, g in Z[Z.event == "post"].sort_values(["rep", "gen"]).groupby("rep"):
    g = g.reset_index(drop=True)
    pp = (g.mS * g.mV).values
    esc = pp > 1
    i = 0
    while i < len(g):
        if esc[i]:
            j = i
            while j + 1 < len(g) and esc[j + 1]: j += 1
            on, last = g.iloc[i], g.iloc[j]
            exit_by = g.iloc[j + 1].mutator if j + 1 < len(g) else None
            prev_p = pp[i - 1] if i > 0 else np.nan
            eps.append(dict(rep=rep, length=j - i + 1, on_mS=on.mS, on_mV=on.mV,
                            end_c=last.s, end_v=last.v, exit_by=exit_by, onset_by=on.mutator,
                            prev_p=prev_p, max_dc=np.nan))
            i = j + 1
        else:
            i += 1
E = pd.DataFrame(eps)
zpost = Z[Z.event == "post"]
print(f"  [zoom] episodes n={len(E)}; median length {E.length.median():.0f} substitution(s), "
      f"mean {E.length.mean():.2f}, 90th pct {E.length.quantile(.9):.0f}, max {E.length.max()}")
print(f"      share of substitutions escalating {pc((zpost.mS * zpost.mV > 1).mean())}")
print(f"      onset both slopes negative {pc(((E.on_mS < 0) & (E.on_mV < 0)).mean())}, both positive {pc(((E.on_mS > 0) & (E.on_mV > 0)).mean())}")
print(f"      end state c<0.5 {pc((E.end_c < .5).mean())}; c<0.5&v>0.5 {pc(((E.end_c < .5) & (E.end_v > .5)).mean())}; "
      f"c>0.5&v>0.5 {pc(((E.end_c > .5) & (E.end_v > .5)).mean())}; c>0.5&v<0.5 {pc(((E.end_c > .5) & (E.end_v < .5)).mean())}; c<0.5&v<0.5 {pc(((E.end_c < .5) & (E.end_v < .5)).mean())}")
ex = E.exit_by.dropna()
print(f"      exits by pathogen substitution {pc((ex == 'path').mean())} (base: pathogen share of substitutions {pc((zpost.mutator == 'path').mean())})")
print(f"      onsets by host substitution {pc((E.onset_by == 'host').mean())} (base host share {pc((zpost.mutator == 'host').mean())})")
print(f"      onset from mc mv < -1 on previous substitution {pc((E.prev_p < -1).mean())}")

# ---------------------------------------------------------------------------
hdr("4. Equilibrium multiplicity under the exact enumerator (ER/ER main post rows)")
import simulation as sim
sim.set_fitness_model("minimal")
for nm, pat in (("anchored", "data/main/minimal/ERhost_ERpath_sigma0.01_diploid_rep*"),
                ("best for the mutant", "data/rule_mutant/minimal/ERhost_ERpath_sigma0.01_diploid_rep*")):
    d = read_runs(pat); d = d[d.event == "post"]
    n_eq = np.array([len(sim.find_all_equilibria(a, b, c, e)) for a, b, c, e in zip(d.bS, d.mS, d.bV, d.mV)])
    pp = (d.mS * d.mV).values
    print(f"  {nm:20s}: >=2 fixed points {pc(np.mean(n_eq >= 2))}; among mc mv>1 {pc(np.mean(n_eq[pp > 1] >= 2))}; "
          f"any with mc mv<=1 {int(np.sum((n_eq >= 2) & (pp <= 1)))}")

# ---------------------------------------------------------------------------
hdr("5. Old rule vs anchored rule, ER/ER minimal (SI robustness arm)")
for nm, pat in (("anchored", "data/main/minimal/ERhost_ERpath_sigma0.01_diploid_rep*"),
                ("best for the mutant", "data/rule_mutant/minimal/ERhost_ERpath_sigma0.01_diploid_rep*")):
    d = read_runs(pat); P = pairs(d); d = d[d.event == "post"]; w = P.dwell1.values
    pp = d.mS * d.mV
    oh = d.omegaHost[d.omegaHost > 0]; op = d.omegaPath[d.omegaPath > 0]
    bt = ((P.s0 < BND) | (P.s0 > 1 - BND))
    print(f"  {nm:20s} reps={d.rep.nunique()} escalating {pc((pp > 1).mean())}  |p|<1 {pc((pp.abs() < 1).mean())}  "
          f"any trait exactly at bound {pc(((d.s <= 1e-9) | (d.s >= 1 - 1e-9) | (d.v <= 1e-9) | (d.v >= 1 - 1e-9)).mean())}  "
          f"any trait within .02 {pc(((d.s < BND) | (d.s > 1 - BND) | (d.v < BND) | (d.v > 1 - BND)).mean())}  "
          f"c at bound (time) {pc(wmean(bt, w))}")
    print(f"  {'':20s} W_H {d.hostFit.mean():.4f}  W_P {d.pathFit.mean():.4f}  SD(c) {d.s.std():.3f}  SD(v) {d.v.std():.3f}  "
          f"median wH {oh.median():.0f}  median wP {op.median():.0f}  mean c {d.s.mean():.3f}  mean v {d.v.mean():.3f}")

# ---------------------------------------------------------------------------
hdr("6. Tempo sweep (results_gamma, N_H = N_P = 1e4)")
rows = []
for cond, lab in SCEN[1:]:
    for g, tag in ((1e-4, "_gamma0.0001"), (0.01, ""), (0.5, "_gamma0.5")):
        d = read_runs(f"data/tempo/minimal/{cond}_sigma0.01_diploid{tag}_N10000_NP10000_rep*")
        if d is None: continue
        for rep, dd in d.groupby("rep"):
            P = pairs(dd); w = P.dwell1.values
            rows.append(dict(scen=lab, gamma=g, R=(1 - g) / g, rep=rep,
                             relWH=wmean(fH(P.v0, P.s0), w) / .125, relWP=wmean(fP(P.v0, P.s0), w) / .125,
                             cbound=wmean((P.s0 < BND) | (P.s0 > 1 - BND), w)))
T = pd.DataFrame(rows)
print(T.groupby(["scen", "R"])[["relWH", "relWP", "cbound"]].agg(["mean", "min", "max"]).round(3).to_string())
def ks(a, b):
    a, b = np.sort(a), np.sort(b); x = np.concatenate([a, b])
    return float(np.max(np.abs(np.searchsorted(a, x, side="right") / len(a) - np.searchsorted(b, x, side="right") / len(b))))
A = read_runs("data/tempo/minimal/EThost_ERpath_sigma0.01_diploid_gamma0.5_N10000_NP10000_rep*"); A = A[A.event == "post"]
B = read_runs("data/tempo/minimal/ERhost_ETpath_sigma0.01_diploid_gamma0.5_N10000_NP10000_rep*"); B = B[B.event == "post"]
print(f"  R=1 KS D: rule trait {ks(A.v, B.s):.3f}, fixed trait {ks(A.s, B.v):.3f}, rule fitness {ks(A.pathFit, B.hostFit):.3f}, fixed fitness {ks(A.hostFit, B.pathFit):.3f}")
for nm, rule_fit, fixed_fit, rule_tr in (("ETh/ERp", A.pathFit, A.hostFit, A.v), ("ERh/ETp", B.hostFit, B.pathFit, B.s)):
    print(f"  R=1 {nm}: rule player below Nash fitness {pc((rule_fit < .125).mean())}; fixed partner below {pc((fixed_fit < .125).mean())}, "
          f"fixed min {fixed_fit.min():.4f} ({fixed_fit.min() / .125:.3f} of Nash); rule trait max {rule_tr.max():.4f}; "
          f"rule mean fitness/Nash {rule_fit.mean() / .125:.3f}; fixed mean/Nash {fixed_fit.mean() / .125:.3f}")

# ---------------------------------------------------------------------------
hdr("7. Tracking sweep (data/main/tracking)")
rows = []
for cond, lab in SCEN:
    for d_ in run_dirs(f"data/main/tracking/{cond}_sigma0.01_diploid_k*_rep*"):
        k = float(re.search(r"_k([0-9.]+)_rep", d_).group(1))
        dd = pd.read_csv(f"{d_}/simulation.csv", low_memory=False); dd = dd[dd.event == "post"]
        rows.append(dict(scen=lab, k=k, v_sd=dd.v.std(), r=np.corrcoef(dd.hostFit, dd.pathFit)[0, 1] if dd.hostFit.std() > 0 else np.nan,
                         w_h=dd.hostFit.mean(), v_mean=dd.v.mean(), c_mean=dd.s.mean(),
                         v_hi=(dd.v > 1 - BND).mean(), c_hi=(dd.s > 1 - BND).mean(),
                         v_bnd=((dd.v < BND) | (dd.v > 1 - BND)).mean(), c_bnd=((dd.s < BND) | (dd.s > 1 - BND)).mean(), n=len(dd)))
K = pd.DataFrame(rows)
ref = K[K.scen == "ET/ET"].groupby("k").w_h.mean()
K["rel_wh"] = K.w_h / K.k.map(ref)
pd.set_option("display.width", 200)
print(K.groupby(["scen", "k"])[["v_sd", "r", "w_h", "rel_wh", "v_mean", "c_mean", "v_hi", "c_hi", "v_bnd", "c_bnd", "n"]].mean().round(3).to_string())

# ---------------------------------------------------------------------------
hdr("8. Acute model (data/main/acute)")
aH = lambda v, c: sim._host_acute(v, c)
sim.set_fitness_model("acute")
vH = np.vectorize(sim._host_acute); vP = np.vectorize(sim._path_acute)
def golden(f, lo=0.0, hi=1.0, tol=1e-12):
    g = (math.sqrt(5) - 1) / 2; a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    while b - a > tol:
        if f(c) > f(d): b = d
        else: a = c
        c, d = b - g * (b - a), a + g * (b - a)
    return (a + b) / 2
def nash_and_product(model):
    sim.set_fitness_model(model)
    H, Pf = sim.host_fitness, sim.path_fitness
    v, c = 0.5, 0.5
    for _ in range(500):
        c2 = golden(lambda x: H(v, x)); v2 = golden(lambda x: Pf(x, c2))
        if abs(c2 - c) + abs(v2 - v) < 1e-12: c, v = c2, v2; break
        c, v = c2, v2
    h = 1e-4
    Hcc = (H(v, c + h) - 2 * H(v, c) + H(v, c - h)) / h ** 2
    Hcv = (H(v + h, c + h) - H(v + h, c - h) - H(v - h, c + h) + H(v - h, c - h)) / (4 * h * h)
    Pvv = (Pf(v + h, c) - 2 * Pf(v, c) + Pf(v - h, c)) / h ** 2
    Pvc = (Pf(v + h, c + h) - Pf(v + h, c - h) - Pf(v - h, c + h) + Pf(v - h, c - h)) / (4 * h * h)
    s1, s2 = -Hcv / Hcc, -Pvc / Pvv
    return v, c, H(v, c), Pf(v, c), s1, s2, s1 * s2
for m in ("minimal", "acute", "chronic"):
    v, c, wh, wp, s1, s2, pr = nash_and_product(m)
    print(f"  {m:8s} Nash v*={v:.4f} c*={c:.4f}  W_H*={wh:.4f} W_P*={wp:.4f}  dc*/dv={s1:.4f} dv*/dc={s2:.4f}  product={pr:.4f}")
sim.set_fitness_model("acute")
v_, c_, whN, wpN, *_ = nash_and_product("acute")
sim.set_fitness_model("acute")
AC = {}
for cond, lab in SCEN:
    res = scenario_block("acute", "results", cond, fH_=vH, fP_=vP, nashH=whN, nashP=wpN)
    if res is None: continue
    AC[cond] = res
    print(f"\n--- acute {lab}")
    show(res[0], keys=["n_reps", "v_mean", "c_mean", "v_sd", "c_sd", "WH", "WP", "frac_W", "c_time", "v_time", "omega", "step_c_own_median", "step_v_own_median", "step_c_opp_q99", "step_v_opp_q99"])
d = AC["ERhost_ERpath"][2]; pp = d.mS * d.mV
print(f"\n  acute ER/ER escalating {pc((pp > 1).mean())}, |p|<1 {pc((pp.abs() < 1).mean())}, "
      f"any trait within .02 of bound {pc(((d.s < BND) | (d.s > 1 - BND) | (d.v < BND) | (d.v > 1 - BND)).mean())}")
sim.set_fitness_model("minimal")

# ---------------------------------------------------------------------------
hdr("9. Time-shift assay (data/timeshift)")
tz = pd.read_csv("data/timeshift/timeshift_zoom.csv"); tm = pd.read_csv("data/timeshift/timeshift_main.csv")
tz = tz[tz.fitness_model == "minimal"]; tm = tm[tm.fitness_model == "minimal"]
def rule(d): return d[(d.protocol == "rule") & (d.settle == "flow-contemporary")]
for cond, lab in SCEN:
    rz = rule(tz[tz.condition == cond]); sym = rz[rz.sympatric == 1].sort_values("delta"); allo = rz[rz.sympatric == 0]
    peak = sym.W_P_mean[sym.delta == 0].iloc[0]; base = allo.W_P_mean.mean()
    ex = sym.set_index("delta").W_P_mean - base; half = (peak - base) / 2
    left = [dl for dl in sorted(ex.index[ex.index < 0], reverse=True) if ex[dl] < half]
    right = [dl for dl in sorted(ex.index[ex.index > 0]) if ex[dl] < half]
    ph = tz[(tz.condition == cond) & (tz.protocol == "phenotype")]
    phs = ph[ph.sympatric == 1]; pha = ph[ph.sympatric == 0]
    print(f"  {lab:8s} rule: sym peak(0) {peak:.4f}  sym(-200) {sym.W_P_mean.iloc[0]:.4f}  sym(+200) {sym.W_P_mean.iloc[-1]:.4f}  allo mean {base:.4f}  "
          f"half-max at {left[0] if left else 'none'} / +{right[0] if right else 'none'}")
    print(f"  {'':8s} phenotype: sym(0) {phs.W_P_mean[phs.delta == 0].iloc[0]:.4f}  sym mean {phs.W_P_mean.mean():.4f} (range {phs.W_P_mean.min():.4f}-{phs.W_P_mean.max():.4f})  allo mean {pha.W_P_mean.mean() if len(pha) else float('nan'):.4f}")
rm = rule(tm[tm.delta == 0])
for cond, lab in SCEN:
    s1 = rm[(rm.condition == cond) & (rm.sympatric == 1)].iloc[0]; s0 = rm[(rm.condition == cond) & (rm.sympatric == 0)].iloc[0]
    print(f"  main partition {lab:8s} pathogen position {(s1.E_f - s0.E_f) * s0.E_g:+.4f}  host level {s1.E_f * (s1.E_g - s0.E_g):+.4f}  "
          f"matching {s1.cov_fg - s0.cov_fg:+.4f}  total {s1.W_P_mean - s0.W_P_mean:+.4f}")
TT = {2: 4.303, 7: 2.365}
for which in ("main", "zoom"):
    pr = pd.read_csv(f"data/timeshift/timeshift_{which}_pairs.csv")
    pr = pr[(pr.fitness_model == "minimal") & (pr.protocol == "phenotype") & (pr.delta == 0)]
    for cond, lab in SCEN:
        q = pr[pr.condition == cond]
        g = q.groupby("rep_path").apply(lambda x: x.W_P_mean[x.sympatric == 1].mean() - x.W_P_mean[x.sympatric == 0].mean(), include_groups=False)
        n = len(g); est, se = g.mean(), g.std() / math.sqrt(n)
        t = TT.get(n - 1, 1.96)
        print(f"  matching ({which}, phenotype, delta 0) {lab:8s} {est:+.5f} +/- {se:.5f} SE (n={n}); 95% CI [{est - t * se:+.5f}, {est + t * se:+.5f}]")
for cond, lab in SCEN:
    post = MAIN[cond][2]
    f, g = post.v * (1 - post.v), 1 - post.s
    print(f"  sd(f) sd(g) bound {lab:8s} {f.std() * g.std():.5f}   sd(f) {f.std():.4f} sd(g) {g.std():.4f}")

# ---------------------------------------------------------------------------
hdr("10. Extras quoted in the text")
for cond, lab in SCEN:
    d = AC[cond][2]; P = AC[cond][3]; w = P.dwell1.values
    print(f"  acute {lab:8s} any trait exactly at a bound (subs) {pc(((d.s <= 1e-9) | (d.s >= 1 - 1e-9) | (d.v <= 1e-9) | (d.v >= 1 - 1e-9)).mean())}; "
          f"c>0.5&v>0.5 subs {pc(((d.s > .5) & (d.v > .5)).mean())} time {pc(wmean((P.s0.values > .5) & (P.v0.values > .5), w))}; "
          f"c>0.75&v>0.75 subs {pc(((d.s > .75) & (d.v > .75)).mean())} time {pc(wmean((P.s0.values > .75) & (P.v0.values > .75), w))}")
ph = tz[(tz.condition == "ERhost_ERpath") & (tz.protocol == "phenotype") & (tz.sympatric == 1) & (tz.delta != 0)]
print(f"  ER/ER phenotype-protocol sympatric level away from delta 0: {ph.W_P_mean.min():.4f}-{ph.W_P_mean.max():.4f} (mean {ph.W_P_mean.mean():.4f})")
for cond, lab in SCEN[1:]:
    post = MAIN[cond][2]
    print(f"  {lab:8s} c<0.25 {pc((post.s < .25).mean())}  v<0.25 {pc((post.v < .25).mean())}  v<0.1 {pc((post.v < .1).mean())}")
