#!/usr/bin/env python3
"""
missing_runs.py — print the run_experiments.py commands for every run in the
manuscript design that is not already complete on disk.

Resume-safe: a run counts as done when its simulation.csv has the expected row
count (2 rows per recorded substitution, plus a header). Use it after a code fix
that only invalidates part of the design — archive the affected runs, then let
this fill the gaps.

    python3 scripts/missing_runs.py | tr '\\n' '\\0' | xargs -0 -n 1 -P 10 bash -c
    python3 scripts/missing_runs.py --count
"""
import argparse, os, sys

CONDS = ["ERhost_ERpath", "ERhost_ETpath", "EThost_ERpath", "EThost_ETpath"]
MIXED = ["ERhost_ETpath", "EThost_ERpath"]
REPS = 4
ZOOM_REPS = 8
KS = ["0", "0.5", "1", "1.5", "2", "3", "4"]
GAMMAS = ["0.0001", "0.01", "0.5", "0.99"]      # 0.99 added: tests the slower player in both directions


def rows_expected(gens, write_every):
    return 2 * (gens // write_every) + 1


def done(path, gens, write_every):
    try:
        with open(path) as fh:
            return sum(1 for _ in fh) == rows_expected(gens, write_every)
    except OSError:
        return False


def dirname(cond, sigma="0.01", diploid=True, gamma=None, N=None, NP=None,
            k=None, rep=None, tag=None, extra=()):
    parts = [cond, f"sigma{sigma}"]
    if diploid: parts.append("diploid")
    if gamma is not None and abs(float(gamma) - 0.01) > 1e-9: parts.append(f"gamma{gamma}")
    if N is not None: parts.append(f"N{N}")
    if NP is not None: parts.append(f"NP{NP}")
    if k is not None: parts.append(f"k{float(k)}")
    parts += list(extra)
    if rep is not None: parts.append(f"rep{rep}")
    if tag is not None: parts.append(tag)
    return "_".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", action="store_true", help="print a summary instead of commands")
    ap.add_argument("--only", nargs="*", default=None,
                    help="restrict to these sets: main zoom tempo tracking random oldrule")
    args = ap.parse_args()
    cmds, summary = [], {}

    def add(set_name, out, model, cond, gens, write, rep, log, **kw):
        # only these reach the directory name; burn/raw are command-line only
        name_kw = {k: kw[k] for k in ("gamma", "N", "NP", "k", "tag", "extra") if k in kw}
        d = os.path.join(out, model if model != "tracking" else "tracking",
                         dirname(cond, rep=rep, **name_kw))
        if done(os.path.join(d, "simulation.csv"), gens, write):
            return
        flags = [f"-f {model}", f"-c {cond}", "--diploid", "--step-sizes 0.01",
                 f"--rep {rep}", f"--write-every {write}", f"--gens {gens}"]
        if kw.get("k") is not None: flags.insert(1, f"--k {kw['k']}")
        if kw.get("gamma") is not None: flags.append(f"--gamma {kw['gamma']}")
        if kw.get("N") is not None: flags.append(f"--pop-size {kw['N']}")
        if kw.get("NP") is not None: flags.append(f"--path-pop-size {kw['NP']}")
        if kw.get("tag") is not None: flags.append(f"--tag {kw['tag']}")
        if kw.get("burn") is not None: flags.append(f"--burn-in {kw['burn']}")
        for f in kw.get("raw", []): flags.append(f)
        cmds.append(f"python3 src/run_experiments.py {' '.join(flags)} -o {out} "
                    f"> logs/rebuild/{log}.log 2>&1")
        summary[set_name] = summary.get(set_name, 0) + 1

    want = lambda s: args.only is None or s in args.only

    if want("main"):
        for model in ("minimal", "acute"):
            for cond in CONDS:
                for rep in range(1, REPS + 1):
                    add("main", "data/main", model, cond, 100000, 10, rep,
                        f"main_{model}_{cond}_rep{rep}", burn=10000)
    if want("zoom"):
        for cond in CONDS:
            for rep in range(1, ZOOM_REPS + 1):
                add("zoom", "data/zoom", "minimal", cond, 10000, 1, rep,
                    f"zoom_{cond}_rep{rep}", tag="zoom")
    if want("tempo"):
        for g in GAMMAS:
            for cond in ("ERhost_ERpath", "ERhost_ETpath", "EThost_ERpath"):
                for rep in range(1, REPS + 1):
                    add("tempo", "data/tempo", "minimal", cond, 10000, 10, rep,
                        f"tempo_g{g}_{cond}_rep{rep}", gamma=g, N=10000, NP=10000, burn=10000)
    if want("tracking"):
        for k in KS:
            for cond in CONDS:
                for rep in range(1, REPS + 1):
                    add("tracking", "data/main", "tracking", cond, 10000, 10, rep,
                        f"track_k{k}_{cond}_rep{rep}", k=k, burn=10000)
    if want("random"):
        for cond in CONDS:
            for rep in range(1, REPS + 1):
                add("random", "data/random", "minimal", cond, 100000, 10, rep,
                    f"random_{cond}_rep{rep}", burn=10000, extra=("random",),
                    raw=["--proposals random"])
    if want("oldrule"):
        for rep in range(1, REPS + 1):
            # note: eq_selection == "mutant" adds no directory tag
            add("oldrule", "data/rule_mutant", "minimal", "ERhost_ERpath", 100000, 10, rep,
                f"oldrule_rep{rep}", burn=10000, raw=["--eq-selection mutant"])

    if args.count:
        for k, v in sorted(summary.items()):
            print(f"{k:10s} {v:4d} runs missing")
        print(f"{'total':10s} {len(cmds):4d}")
    else:
        print("\n".join(cmds))


if __name__ == "__main__":
    main()
