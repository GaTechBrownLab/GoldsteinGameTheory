#!/usr/bin/env python3
"""
Goldstein Host-Pathogen Coevolution - Clean Experiment Runner
==============================================================

Runs all 4 experimental conditions with clear output naming:
  - EThost_ETpath: Both fixed traits (EI mode)
  - EThost_ERpath: Host fixed, pathogen reactive  
  - ERhost_ETpath: Host reactive, pathogen fixed
  - ERhost_ERpath: Both reactive (ES mode)

Output structure:
  results/{FITNESS_MODEL}/{CONDITION}/
    simulation.csv      # Full trajectory
    config.json         # Run parameters
    
Usage:
  python run_experiments.py                           # All conditions, acute model
  python run_experiments.py --fitness minimal         # All conditions, minimal model  
  python run_experiments.py --condition ERhost_ERpath # Single condition
  python run_experiments.py --quick                   # Quick test (10K gens)
  python run_experiments.py --list                    # Show available options
"""

import os
import sys
import json
import random
import datetime
import argparse
import importlib

# =============================================================================
# CONFIGURATION
# =============================================================================

# Experimental conditions: (host_reactive, path_reactive)
CONDITIONS = {
    "EThost_ETpath": (False, False),  # EI mode - both fixed traits
    "EThost_ERpath": (False, True),   # Host fixed (mS=0), path reactive
    "ERhost_ETpath": (True,  False),  # Host reactive, path fixed (mV=0)  
    "ERhost_ERpath": (True,  True),   # ES mode - both reactive
}

FITNESS_MODELS = ["acute", "chronic", "minimal", "taylor"]

DEFAULT_PARAMS = {
    "max_gens": 1_000_000,
    "burn_in_gens": 10_000,
    "seed_base": 3248232,
    "write_every": 100,
}

QUICK_PARAMS = {
    "max_gens": 10_000,
    "burn_in_gens": 1_000,
    "seed_base": 3248232,
    "write_every": 10,
}

# =============================================================================
# SIMULATION RUNNER
# =============================================================================

def run_single_condition(
    fitness_model: str,
    host_reactive: bool,
    path_reactive: bool,
    params: dict,
    output_base: str = "results",
    sigma: float = None,
    diploid: bool = False,
    fix_host_trait: float = None,
    fix_path_trait: float = None,
    gamma: float = None,
) -> str:
    """Run a single experimental condition.
    
    Parameters
    ----------
    gamma : float, optional
        Mutation rate asymmetry (prob_host_mutate).
        0.01 = pathogen-fast (default in Goldstein 2020)
        0.50 = equal rates
        0.99 = host-fast (reversed asymmetry)
        None = use simulation.py default (0.01).
    """
    

    # Import the unified simulation module
    sim = importlib.import_module("simulation")
    sim.set_fitness_model(fitness_model)

    
    # Determine condition name
    host_str = "ERhost" if host_reactive else "EThost"
    path_str = "ERpath" if path_reactive else "ETpath"
    condition = f"{host_str}_{path_str}"
    
    # Build output directory with tags for non-default settings
    dir_name = condition
    tags = []
    if sigma is not None:
        tags.append(f"sigma{sigma}")
    if diploid:
        tags.append("diploid")
    if gamma is not None and abs(gamma - 0.01) > 1e-6:
        tags.append(f"gamma{gamma}")
    if fix_host_trait is not None:
        tags.append(f"fixH{fix_host_trait}")
    if fix_path_trait is not None:
        tags.append(f"fixP{fix_path_trait}")
    if tags:
        dir_name += "_" + "_".join(tags)
    
    # Setup output
    output_dir = os.path.join(output_base, fitness_model, dir_name)
    os.makedirs(output_dir, exist_ok=True)
    output_csv = os.path.join(output_dir, "simulation.csv")
    config_json = os.path.join(output_dir, "config.json")
    
    # Determine mode: EI (both ET) or ES (any reactive)
    if not host_reactive and not path_reactive:
        evolved_strategy = False
        fix_host_react = False
        fix_path_react = False
    else:
        evolved_strategy = True
        fix_host_react = not host_reactive  # Fix reactivity if NOT reactive
        fix_path_react = not path_reactive
    
    # === SET SIMULATION GLOBALS ===
    sim.USE_GAUSSIAN = False
    sim.FIX_HOST_REACTIVITY = fix_host_react
    sim.FIX_PATH_REACTIVITY = fix_path_react
    sim.DIPLOID_KIMURA = diploid
    sim.FIX_HOST_TRAIT = fix_host_trait
    sim.FIX_PATH_TRAIT = fix_path_trait
    if sigma is not None:
        sim.std_dev_move = sigma
    if gamma is not None:
        sim.prob_host_mutate = gamma
    sim.max_gens = params["max_gens"]
    sim.burn_in_gens = params["burn_in_gens"]
    sim.seed = params["seed_base"]
    sim.WRITE_EVERY_FAST = params["write_every"]
    sim.RUNTIME_MODE = "fast"
    
    # Build config record
    config = {
        "condition": condition,
        "fitness_model": fitness_model,
        "host_reactive": host_reactive,
        "path_reactive": path_reactive,
        "evolved_strategy": evolved_strategy,
        "FIX_HOST_REACTIVITY": fix_host_react,
        "FIX_PATH_REACTIVITY": fix_path_react,
        "USE_BOUNDED_TRAITS": sim.USE_BOUNDED_TRAITS,
        "TRAIT_MIN": sim.TRAIT_MIN,
        "TRAIT_MAX": sim.TRAIT_MAX,
        "std_dev_move": sim.std_dev_move,
        "DIPLOID_KIMURA": diploid,
        "HOST_POP_N": sim.HOST_POP_N,
        "PATH_POP_N": sim.PATH_POP_N,
        "prob_host_mutate": sim.prob_host_mutate,
        "FIX_HOST_TRAIT": fix_host_trait if fix_host_trait is not None else False,
        "FIX_PATH_TRAIT": fix_path_trait if fix_path_trait is not None else False,
        "timestamp": datetime.datetime.now().isoformat(),
        "parameters": params.copy(),
    }
    
    with open(config_json, "w") as f:
        json.dump(config, f, indent=2)
    
    # Print info
    print(f"\n{'='*60}")
    print(f"Running: {fitness_model} / {dir_name}")
    print(f"{'='*60}")
    print(f"  Host: {'ER (reactive, mS evolves)' if host_reactive else 'ET (fixed trait)'}")
    print(f"  Path: {'ER (reactive, mV evolves)' if path_reactive else 'ET (fixed trait)'}")
    print(f"  Internal mode: {'ES' if evolved_strategy else 'EI'}")
    if evolved_strategy:
        print(f"  FIX_HOST_REACTIVITY={fix_host_react}, FIX_PATH_REACTIVITY={fix_path_react}")
    print(f"  Fitness model: {fitness_model}")
    print(f"  Trait domain: [{sim.TRAIT_MIN}, {sim.TRAIT_MAX}]")
    print(f"  Step size (σ): {sim.std_dev_move}")
    print(f"  Mutation asymmetry (γ): {sim.prob_host_mutate}"
          f"  ({'pathogen-fast' if sim.prob_host_mutate < 0.5 else 'host-fast' if sim.prob_host_mutate > 0.5 else 'equal'})")
    if diploid:
        print(f"  Kimura: DIPLOID (4Ns)")
    if fix_host_trait is not None:
        print(f"  Host trait PINNED at s={fix_host_trait}")
    if fix_path_trait is not None:
        print(f"  Path trait PINNED at v={fix_path_trait}")
    print(f"  Output: {output_csv}")
    print(f"  Generations: {params['max_gens']:,}")
    
    # === RUN SIMULATION ===
    rng = random.Random(params["seed_base"])
    simulation = sim.Simulation(evolved_strategy=evolved_strategy, rng=rng)
    sim.run_with_runtime_modes(simulation, output_csv)
    
    # Print stats
    if hasattr(simulation, 'interior_count_total'):
        total = simulation.interior_count_total + simulation.boundary_count_total
        if total > 0:
            pct = 100 * simulation.interior_count_total / total
            print(f"  Interior equilibria: {pct:.1f}%")
    
    print(f"  Config saved: {config_json}")
    
    return output_csv


def run_all_conditions(
    fitness_model: str,
    params: dict,
    output_base: str = "results",
    sigma: float = None,
    diploid: bool = False,
    fix_host_trait: float = None,
    fix_path_trait: float = None,
    gamma: float = None,
) -> dict:
    """Run all 4 conditions for a fitness model."""
    results = {}
    
    for condition, (host_reactive, path_reactive) in CONDITIONS.items():
        output_csv = run_single_condition(
            fitness_model=fitness_model,
            host_reactive=host_reactive,
            path_reactive=path_reactive,
            params=params,
            output_base=output_base,
            sigma=sigma,
            diploid=diploid,
            fix_host_trait=fix_host_trait,
            fix_path_trait=fix_path_trait,
            gamma=gamma,
        )
        results[condition] = output_csv
    
    return results


def run_gamma_sweep(
    fitness_model: str,
    params: dict,
    gammas: list = None,
    conditions: list = None,
    output_base: str = "results",
    sigma: float = 0.1,
    diploid: bool = True,
) -> dict:
    """Sweep mutation rate asymmetry (gamma) across conditions.
    
    This probes whether fast/slow dynamics shape the distribution of 
    realized phenotypes, or whether volatility is intrinsic to ER status.
    
    Parameters
    ----------
    gammas : list of float
        Values of prob_host_mutate to sweep.
        Default: [0.01, 0.1, 0.5, 0.9, 0.99]
    conditions : list of str, optional
        Which conditions to run. Default: all 4.
        For a focused experiment, use ["ERhost_ERpath"] or
        ["EThost_ERpath", "ERhost_ETpath", "ERhost_ERpath"].
    """
    if gammas is None:
        gammas = [0.01, 0.1, 0.5, 0.9, 0.99]
    if conditions is None:
        conditions = list(CONDITIONS.keys())
    
    results = {}
    total = len(gammas) * len(conditions)
    done = 0
    
    print(f"\n{'#'*60}")
    print(f"# Gamma Sweep: {fitness_model}")
    print(f"# γ values: {gammas}")
    print(f"# Conditions: {conditions}")
    print(f"# Total runs: {total}")
    print(f"{'#'*60}")
    
    for gamma in gammas:
        for cond_name in conditions:
            if cond_name not in CONDITIONS:
                print(f"  WARNING: Unknown condition '{cond_name}', skipping")
                continue
            host_reactive, path_reactive = CONDITIONS[cond_name]
            done += 1
            print(f"\n  [{done}/{total}] γ={gamma}, {cond_name}")
            
            output_csv = run_single_condition(
                fitness_model=fitness_model,
                host_reactive=host_reactive,
                path_reactive=path_reactive,
                params=params,
                output_base=output_base,
                sigma=sigma,
                diploid=diploid,
                gamma=gamma,
            )
            results[(gamma, cond_name)] = output_csv
    
    print(f"\n{'='*60}")
    print(f"GAMMA SWEEP COMPLETE — {done} runs")
    print(f"{'='*60}")
    for (g, c), path in results.items():
        print(f"  γ={g:5.3f}  {c:20s}  →  {path}")
    
    return results


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run Goldstein coevolution experiments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_experiments.py                           # All 4 conditions, acute
  python run_experiments.py --fitness minimal         # All 4, minimal model
  python run_experiments.py --condition ERhost_ERpath # Single condition  
  python run_experiments.py --quick                   # Quick test (10K gens)
  python run_experiments.py --list                    # Show options

Output:
  results/{fitness}/{condition}[_tags]/simulation.csv
  results/{fitness}/{condition}[_tags]/config.json

Step-size sweep:
  python run_experiments.py -f acute -c EThost_ETpath --step-sizes 0.01,0.05,0.2,0.5,1
Diploid:
  python run_experiments.py -f acute --diploid 
Pin host trait:
  python run_experiments.py -f acute -c EThost_ETpath --fix-host 0.5
Single gamma (host-fast):
  python run_experiments.py -f acute --gamma 0.99 --diploid
Gamma sweep (all conditions):
  python run_experiments.py -f acute --gamma-sweep 0.01,0.5,0.99 --diploid
Gamma sweep (ER-ER only):
  python run_experiments.py -f acute -c ERhost_ERpath --gamma-sweep 0.01,0.1,0.5,0.9,0.99 --diploid
        """
    )
    
    parser.add_argument("--fitness", "-f", choices=FITNESS_MODELS, default="acute",
                        help="Fitness model (default: acute)")
    parser.add_argument("--condition", "-c", choices=list(CONDITIONS.keys()),
                        help="Single condition (default: all 4)")
    parser.add_argument("--quick", "-q", action="store_true",
                        help="Quick test (10K gens)")
    parser.add_argument("--output", "-o", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--seed", type=int, help="Random seed override")
    parser.add_argument("--list", "-l", action="store_true",
                        help="List available options")
    parser.add_argument("--step-sizes", type=str, default=None,
                        help="Comma-separated step sizes to sweep, e.g. 0.01,0.05,0.2,0.5")
    parser.add_argument("--diploid", action="store_true",
                        help="Use diploid semi-dominant Kimura (4Ns denominator)")
    parser.add_argument("--fix-host", type=float, default=None,
                        help="Pin host clearance at this value (e.g. 0.5)")
    parser.add_argument("--fix-path", type=float, default=None,
                        help="Pin pathogen virulence at this value (e.g. 0.3)")
    parser.add_argument("--gamma", type=float, default=None,
                        help="Mutation rate asymmetry (prob_host_mutate). "
                             "0.01=pathogen-fast (default), 0.5=equal, 0.99=host-fast")
    parser.add_argument("--gamma-sweep", type=str, default=None,
                        help="Comma-separated gamma values to sweep, "
                             "e.g. 0.01,0.1,0.5,0.9,0.99")
   
    args = parser.parse_args()
    
    if args.list:
        print("\nFitness models:")
        for fm in FITNESS_MODELS:
            print(f"  {fm}")
        print("\nConditions:")
        for cond, (hr, pr) in CONDITIONS.items():
            h = "ER" if hr else "ET"
            p = "ER" if pr else "ET"
            print(f"  {cond}: Host={h}, Path={p}")
        print("\nOutput: results/{fitness}/{condition}/simulation.csv")
        return
    
    params = QUICK_PARAMS.copy() if args.quick else DEFAULT_PARAMS.copy()
    if args.seed:
        params["seed_base"] = args.seed
    
    # Parse step sizes
    sigmas = [None]  # default: use simulation.py default (0.1)
    if args.step_sizes:
        sigmas = [float(x.strip()) for x in args.step_sizes.split(",")]
    
    # Parse gamma sweep
    gammas = None
    if args.gamma_sweep:
        gammas = [float(x.strip()) for x in args.gamma_sweep.split(",")]
    
    print(f"\n{'#'*60}")
    print(f"# Goldstein Coevolution Experiment Runner")
    print(f"{'#'*60}")
    print(f"Fitness: {args.fitness}")
    print(f"Mode: {'QUICK (10K gens)' if args.quick else 'FULL (1M gens)'}")
    if args.diploid:
        print(f"Kimura: DIPLOID (4Ns)")
    if len(sigmas) > 1 or sigmas[0] is not None:
        print(f"Step sizes: {sigmas}")
    if args.gamma is not None:
        print(f"Mutation asymmetry (γ): {args.gamma}")
    if gammas is not None:
        print(f"Gamma sweep: {gammas}")
    if args.fix_host is not None:
        print(f"Host pinned at s={args.fix_host}")
    if args.fix_path is not None:
        print(f"Path pinned at v={args.fix_path}")
    
    # === Gamma sweep mode ===
    if gammas is not None:
        conditions_to_run = [args.condition] if args.condition else None
        for sigma in sigmas:
            run_gamma_sweep(
                fitness_model=args.fitness,
                params=params,
                gammas=gammas,
                conditions=conditions_to_run,
                output_base=args.output,
                sigma=sigma if sigma is not None else 0.1,
                diploid=args.diploid,
            )
    else:
        # === Standard mode (with optional single gamma) ===
        for sigma in sigmas:
            if args.condition:
                host_reactive, path_reactive = CONDITIONS[args.condition]
                run_single_condition(
                    fitness_model=args.fitness,
                    host_reactive=host_reactive,
                    path_reactive=path_reactive,
                    params=params,
                    output_base=args.output,
                    sigma=sigma,
                    diploid=args.diploid,
                    fix_host_trait=args.fix_host,
                    fix_path_trait=args.fix_path,
                    gamma=args.gamma,
                )
            else:
                results = run_all_conditions(
                    fitness_model=args.fitness,
                    params=params,
                    output_base=args.output,
                    sigma=sigma,
                    diploid=args.diploid,
                    fix_host_trait=args.fix_host,
                    fix_path_trait=args.fix_path,
                    gamma=args.gamma,
                )
                print(f"\n{'='*60}")
                print("COMPLETE - Output files:")
                print(f"{'='*60}")
                for cond, path in results.items():
                    print(f"  {cond}: {path}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
