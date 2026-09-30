"""
test_dimensional_weight.py — every player carries the same dimensional weight.

A player's proposals are a quadrature of one mutation distribution, so the total
neutral weight must be num_step_bins^2 for both players in every scenario:
  ET/ET            51 trait steps x factor 51
  pinned player    51 trait steps x factor 51   (mixed scenarios)
  plastic player   51 x 51 steps  x factor 1

Before the fix the pinned player in a mixed scenario carried 51 instead of 2601,
i.e. a 51-fold shortfall in mutational supply that scaled with the grid size.

    python3 test_dimensional_weight.py
"""
import os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import simulation as sim

CONDS = {"EThost_ETpath": (False, False), "EThost_ERpath": (False, True),
         "ERhost_ETpath": (True, False),  "ERhost_ERpath": (True, True)}

def main():
    sim.set_fitness_model("minimal")
    sim.std_dev_move, sim.DIPLOID_KIMURA = 0.01, True
    expected = sim.num_step_bins ** 2
    ok = True
    print(f"expected total neutral weight per player: {expected}\n")
    print(f"{'scenario':16s} {'host':>10} {'pathogen':>10}   result")
    for cond, (hr, pr) in CONDS.items():
        sim.FIX_HOST_REACTIVITY = not hr
        sim.FIX_PATH_REACTIVITY = not pr
        s = sim.Simulation(evolved_strategy=(hr or pr), rng=random.Random(1))
        for _ in range(20):
            s.step_generation()
        _, _, nh = s._evaluate_host_mutations()
        _, _, np_ = s._evaluate_path_mutations()
        good = abs(nh - expected) < 1e-9 and abs(np_ - expected) < 1e-9
        ok &= good
        print(f"{cond:16s} {nh:10.0f} {np_:10.0f}   {'ok' if good else 'MISMATCH'}")
    sim.FIX_HOST_REACTIVITY = sim.FIX_PATH_REACTIVITY = False
    print("\nPASS" if ok else "\nFAIL")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
