"""minimize.py d D: greedy deletion of units from U_D keeping the relation-space Lemma W test INFEASIBLE."""
import sys, json, time, random
from fractions import Fraction as Fr
from kappaD import units, halve
from kapparel import relation_basis, solve_rel
d, D = int(sys.argv[1]), int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
V = halve(units(d, D))
def infeasible(W):
    ker = relation_basis(W)
    name, _ = solve_rel(len(W), ker, Fr(1, 3), 120)
    return name == "INFEASIBLE", name
ok, nm = infeasible(V); print("start", len(V), nm, flush=True); assert ok
order = list(range(len(V))); random.Random(seed).shuffle(order)
keep = list(V)
for idx in order:
    u = V[idx]
    trial = [w for w in keep if w != u]
    ok, nm = infeasible(trial)
    if ok:
        keep = trial
    print(f"  drop {u}: {nm} -> keep {len(keep)}", flush=True)
print("MINIMAL", len(keep), keep)
json.dump({"d": d, "D": D, "units": keep}, open(f"min_{d}_{D}_s{seed}.json", "w"))
