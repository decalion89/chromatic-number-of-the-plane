"""core_min.py CORE.json OUT.json [seed]: deletion-based minimal infeasible subset of the period system (Lemma F17
ranges on p in Hom(Rel, Z)), longest relations tried first."""
import json, sys, random
from fractions import Fraction as Fr
from math import floor, ceil
from ortools.sat.python import cp_model
core = json.load(open(sys.argv[1])); rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
n = len(core[0]["a"])
def feasible(rels):
    md = cp_model.CpModel()
    p = [md.NewIntVar(-200, 200, f"p{j}") for j in range(n)]
    for r in rels:
        P = sum(x for x in r["rho"] if x > 0); N = -sum(x for x in r["rho"] if x < 0)
        e = sum(a * p[j] for j, a in enumerate(r["a"]) if a)
        md.Add(e >= floor(Fr(P - 3 * N, 4)) + 1); md.Add(e <= ceil(Fr(3 * P - N, 4)) - 1)
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = 120
    st = s.Solve(md)
    assert st in (cp_model.OPTIMAL, cp_model.FEASIBLE, cp_model.INFEASIBLE), s.StatusName(st)
    return st != cp_model.INFEASIBLE
assert not feasible(core)
order = sorted(range(len(core)), key=lambda i: (-core[i]["l1"], rng.random()))
keep = set(range(len(core)))
for i in order:
    trial = [core[j] for j in keep if j != i]
    if not feasible(trial):
        keep.discard(i)
res = [core[j] for j in sorted(keep)]
json.dump(res, open(sys.argv[2], "w"))
from collections import Counter
print(len(res), "relations;", sorted(Counter(r["l1"] for r in res).items()))
