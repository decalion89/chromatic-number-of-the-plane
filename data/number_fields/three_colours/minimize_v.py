"""Greedy deletion: smallest subset of a saved infeasible set V (keeping v = 1) that stays infeasible in the limit test."""
import json, sys
from fractions import Fraction as Fr
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from nf import Field
from limit_test import feasible
d = json.load(open(sys.argv[1]))
Fd = Field([[Fr(c) for c in p] for p in d["polys"]])
dec = lambda a: {tuple(int(x) for x in k.split(",")): Fr(v) for k, v in a.items()}
V = [(dec(X), dec(Y)) for X, Y in d["V"]]
assert not feasible(Fd, V)
keep = list(range(len(V)))
# try deleting the vectors with the largest denominators first
def size(k):
    X, Y = V[k]
    return max([v.denominator for v in list(X.values()) + list(Y.values())] + [1])
order = sorted(range(1, len(V)), key=lambda k: -size(k))
for k in order:
    trial = [j for j in keep if j != k]
    if not feasible(Fd, [V[j] for j in trial]):
        keep = trial
        print("dropped", k, "->", len(keep), flush=True)
print("minimal:", keep, [size(k) for k in keep])
enc = lambda a: {",".join(map(str, k)): str(v) for k, v in a.items()}
json.dump({"field": d["field"], "polys": d["polys"], "V": [[enc(V[j][0]), enc(V[j][1])] for j in keep]},
          open(sys.argv[2], "w"))
