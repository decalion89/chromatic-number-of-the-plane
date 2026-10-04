"""check_saved.py V.json: rebuild the field and the vectors of a saved set, check that they are unit vectors, and run
the limit test (INFEASIBLE means: for large N, G_N V has no character into [1/3, 2/3], so chi(F^2) >= 4)."""
import json, sys
from fractions import Fraction as Fr
from nf import Field
from limit_test import feasible

d = json.load(open(sys.argv[1]))
Fd = Field([[Fr(c) for c in p] for p in d["polys"]])
dec = lambda a: {tuple(int(x) for x in k.split(",")): Fr(v) for k, v in a.items()}
V = [(dec(X), dec(Y)) for X, Y in d["V"]]
for X, Y in V:
    assert Fd.add(Fd.add(Fd.mul(X, X), Fd.mul(Y, Y)), Fd.one(), -1) == {}, "not a unit vector"
ok, pat = feasible(Fd, V, want=True)
print(d["field"], len(V), "vectors:", ("feasible " + pat) if ok else "INFEASIBLE")
