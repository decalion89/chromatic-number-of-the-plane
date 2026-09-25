"""The full-quotient gate: 5-colour the Cayley graph of M/qM on the unit vectors.

A proper 5-colouring of Cay(M/qM, U) pulls back to a qM-periodic 5-colouring
of the whole unit-distance graph on M.  Its consequences are automatic:
  q = 2       : c(x + 2e) = c(x) for every unit e  -> no 2e gadget anywhere in M,
                and c(x + 5e) = c(x + e) != c(x)    -> no forced 5e pair in M;
  q = 3, 4, 6 : 5 = +-1 mod q, so c(x + 5e) = c(x +- e) != c(x) -> no forced 5e pair.
For other q (and to see what is left) the colouring is tested pair by pair.
The cyclic Z/n gate only sees quotients M -> Z/n; this one sees all of M/qM.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, time, os
from fractions import Fraction as Fr
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
ROOT = HN_DIR
name = sys.argv[1]; QS = [int(x) for x in sys.argv[2:]]
d = json.load(open(f"{ROOT}/data/{name}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
E = edge_vectors(build_graph(P))
B = echelon(E); C = [coords(B, v) for v in E]; r = len(B)
print(f"{name}: {len(E)} directions, rank {r}", flush=True)
t0 = time.time()
for q in QS:
    N = q ** r
    enc = lambda v: sum((x % q) * q ** i for i, x in enumerate(v))
    img = [enc(c) for c in C]
    if 0 in img:
        print(f"  M/{q}M: a unit lies in {q}M -> loop, no {q}M-periodic colouring", flush=True); continue
    S = set()
    for c in C:
        S.add(tuple(x % q for x in c)); S.add(tuple((-x) % q for x in c))
    S = sorted(S)
    print(f"  M/{q}M: {N} vertices, {len(S)} generators, {N*len(S)//2} edges", flush=True)
    if N * len(S) > 4e7:
        print("    too large, skipped", flush=True); continue
    digits = list(itertools.product(range(q), repeat=r))
    index = {v: i for i, v in enumerate(digits)}
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(N)]
    for i, v in enumerate(digits):
        for s in S:
            w = index[tuple((a + b) % q for a, b in zip(v, s))]
            if i < w:
                for c in range(5): cl.append([-X(i, c), -X(w, c)])
    cl.append([X(0, 0)])                                   # colour symmetry: vertex 0 gets colour 0
    s0 = index[S[0]]; cl.append([X(s0, 1)])               # and its first neighbour colour 1
    sol = Solver(name="cd19", bootstrap_with=cl)
    ok = sol.solve()
    print(f"    5-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
    if not ok: continue
    m = sol.get_model(); col = {}
    for i in range(N):
        for c in range(5):
            if m[X(i, c) - 1] > 0: col[i] = c; break
    # which pairs does this colouring keep alike / split, direction by direction?
    alike2 = sum(1 for c in C if col[index[tuple((2 * x) % q for x in c)]] == col[0])
    split5 = sum(1 for c in C if col[index[tuple((5 * x) % q for x in c)]] != col[0])
    print(f"    one colouring: keeps 2e alike at the origin in {alike2}/{len(C)} directions, "
          f"splits 5e at the origin in {split5}/{len(C)}", flush=True)
