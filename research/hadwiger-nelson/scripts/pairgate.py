"""The pair gate along one direction: can a periodic 5-colouring split 5e?

For each modulus n, sample functionals psi : M -> Z/n nonzero on every unit
direction; for each with a 5-colourable Cayley graph, ask for a 5-colouring
with c(0) != c(psi(5e)).  A hit refutes forcing (A, A + 5e) on this module.
"""
import sys, json, random, time
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
NAME = sys.argv[1]; GDIV = int(sys.argv[2]); MODS = [int(x) for x in sys.argv[3:]]
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
E = edge_vectors(g); B = echelon(E); C = [coords(B, v) for v in E]; r = len(B)
ks = [k for k, c in enumerate(C) if reduce(gcd, [abs(t) for t in c]) == GDIV]
print(f"{NAME}: rank {r}, {len(C)} directions; target directions with e in {GDIV}M: {ks}", flush=True)
def colour(n, S, t):
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s in S:
            w = (v + s) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    if t is not None:
        for c in range(5): cl.append([-X(0, c), -X(t % n, c)])
    return Solver(name="cd19", bootstrap_with=cl).solve()
random.seed(7); t0 = time.time()
for n in MODS:
    adm = col = 0; hit = None; seen = {}
    for trial in range(200000):
        psi = [random.randrange(n) for _ in range(r)]
        img = [sum(a * b for a, b in zip(psi, c)) % n for c in C]
        if 0 in img: continue
        adm += 1
        S = frozenset(img) | frozenset((-x) % n for x in img)
        if S not in seen: seen[S] = colour(n, S, None)
        if not seen[S]: continue
        col += 1
        for k in ks:
            t = (5 * img[k]) % n
            if t and colour(n, S, t): hit = (psi, k, t); break
        if hit: break
    print(f"  Z/{n}: {adm} admissible, {col} with 5-colourable Cayley graph: "
          f"{'SPLIT ' + str(hit) if hit else 'no split along the target'}   [{time.time()-t0:.0f}s]", flush=True)
