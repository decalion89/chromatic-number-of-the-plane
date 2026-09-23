"""Can a periodic 5-colouring split a 5M pair?  If so, no integral growth along
these directions can ever force it.

Coset colourings mod 5 never split u, u + 5e.  But a colouring through a larger
quotient psi : M -> Z/n (n = 10, 15, 20, 25, ...) with psi(5e) != 0, whose
Cayley graph on psi(directions) has a 5-colouring with c(0) != c(psi(5e)),
colours every finite graph on these directions and splits the pair.
"""
import sys, json, itertools, random, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
from math import gcd
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
NAME = sys.argv[1]; MODS = [int(x) for x in sys.argv[2:]]
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
E = edge_vectors(g); B = echelon(E); C = [coords(B, v) for v in E]
r = len(B)
# the 5M target: 5 times each edge direction
print(f"{NAME}: {len(E)} directions, rank {r}", flush=True)
def cayley_split(n, S, t):
    """5-colouring of Cay(Z/n, S) with c(0) != c(t)?"""
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s in S:
            w = (v + s) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    for c in range(5): cl.append([-X(0, c), -X(t % n, c)])
    return Solver(name="cd19", bootstrap_with=cl).solve()
random.seed(0)
t0 = time.time()
for n in MODS:
    found = None; tried = 0; admissible = 0; seen = {}
    space = n ** r
    it = itertools.product(range(n), repeat=r) if space <= 400000 else (tuple(random.randrange(n) for _ in range(r)) for _ in range(400000))
    for psi in it:
        tried += 1
        img = [sum(a * b for a, b in zip(psi, c)) % n for c in C]
        if 0 in img: continue
        admissible += 1
        S = frozenset(img) | frozenset((-x) % n for x in img)
        for k, c in enumerate(C):
            t = (5 * img[k]) % n
            if t == 0: continue
            key = (S, t)
            if key not in seen: seen[key] = cayley_split(n, S, t)
            if seen[key]:
                found = (psi, k, t); break
        if found: break
    print(f"  Z/{n}: {tried} functionals tried ({'all' if space <= 400000 else 'random'}), {admissible} admissible: "
          f"{'a periodic 5-colouring SPLITS a 5M pair: ' + str(found) if found else 'no periodic 5-colouring splits any 5e pair'}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
