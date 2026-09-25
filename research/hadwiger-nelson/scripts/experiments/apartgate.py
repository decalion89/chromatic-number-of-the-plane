"""The apart gate: can a periodic 5-colouring colour a 2e pair ALIKE?

If some finite quotient psi : M -> Z/n, nonzero on every unit vector, has a
proper 5-colouring of its Cayley graph with c(0) = c(psi(2e)), then no
unit-distance graph on M can force u, u + 2e apart -- the Exoo-Ismailescu
gadget cannot live on M.  Coset colourings mod 5 always split 2e, so only
larger quotients can refute it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, random, itertools, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
ROOT = HN_DIR
MODE = sys.argv[1]; MODS = [int(x) for x in sys.argv[2:]]
if MODE == "ei":
    sols = [(a,b,c,d) for a in range(-7,8) for b in range(-4,5) for c in range(-12,13) for d in range(-3,4)
            if 3*a*a + 11*b*b + c*c + 33*d*d == 144 and a*b == -c*d]
    E = list({max(v, tuple(-x for x in v)) for v in sols})
else:
    d = json.load(open(f"{ROOT}/data/{MODE}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    E = edge_vectors(build_graph(P))
B = echelon(E); C = [coords(B, v) for v in E]; r = len(B)
print(f"{MODE}: {len(E)} directions, rank {r}", flush=True)
def colour(n, S, t):
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s in S:
            w = (v + s) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    if t is not None:
        for c in range(5): cl.append([-X(0, c), X(t % n, c)])   # c(0) == c(t)
    return Solver(name="cd19", bootstrap_with=cl).solve()
def colour2(n, S, t):
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s_ in S:
            w = (v + s_) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    for c in range(5): cl.append([-X(0, c), -X(t % n, c)])   # c(0) != c(t)
    return Solver(name="cd19", bootstrap_with=cl).solve()
random.seed(11); t0 = time.time()
MULT = int(__import__("os").environ.get("MULT", "2")); SAME = __import__("os").environ.get("SAME", "1") == "1"
alive = set(int(x) for x in __import__("os").environ["DIRS"].split(",")) if __import__("os").environ.get("DIRS") else set(range(len(C)))
for n in MODS:
    adm = col = 0; hit = None; seen = {}
    space = n ** r
    it = itertools.product(range(n), repeat=r) if space <= 300000 else (tuple(random.randrange(n) for _ in range(r)) for _ in range(300000))
    for psi in it:
        img = [sum(a * b for a, b in zip(psi, c)) % n for c in C]
        if 0 in img: continue
        adm += 1
        S = frozenset(img) | frozenset((-x) % n for x in img)
        if S not in seen: seen[S] = colour(n, S, None)
        if not seen[S]: continue
        col += 1
        for k in list(alive):
            t = (MULT * img[k]) % n
            if SAME:
                bad = (t == 0) or colour(n, S, t)
            else:
                bad = (t != 0) and colour2(n, S, t)
            if bad: alive.discard(k)
        if not alive: break
    print(f"  Z/{n}: {adm} admissible, {col} 5-colourable; directions still unrefuted: {sorted(alive)}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
