"""Periodic gates on a rich unit module: can a finite quotient refute the
directions the twisted colourings could not?

Same unit set as richmod.py.  For each n: sample functionals psi : M -> Z/n
nonzero on every unit vector; for each with a 5-colourable Cayley graph, try to
keep a 2e pair alike (apart gate) / split a 5e pair (pair gate), per direction.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, itertools, random, time
from fractions import Fraction as Fr
from math import gcd
import numpy as np
from pysat.solvers import Solver
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "richmod.py")).read().split("PS = np.array")[0])
Cm = np.array(Cm, dtype=object)
MODS = [int(x) for x in sys.argv[3].split(",")]
MODE = sys.argv[4] if len(sys.argv) > 4 else "apart"
def colour(n, S, t, same):
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s in S:
            w = (v + s) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    if t is not None:
        for c in range(5):
            cl.append([-X(0, c), X(t % n, c)] if same else [-X(0, c), -X(t % n, c)])
    return Solver(name="cd19", bootstrap_with=cl).solve()
random.seed(5)
alive = set(range(len(E)))
for n in MODS:
    adm = col = 0; seen = {}
    for trial in range(60000):
        psi = [random.randrange(n) for _ in range(r)]
        img = [int(sum(int(a) * int(b) for a, b in zip(psi, c)) % n) for c in Cm]
        if 0 in img: continue
        adm += 1
        S = frozenset(img) | frozenset((-x) % n for x in img)
        if S not in seen: seen[S] = colour(n, S, None, True)
        if not seen[S]: continue
        col += 1
        for k in list(alive):
            if MODE == "apart":
                t = (2 * img[k]) % n
                if t == 0 or colour(n, S, t, True): alive.discard(k)
            else:
                t = (5 * img[k]) % n
                if t != 0 and colour(n, S, t, False): alive.discard(k)
        if not alive: break
    print(f"  Z/{n}: {adm} admissible, {col} with 5-colourable Cayley graph; directions still unrefuted: {len(alive)} of {len(E)}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
