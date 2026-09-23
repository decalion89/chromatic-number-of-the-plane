"""Periodic splitting on E-I's FULL coordinate module: all 30 unit vectors."""
import sys, itertools, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
from pysat.solvers import Solver
sols = [(a,b,c,d) for a in range(-7,8) for b in range(-4,5) for c in range(-12,13) for d in range(-3,4)
        if 3*a*a + 11*b*b + c*c + 33*d*d == 144 and a*b == -c*d]
# integer coordinates in the basis (sqrt3/12, sqrt11/12, 1/12, sqrt33/12): the vector itself
half = list({max(v, tuple(-x for x in v)) for v in sols})
Bm = echelon(half); C = [coords(Bm, v) for v in half]
print(len(sols), "unit vectors,", len(half), "directions, module rank", len(Bm))
def cayley_split(n, S, t):
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s in S:
            w = (v + s) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    for c in range(5): cl.append([-X(0, c), -X(t % n, c)])
    return Solver(name="cd19", bootstrap_with=cl).solve()
def cayley_col(n, S):
    X = lambda v, c: 1 + v * 5 + c
    cl = [[X(v, c) for c in range(5)] for v in range(n)]
    for v in range(n):
        for s in S:
            w = (v + s) % n
            if v < w:
                for c in range(5): cl.append([-X(v, c), -X(w, c)])
    return Solver(name="cd19", bootstrap_with=cl).solve()
for n in (5, 10, 15, 20, 25, 30, 40, 50):
    r = len(Bm); adm = 0; col5 = 0; split = None; seen = {}
    for psi in itertools.product(range(n), repeat=r):
        img = [sum(a * b for a, b in zip(psi, c)) % n for c in C]
        if 0 in img: continue
        adm += 1
        S = frozenset(img) | frozenset((-x) % n for x in img)
        if S not in seen: seen[S] = cayley_col(n, S)
        if not seen[S]: continue
        col5 += 1
        for k in range(len(C)):
            t = (5 * img[k]) % n
            if t and cayley_split(n, S, t): split = (psi, k, t); break
        if split: break
    print(f"  Z/{n}: {adm} admissible functionals, {col5} with 5-colourable Cayley graph; "
          f"{'SPLITS a 5M pair ' + str(split) if split else 'none splits a 5e pair'}", flush=True)
