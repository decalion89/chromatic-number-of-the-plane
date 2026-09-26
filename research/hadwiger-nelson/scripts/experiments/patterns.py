"""The colour patterns a circle can carry.

A spindle glues H to rho(H) along a circle R.  Whether the glued graph is
k-colourable is decided entirely by *which colour patterns R can carry* inside
H -- nothing else about H matters at the interface.  So the object to compute
is not a cap (one number) but the whole set

    P_k(R, H) = { colouring restricted to R : over all proper k-colourings of H }

up to permutation of the colours.  The cap is max |image| over P; the set is
strictly more information, and it is what decides every possible gluing.

Enumeration: solve, read the pattern on R, then block *all k! relabellings* of
it at once, so each solve returns a genuinely new pattern class and the loop
runs once per class rather than once per raw colouring.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools, pickle
from fractions import Fraction as Fr
from collections import Counter

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 5, 7, 11))
ORI = Point(F.zero(), F.zero())

def circle(pts, centre, d2):
    return [i for i, p in enumerate(pts) if (p - centre).norm2() == d2]

def patterns(g, R, k, cap_classes=400000, verbose=True):
    n = g.n
    X = lambda v, c: 1 + v * k + c
    cnf = [[X(v, c) for c in range(k)] for v in range(n)]
    for u, v in g.edges():
        for c in range(k):
            cnf.append([-X(u, c), -X(v, c)])
    for v in R:                       # at-most-one, only where the pattern is read
        for c in range(k):
            for d in range(c + 1, k):
                cnf.append([-X(v, c), -X(v, d)])
    s = Solver(name="m22", bootstrap_with=cnf)
    perms = list(itertools.permutations(range(k)))
    out, t0 = [], time.time()
    while s.solve():
        pos = set(l for l in s.get_model() if l > 0)
        pat = []
        for v in R:
            pat.append(next(c for c in range(k) if X(v, c) in pos))
        # canonical form: relabel by first occurrence
        seen, ren = {}, []
        for c in pat:
            if c not in seen:
                seen[c] = len(seen)
            ren.append(seen[c])
        out.append(tuple(ren))
        for sg in perms:
            s.add_clause([-X(v, sg[pat[i]]) for i, v in enumerate(R)])
        if len(out) >= cap_classes:
            if verbose: print(f"    stopped at cap {cap_classes}", flush=True)
            s.delete(); return out, False
    s.delete()
    if verbose: print(f"    {len(out)} classes in {time.time()-t0:.1f}s", flush=True)
    return out, True

if __name__ == "__main__":
    which = sys.argv[1]
    k = int(sys.argv[2])
    d2 = Fr(sys.argv[3])
    pts = {"Sa": build_Sa, "Y": build_Y}[which](F)
    g = build_graph(pts)
    R = circle(pts, ORI, d2)
    print(f"{which}: n={g.n}  circle r^2={d2}: {len(R)} points  k={k}", flush=True)
    if not R:
        sys.exit("empty circle")
    pat, complete = patterns(g, R, k)
    sizes = Counter(len(set(p)) for p in pat)
    print(f"  palette sizes over the {len(pat)} classes: "
          + "  ".join(f"{s}->{c}" for s, c in sorted(sizes.items())))
    print(f"  CAP (max colours the circle can show) = {max(len(set(p)) for p in pat)}"
          + ("" if complete else "   [lower bound: enumeration truncated]"))
    # the distinguishing structure: which pairs of the circle can ever agree
    m = len(R)
    agree = [[False]*m for _ in range(m)]
    for p in pat:
        for i in range(m):
            for j in range(i+1, m):
                if p[i] == p[j]: agree[i][j] = True
    forced = [(i, j) for i in range(m) for j in range(i+1, m) if not agree[i][j]]
    print(f"  pairs of the circle FORCED to differ: {len(forced)} of {m*(m-1)//2}")
    pickle.dump({"which":which,"k":k,"d2":str(d2),"R":R,"patterns":pat,
                 "complete":complete}, open(f"/tmp/hn/pat_{which}_{k}_{d2.numerator}_{d2.denominator}.pkl","wb"))
