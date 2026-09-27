"""Two pivots, because one cannot reach five colours.

rho(Sa,4) = 7 has a mechanism behind it: the pressure at the pivot is 3, so
N(p) always carries three colours, and {p} together with a small piece of its
circle already forces all four.  One pivot suffices at four colours.

At five it cannot.  Pressure is 2 everywhere, so {p} u N(p) forces only three
colours -- the pivot's own plus the circle's two -- and three is not five.
The natural repair is a SECOND pivot: p and q adjacent, each contributing its
own colour and its circle's two, which could reach five between them.

If {p,q} u N(p) u N(q) is rainbow-forcing on de Grey's G, then rho(G,5) is at
most about 120 rather than 1581, and shrinking it is the next question.  A
core of three at a degree-60 pivot needs 63, so 120 would be within a factor
of two of the target for the first time.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
k = 5


def x(v, c):
    return 1 + v * k + c


base = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        base.append([-x(a, c), -x(b, c)])


def forcing(S):
    """Rainbow-forcing: no proper 5-colouring leaves a colour off S."""
    for c in range(k):
        with Solver(name="cd19",
                    bootstrap_with=base + [[-x(v, c)] for v in S]) as s:
            if s.solve():
                return False, c
    return True, None


deg = g.degrees()
hubs = sorted(range(g.n), key=lambda v: -deg[v])[:40]
print(f"G: {g.n} vertices, max degree {max(deg)}", flush=True)

# one pivot first, to confirm it cannot work
p = hubs[0]
S1 = [p] + sorted(g.adj[p])
t = time.time()
ok, c = forcing(S1)
print(f"  one pivot, {len(S1)} vertices: forcing={ok}"
      + (f" (colour {c} avoidable)" if not ok else "")
      + f"  [{time.time()-t:.0f}s]", flush=True)

# then pairs of adjacent hubs
tried = 0
for p, q in itertools.combinations(hubs, 2):
    if q not in g.adj[p]:
        continue
    S = sorted({p, q} | g.adj[p] | g.adj[q])
    t = time.time()
    ok, c = forcing(S)
    tried += 1
    print(f"  pivots {p},{q}: {len(S)} vertices, forcing={ok}"
          + (f" (colour {c} avoidable)" if not ok else "   *** FORCING ***")
          + f"  [{time.time()-t:.0f}s]", flush=True)
    if ok:
        from hn.forced import ColourRelations, shrink_forcing_set
        rel = ColourRelations(g, k)
        small = shrink_forcing_set(rel, S)
        print(f"    shrunk to {len(small)} vertices, against the 63 a core "
              f"of three needs", flush=True)
        break
    if tried >= 6:
        break
print(f"\n{tried} adjacent hub pairs tested")
