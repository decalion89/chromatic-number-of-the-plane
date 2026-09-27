"""Hunt a forcing 63-set directly, by local search on sets of that exact size.

The lower-bound side cannot reach 63 by exhaustion -- the trajectory says
10^11 rounds -- and the structural route is closed.  What is left is to look
for the object itself: sixty-three vertices that every proper 5-colouring
paints with all five colours.

Local search fits the shape of the problem.  Hold |S| = 63.  Ask for a proper
colouring whose colour-0 class escapes S; if none exists S is FORCING and the
hunt is over.  If one exists, that class names vertices S is missing, so swap
one in -- and swap out whichever member of S has been least useful, measured
by how many of the classes seen so far it is the only cover for.

Unlike the hitting-set loop this never grows the set, so it searches the space
of 63-sets rather than climbing towards it, and every step is one SAT call.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.solvers import Solver
from hn import degrey

g = degrey.build_G()
n, k, SIZE = g.n, 5, 63


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
solver = Solver(name="cd19", bootstrap_with=cls)

rng = random.Random(11)
deg = g.degrees()
# start spread out: highest degree first, then random
S = sorted(rng.sample(range(n), SIZE))
hits = {v: 0 for v in range(n)}       # how many seen classes each vertex cuts
seen, t0, best_streak, streak = 0, time.time(), 0, 0
try:
    for rnd in range(2000000):
        if not solver.solve(assumptions=[-x(v, 0) for v in S]):
            print(f"  *** FORCING SET OF {len(S)} FOUND at round {rnd} ***",
                  flush=True)
            print(f"  {sorted(S)}", flush=True)
            break
        m = set(solver.get_model())
        C = [v for v in range(n) if x(v, 0) in m]
        seen += 1
        streak += 1
        for v in C:
            hits[v] += 1
        # swap: bring in a vertex of the escaping class that cuts the most
        incoming = max(C, key=lambda v: hits[v])
        if incoming in S:
            incoming = rng.choice(C)
        outgoing = min(S, key=lambda v: hits[v] + rng.random())
        S = sorted(set(S) - {outgoing} | {incoming})
        while len(S) < SIZE:
            S = sorted(set(S) | {rng.randrange(n)})
        if streak > best_streak:
            best_streak = streak
        if rnd and rnd % 20000 == 0:
            print(f"  round {rnd}: {seen} escapes seen, longest run without "
                  f"one {best_streak}  [{time.time()-t0:.0f}s]", flush=True)
            streak = 0
finally:
    solver.delete()
