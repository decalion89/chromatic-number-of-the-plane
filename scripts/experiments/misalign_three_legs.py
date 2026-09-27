"""Three legs, with reflections: the freedom the circulant model left out.

Rotations alone put every leg's conflict graph on the same circulant, and a
search over all N <= 30 and all angle triples found no blocking triple at all
-- against 298 blocking pairs for a core of two, so the search works.

Reflections change the picture, and they are free: a reflection fixing the
pivot shares the rotations' discriminant, so any field naming one names the
other. Index a copy by (k, e) with e = 0 a rotation by 2*pi*k/N and e = 1 a
reflection. Leg L sits at its own angle phi_L, so its image under (k, 0) is at
phi_L + 2*pi*k/N and under (k, 1) at -phi_L + 2*pi*k/N. Within one orbit the
conflict is the old circulant, k - k' = +-t_L; ACROSS the orbits it picks up
the leg's own angle,

    k - k'  =  +-t_L - s_L,     s_L = 2*phi_L read in N-ths of a turn,

and s_L differs from leg to leg because the legs sit at different angles.
That per-leg shift is exactly the misalignment counting cannot see, and it is
what made the two-leg block work: the six copies there are three rotations
and the same three composed with a second isometry.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from itertools import combinations_with_replacement, product
from math import cos, pi
from pysat.solvers import Solver

t0, found, tried = time.time(), [], 0
for N in range(3, 16):
    ts = list(range(1, N // 2 + 1))
    for trio in combinations_with_replacement(ts, 3):
        for shifts in product(range(N), repeat=3):
            tried += 1
            def x(k, e, L):
                return 1 + ((k * 2 + e) * 3 + L)
            cls = [[x(k, e, L) for L in range(3)]
                   for k in range(N) for e in (0, 1)]
            for L in range(3):
                t, s = trio[L], shifts[L]
                for k in range(N):
                    for e in (0, 1):
                        for d in (t, -t):
                            w = (k - d) % N
                            if (w, e) != (k, e):
                                cls.append([-x(k, e, L), -x(w, e, L)])
                        for d in (t - s, -t - s):
                            w = (k - d) % N
                            cls.append([-x(k, 0, L), -x(w, 1, L)])
            with Solver(name="cd19", bootstrap_with=cls) as s_:
                if not s_.solve():
                    d2 = [1.0 / (2 * (1 - cos(2 * pi * t / N))) for t in trio]
                    found.append((N, trio, shifts))
                    print(f"  *** N={N} t={trio} shifts={shifts} BLOCKS "
                          f"-- d^2 = {[round(v,5) for v in d2]}  "
                          f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"    ... N={N} done, {tried} configurations, {len(found)} block  "
          f"[{time.time()-t0:.0f}s]", flush=True)
print(f"blocking triples with reflections: {len(found)}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
