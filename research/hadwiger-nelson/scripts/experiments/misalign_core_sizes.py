"""Does misalignment block cores larger than three?

Pressure(p) >= k - r is the whole constraint on core size, so with pressure 2
at five colours every r >= 3 is permitted -- and a LARGER core is easier to
have, not harder, since it asks less of the graph. Counting blocks none of
them (alpha >= m/3 per leg forces r <= 2), so the only question is whether
misalignment reaches past three.

Same model as before. Copies are (k, e) with e = 0 a rotation by 2*pi*k/N and
e = 1 a reflection; leg L conflicts within an orbit when k - k' = +-t_L and
across the orbits when k - k' = +-t_L - s_L, the shift s_L = 2*phi_L being the
leg's own angle. An escape is an assignment of every copy to a leg with no
conflicting pair sharing a leg.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from math import cos, pi
from pysat.solvers import Solver

random.seed(17)


def blocks(N, ts, ss):
    r = len(ts)
    def x(k, e, L):
        return 1 + ((k * 2 + e) * r + L)
    cls = [[x(k, e, L) for L in range(r)] for k in range(N) for e in (0, 1)]
    for L in range(r):
        t, s = ts[L], ss[L]
        for k in range(N):
            for e in (0, 1):
                for d in (t, -t):
                    w = (k - d) % N
                    if (w, e) != (k, e):
                        cls.append([-x(k, e, L), -x(w, e, L)])
                for d in (t - s, -t - s):
                    cls.append([-x(k, 0, L), -x((k - d) % N, 1, L)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        return not sv.solve()


t0 = time.time()
for r in (3, 4, 5, 6):
    found, tried = [], 0
    for N in range(3, 13):
        pool = list(range(1, N // 2 + 1))
        for _ in range(4000):
            tried += 1
            ts = [random.choice(pool) for _ in range(r)]
            ss = [random.randrange(N) for _ in range(r)]
            if blocks(N, ts, ss):
                found.append((N, tuple(ts), tuple(ss)))
                if len(found) == 1:
                    d2 = [1.0 / (2 * (1 - cos(2 * pi * t / N))) for t in ts]
                    print(f"  r={r}: FIRST BLOCK at N={N} t={ts} s={ss}, "
                          f"d^2 = {[round(v,5) for v in d2]}  "
                          f"[{time.time()-t0:.0f}s]", flush=True)
                break
        if found:
            break
    print(f"r={r}: {len(found)} blocking configuration(s) found in {tried} "
          f"random draws  [{time.time()-t0:.0f}s]", flush=True)
