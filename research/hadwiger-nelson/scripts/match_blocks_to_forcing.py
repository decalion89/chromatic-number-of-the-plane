"""Cross the blocking configurations against the forceable distances.

Blocking says: IF the three legs are a core of the pivot, THEN the union of
the copies has no colouring. Supplying the core is the other half, and at
three colours it has a clean form -- a core of three needs pressure k-3, so
the pivot has no neighbours at all and the three legs must be pairwise unable
to share. On the triangular lattice that is decided: the 3-colouring is the
residue modulo (1 - omega), of norm 3, so two points are forced to differ
exactly when the norm of their difference is prime to 3. Measured: 1, 4 and
7, every pair of them, and nothing else.

So a blocking configuration is usable when its three legs sit pairwise at
such distances AND that triangle embeds in the lattice. The legs' radii are
fixed absolutely by the blocking pattern -- d_L = 1/(2 sin(pi t_L/N)) -- so
the leg triangle has fixed side lengths and there is nothing to scale.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from itertools import combinations_with_replacement, product
from math import cos, pi, sin
from pysat.solvers import Solver

# Eisenstein norms prime to 3, which are the forceable squared distances
def norms(limit):
    out = set()
    for a in range(-40, 41):
        for b in range(-40, 41):
            n = a * a + a * b + b * b
            if 0 < n <= limit and n % 3:
                out.add(n)
    return out

FORCEABLE = norms(60)
print(f"forceable squared distances (Eisenstein norms prime to 3): "
      f"{sorted(FORCEABLE)[:14]} ...", flush=True)

def blocks(N, trio, shifts):
    def x(k, e, L):
        return 1 + ((k * 2 + e) * 3 + L)
    cls = [[x(k, e, L) for L in range(3)] for k in range(N) for e in (0, 1)]
    for L in range(3):
        t, s = trio[L], shifts[L]
        for k in range(N):
            for e in (0, 1):
                for d in (t, -t):
                    w = (k - d) % N
                    if (w, e) != (k, e):
                        cls.append([-x(k, e, L), -x(w, e, L)])
                for d in (t - s, -t - s):
                    cls.append([-x(k, 0, L), -x((k - d) % N, 1, L)])
    with Solver(name="cd19", bootstrap_with=cls) as s_:
        return not s_.solve()

t0, hits, tested = time.time(), [], 0
for N in range(3, 37):
    for trio in combinations_with_replacement(range(1, N // 2 + 1), 3):
        d = [1.0 / (2 * sin(pi * t / N)) for t in trio]
        for shifts in product(range(N), repeat=3):
            for flip in product((0, 1), repeat=3):
                phi = [pi * shifts[L] / N + pi * flip[L] for L in range(3)]
                pair = [d[i]**2 + d[j]**2 - 2*d[i]*d[j]*cos(phi[i]-phi[j])
                        for i, j in ((0, 1), (0, 2), (1, 2))]
                if any(p < 0.5 for p in pair):
                    continue                      # coincident or too close
                r = [round(p) for p in pair]
                if any(abs(p - q) > 1e-9 for p, q in zip(pair, r)):
                    continue
                if not all(q in FORCEABLE for q in r):
                    continue
                tested += 1
                if blocks(N, trio, shifts):
                    hits.append((N, trio, shifts, flip, r,
                                 [round(v*v, 6) for v in d]))
                    print(f"  *** N={N} t={trio} s={shifts} flip={flip}: legs "
                          f"pairwise {r}, radii^2 {[round(v*v,5) for v in d]}"
                          f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"    ... N={N}, {tested} candidates with forceable legs, "
          f"{len(hits)} block  [{time.time()-t0:.0f}s]", flush=True)
print(f"usable configurations: {len(hits)}  [{time.time()-t0:.0f}s]", flush=True)
