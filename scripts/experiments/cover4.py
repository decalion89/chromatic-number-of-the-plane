"""The orbit cover that blocks at 2, 3, 4 AND 5.

Untargeted growth reached 25 orbits blocking at 5 and was still short at 4,
with rounds costing nine minutes.  Worth knowing the target first: which
orbits, and how many, block at every modulus up to five.

CEGAR again, now over four moduli at once.  Hold a set of chosen orbits; for
each n in turn ask has_homomorphism for a phi the chosen set misses; when one
turns up, add the orbit that kills the most of the escaping functionals
collected so far.  The loop ends when no n produces an escape, which is
blocking at all four, exactly.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
from fractions import Fraction as Fr
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
exec(open("/tmp/hn/orbblock.py")
     .read().split("reps = sorted(set(orbit.values()))")[0])
from hn.homcol import has_homomorphism

t0 = time.time()
byorb = {}
for u in steps_all:
    byorb.setdefault(orbit[u], []).append(u)
reps = sorted(byorb)


def ints(elts):
    raw = [flat(u) for u in elts]
    dn = 1
    for v in raw:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    out = set()
    for v in raw:
        w = tuple(int(q * dn) for q in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w) if g > 1 else w)
    return sorted(out)


orbvec = {r: ints(byorb[r]) for r in reps}
print(f"{len(reps)} orbits  [{time.time()-t0:.0f}s]", flush=True)

escapes = []          # (n, phi) pairs the chosen orbits must kill
chosen = []


def kills(r, n, c):
    return any(sum(x * y for x, y in zip(c, d)) % n == 0 for d in orbvec[r])


rounds = 0
while True:
    acc = []
    for r in chosen:
        acc.extend(orbvec[r])
    acc = sorted(set(acc))
    bad = None
    for n in (5, 4, 3, 2):
        if not acc:
            bad = (n, None)
            break
        ph, _ = has_homomorphism(acc, n)
        if ph is not None:
            bad = (n, tuple(ph))
            break
    rounds += 1
    if bad is None:
        print(f"  round {rounds}: {len(chosen)} orbits, {len(acc)} directions "
              f"-- BLOCKS AT 2, 3, 4 AND 5  [{time.time()-t0:.0f}s]",
              flush=True)
        import pickle
        with open("/tmp/hn/"
                  "orbs4.pkl", "wb") as fh:
            pickle.dump(chosen, fh)
        break
    n, ph = bad
    if ph is not None:
        escapes.append((n, ph))
    best, bk = None, -1
    for r in reps:
        if r in chosen:
            continue
        k = sum(1 for (nn, c) in escapes if kills(r, nn, c))
        if k > bk:
            best, bk = r, k
    chosen.append(best)
    if rounds % 5 == 0 or len(chosen) < 6:
        print(f"  round {rounds}: {len(chosen)} orbits, escape at n = {n}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
