"""How few orbits does it take to block in this field?

The unguarded chain needed 25 zeta_6-orbits, 150 directions.  But that was
whatever the growth happened to pick up, not a minimum.  A chain is more
likely to be critical the shorter it is -- every rhombus has to matter -- so
it is worth knowing the target before growing towards it.

Two numbers: the minimum blocking subset of the 450 directions, and the
minimum number of whole ORBITS, since a chain adds orbits rather than
individual directions.  The second is a set cover over the hyperplanes each
orbit kills, which is too large to enumerate at rank 12 -- so it is done
greedily against sampled functionals and then verified exactly with
has_homomorphism.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, random, time
from fractions import Fraction as Fr
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
exec(open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/orbblock.py")
     .read().split("reps = sorted(set(orbit.values()))")[0])
from hn.homcol import has_homomorphism

t0 = time.time()
byorb = {}
for u in steps_all:
    byorb.setdefault(orbit[u], []).append(u)
reps = sorted(byorb)
print(f"{len(reps)} orbits  [{time.time()-t0:.0f}s]", flush=True)


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


# Greedy set cover, driven by CEGAR rather than by a fixed sample.  A random
# sample of 30000 functionals was exhausted by 11 orbits while the exact check
# still found one escaping -- there are 61 million of them at rank 12, so a
# sample says nothing on its own.  Instead: cover the sample greedily, then
# ASK has_homomorphism for a functional the chosen orbits miss, add it, and
# cover again.  The loop ends only when no escaping functional exists, which
# is blocking, exactly.
import random

rng = random.Random(9)
D = 12
sample = []
while len(sample) < 4000:
    c = tuple(rng.randrange(5) for _ in range(D))
    if any(c):
        sample.append(c)
orbvec = {r: ints(byorb[r]) for r in reps}


def kills(r, c):
    return any(sum(x * y for x, y in zip(c, d)) % 5 == 0 for d in orbvec[r])


chosen, rounds = [], 0
while True:
    alive = [c for c in sample
             if not any(kills(r, c) for r in chosen)]
    while alive:
        best, bestkill = None, -1
        for r in reps:
            if r in chosen:
                continue
            k = sum(1 for c in alive if kills(r, c))
            if k > bestkill:
                best, bestkill = r, k
        chosen.append(best)
        alive = [c for c in alive if not kills(best, c)]
    acc = []
    for r in chosen:
        acc.extend(orbvec[r])
    acc = sorted(set(acc))
    phi, _ = has_homomorphism(acc, 5)
    rounds += 1
    print(f"  round {rounds}: {len(chosen):2d} orbits, {len(acc):3d} "
          f"directions, "
          + ("*** BLOCKS ***" if phi is None else "one more escapes")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
    if phi is None:
        print(f"  a blocking set of {len(chosen)} orbits / {len(acc)} "
              f"directions -- the chain has to reach that, not 25 orbits",
              flush=True)
        import pickle
        with open("/tmp/claude-0/-home-user-darwin-50/"
                  "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/orbs.pkl",
                  "wb") as fh:
            pickle.dump(chosen, fh)
        break
    sample.append(tuple(phi))
