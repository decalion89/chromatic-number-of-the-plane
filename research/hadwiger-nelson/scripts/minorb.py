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
import sys, itertools, random, time
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
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


# Greedy over sampled functionals: pick the orbit killing the most survivors.
rng = random.Random(9)
D = 12
sample = []
while len(sample) < 30000:
    c = tuple(rng.randrange(5) for _ in range(D))
    if any(c):
        sample.append(c)
orbvec = {r: ints(byorb[r]) for r in reps}
alive = list(range(len(sample)))
chosen = []
while alive:
    best, bestkill = None, -1
    for r in reps:
        if r in chosen:
            continue
        vs = orbvec[r]
        kill = 0
        for i in alive:
            c = sample[i]
            if any(sum(x * y for x, y in zip(c, d)) % 5 == 0 for d in vs):
                kill += 1
        if kill > bestkill:
            best, bestkill = r, kill
    chosen.append(best)
    vs = orbvec[best]
    alive = [i for i in alive
             if not any(sum(x * y for x, y in zip(sample[i], d)) % 5 == 0
                        for d in vs)]
    acc = []
    for r in chosen:
        acc.extend(orbvec[r])
    phi, _ = has_homomorphism(sorted(set(acc)), 5)
    print(f"  {len(chosen):2d} orbits, {len(set(acc)):3d} directions, "
          f"{len(alive):5d} sampled functionals still escaping, "
          + ("*** BLOCKS ***" if phi is None else "coset colouring")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
    if phi is None:
        print(f"  minimum found greedily: {len(chosen)} orbits", flush=True)
        break
