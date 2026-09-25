"""De Grey's exhaustion, redone in positions instead of colours.

His lemma came from forcing the centre and its radius-2 hexagon to each of
the 187 available colour partitions and finding that ten survive -- the seven
points take at most two colours and the centre is never alone.  That is the
whole engine, and this project reproduced it earlier in the colour frame.

The same exhaustion has never been done in the position frame, and it is
small enough to do completely: seven points, nine positions, 9^7 = 4.78
million assignments, every one checkable.  If the survivors confine the
seven points to a narrow arc then there is a circular analogue of the bite,
and the spindle's q/2 requirement has something to feed on.  If they spread,
the reason will be visible rather than guessed at.

The configuration is de Grey's own: the origin and the six points at
distance 2 from it that sit a unit apart consecutively -- (+-2, 0) and
(+-1, +-sqrt 3).  Consecutive ring points are adjacent, and every ring point
is at distance 2 from the centre, so the centre is NOT adjacent to them.
The edges are the ring's hexagon only, which is exactly why the colour
lemma needed the bite to add the centre's six edges.  Both versions are run.
"""
import sys, itertools, time
from collections import Counter

p, q = 9, 2
RING = [(i, (i + 1) % 6) for i in range(6)]          # hexagon edges
BITE = [(6, i) for i in range(6)]                     # centre joined to ring


def circ(a, b):
    d = abs(a - b) % p
    return min(d, p - d)


def survivors(edges):
    out = []
    for asg in itertools.product(range(p), repeat=7):
        if all(circ(asg[u], asg[v]) >= q for u, v in edges):
            out.append(asg)
    return out


def describe(name, edges):
    t0 = time.time()
    S = survivors(edges)
    print(f"\n{name}: {len(S)} of {p**7} assignments survive"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not S:
        return
    # how wide an arc does each survivor need?
    def arcspan(a):
        used = sorted(set(a))
        best = p
        for i in range(len(used)):
            # gap after used[i]
            g = (used[(i + 1) % len(used)] - used[i]) % p
            best = min(best, p - g)
        return best
    spans = Counter(arcspan(a) for a in S)
    print(f"   arc span of the seven points: {dict(sorted(spans.items()))}",
          flush=True)
    print(f"   narrowest arc any assignment needs: {min(spans)} of {p}"
          f"   (forced independent would need < 2q = {2*q})", flush=True)
    # how many distinct positions do they use?
    dist = Counter(len(set(a)) for a in S)
    print(f"   distinct positions used: {dict(sorted(dist.items()))}",
          flush=True)
    # is the centre ever forced close to a ring point?
    mind = Counter(min(circ(a[6], a[i]) for i in range(6)) for a in S)
    print(f"   centre's closest distance to a ring point: "
          f"{dict(sorted(mind.items()))}", flush=True)
    # antipodal ring pairs: how close can they be forced?
    anti = Counter(max(circ(a[i], a[i + 3]) for i in range(3)) for a in S)
    print(f"   widest antipodal ring pair: {dict(sorted(anti.items()))}",
          flush=True)


print(f"K({p}/{q}) = {p/q}: adjacency needs circular distance >= {q}",
      flush=True)
describe("ring only (hexagon)", RING)
describe("ring + bite (centre joined to all six)", RING + BITE)
print("\nDONE", flush=True)
