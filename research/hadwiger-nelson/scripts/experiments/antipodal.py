"""The weaker target: some antipodal pair monochromatic, not a palette cap.

His lemma, measured, is a palette cap -- centre and ring held to two colours
out of four -- and that is what every scan here has hunted.  But the cap is
not what the bite consumes.  What the machine needs is only this: in every
colouring, at least one antipodal pair of the ring carries one colour twice.
The cap implies it by pigeonhole and is strictly stronger, so hunting the cap
has been hunting the harder of the two objects all along.

The weaker statement has its own test, and a single call settles it.  Ask for
a colouring in which EVERY antipodal pair of the ring is bichromatic; if none
exists, then in every colouring some pair is monochromatic, which is the
disjunction, and the ring is biteable if the field holds its turn.

A ring about a centre has antipodal structure when p and its reflection
2c - p are both on it, which the order-twelve closures guarantee.  The
statement travels upward like the others: adding vertices removes colourings,
so a disjunction proved on this graph holds in every supergraph.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SEED = sys.argv[1] if len(sys.argv) > 1 else "pivG.pkl"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
CX = Fr(sys.argv[3]) if len(sys.argv) > 3 else None
CY = Fr(sys.argv[4]) if len(sys.argv) > 4 else None
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/hn/")
if SEED.endswith(".pkl"):
    P = pickle.load(open(SC + SEED, "rb"))
else:
    P = {"Sa": build_Sa, "Y": build_Y,
         "G": lambda f: build_G(f, as_graph=False)}[SEED](K)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(int)
for a, c in E:
    adj[a] += 1
    adj[c] += 1
print(f"{SEED}: {n} points, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]",
      flush=True)
base = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        base.append([-(1 + a * k + col), -(1 + c * k + col)])


def biteable(D):
    if D < Fr(1, 4):
        return None
    ct = 1 - Fr(1, 2) / D
    s2 = 1 - ct * ct
    for rr in CLASSES:
        q = s2 / rr
        a, c = q.numerator, q.denominator
        ra, rc = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc * rc == c:
            return f"cos={ct}, sin={Fr(ra,rc)}" + (f"*sqrt{rr}" if rr != 1
                                                   else "")
    return None


idx = {p: i for i, p in enumerate(P)}
dm, D2 = b.dim, b.D * b.D
if CX is not None:
    centres = [Point(K.rational(CX), K.rational(CY))]
else:
    centres = [P[v] for v in sorted(range(n), key=lambda v: -adj[v])[:60]]
found = tested = 0
for C in centres:
    rc = b.rows([C])[0]
    d = r - rc
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        Dv = Fr(int(sq[off, 0]), D2)
        if Dv > 0:
            grp[Dv].append(int(off))
    for Dv, mem in sorted(grp.items(), key=lambda t: -len(t[1])):
        if len(mem) < 4:
            continue
        bt = biteable(Dv)
        if not bt:
            continue
        pairs = []
        seen = set()
        for v in mem:
            if v in seen:
                continue
            q = Point(C.x + C.x - P[v].x, C.y + C.y - P[v].y)
            w = idx.get(q)
            if w is not None and w in mem and w != v:
                pairs.append((v, w))
                seen.add(v)
                seen.add(w)
        if len(pairs) < 2:
            continue
        tested += 1
        cls = list(base)
        for v, w in pairs:              # every antipodal pair bichromatic
            for col in range(k):
                cls.append([-(1 + v * k + col), -(1 + w * k + col)])
        s = Solver(name="cd15", bootstrap_with=cls)
        out = s.solve()
        s.delete()
        if not out:
            found += 1
            print(f"*** DISJUNCTION: centre ({float(C.x):.3f},"
                  f"{float(C.y):.3f}), D={Dv}, {len(mem)} ring pts, "
                  f"{len(pairs)} antipodal pairs -- in every {k}-colouring "
                  f"one of them is monochromatic.  Bite: {bt}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            with open(SC + f"disj_{SEED}_{k}_{str(Dv).replace('/','_')}.pkl",
                      "wb") as f:
                pickle.dump((C, [P[v] for v in mem], pairs), f)
    if tested and tested % 40 == 0:
        print(f"   {tested} biteable antipodal rings tested, {found} "
              f"disjunctions  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tested} biteable rings with antipodal pairs, {found} disjunctions"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
