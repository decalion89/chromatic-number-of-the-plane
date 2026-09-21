"""Bite Gp on every ring the field allows, including the wide ones.

Gp is G closed under the order-twelve group about its own centre, the pivot
(-2, 0): 13873 points, 5-chromatic, and symmetric, which is de Grey's position
at four colours one level up.  Rings about the pivot at every radius, not just
the ones narrower than a unit distance, give twelve turns that exist exactly
in the field -- among them cosine 31/32 with sine 3*sqrt7/32 on the ring of
radius four, which is de Grey's own second-level angle, and which every ring
scan in this work had filtered away.

A bite is the union of the graph with its image under such a turn.  Each point
of the bitten ring lands a unit from itself, so the two copies are stitched
together along that ring, and a disjunction about the ring becomes a statement
about a named pair.  Eleven of the twelve are turns nobody has applied.

Cross edges are counted, because a bite that produces none has not bitten.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
Gp = pickle.load(open(SC + "pivG.pkl", "rb"))
PIV = Point(K.rational(-2), K.zero())
print(f"Gp: {len(Gp)} points  [{time.time()-t0:.0f}s]", flush=True)

# (ring D, points on it, cos, sin as (rational, radicand))
BITES = [
    (Fr(16), 12, Fr(31, 32), (Fr(3, 32), 7)),      # de Grey's own level-2 turn
    (Fr(9), 12, Fr(17, 18), (Fr(1, 18), 35)),
    (Fr(4), 36, Fr(7, 8), (Fr(1, 8), 15)),         # de Grey's Sb turn
    (Fr(7), 24, Fr(13, 14), (Fr(3, 14), 3)),
    (Fr(17, 2), 24, Fr(16, 17), (Fr(1, 17), 33)),
    (Fr(20, 3), 24, Fr(37, 40), (Fr(1, 40), 231)),
    (Fr(31, 3), 24, Fr(59, 62), (Fr(11, 62), 3)),
    (Fr(13, 3), 24, Fr(23, 26), (Fr(7, 26), 3)),
    (Fr(7, 3), 24, Fr(11, 14), (Fr(5, 14), 3)),
    (Fr(3, 2), 24, Fr(2, 3), (Fr(1, 3), 5)),
    (Fr(3), 12, Fr(5, 6), (Fr(1, 6), 11)),         # the Moser hinge
]
base = set(Gp)
for D, cnt, ct, (sc, rad) in BITES:
    c = K.rational(ct)
    s = K.sqrt(rad) * K.rational(sc) if rad != 1 else K.rational(sc)
    assert c * c + s * s == K.rational(1), (D, "not a rotation")
    turn = Rotation(c, s).about(PIV)
    seen, P = set(), []
    for p in Gp:
        for q in (p, turn(p)):
            if q not in seen:
                seen.add(q)
                P.append(q)
    fresh = len(P) - len(Gp)
    b = IntBasis.covering(P)
    r = b.rows(P)
    hr = b.overflow_headroom(r)
    if hr >= 1.0:
        print(f"D={D}: SKIPPED, overflow headroom {hr:.2f}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        continue
    E = sorted(set((min(a, cc), max(a, cc))
                   for a, cc in fast_edges_complete(b, r)))
    n = len(P)
    old = set(range(len(Gp)))
    cross = sum(1 for a, cc in E if (a in old) != (cc in old))
    cls = [[1 + v * k + col for col in range(k)] for v in range(n)]
    for a, cc in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + cc * k + col)])
    for col in range(1, k):
        cls.append([-(1 + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    verdict = (f"{k}-colourable" if ok
               else f"*** NOT {k}-COLOURABLE -- chi > {k} ***")
    print(f"D={D} (rho={float(D)**.5:.3f}, {cnt} ring pts): {n} points "
          f"({fresh} new), {len(E)} edges, {cross} cross -> {verdict}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        with open(SC + f"WITNESS_bite_{str(D).replace('/','_')}.pkl",
                  "wb") as f:
            pickle.dump(P, f)
        print("   witness written", flush=True)
print("DONE", flush=True)
