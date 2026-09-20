"""Test the pair the construction points at, not three hundred million others.

A ring D of the symmetric core is usable twice, and both uses need a radical.
rho, the rotation that makes the ring produce cross edges, needs sqrt(4D-1).
The antipodal pair on that ring sits at squared distance 4D, so spindling it
afterwards needs sqrt(16D-1).  Requiring both at once is a savage filter.

Sa has five rational closable rings and exactly ONE of them -- D = 4 -- passes
both.  That is the ring de Grey used, and the pair he forces, (2,0) and
(-2,0), is its antipodal pair.  The criterion retrodicts his choice uniquely,
which is the reason to trust it going up a level.

G* has twelve, and two of them pass: D = 4 again, and D = 17/2, which is a
ring the four-colour construction never had -- rho by sqrt(33), spindle by
sqrt(135) = 3 sqrt(15).  So there are two candidates, not a search space, and
the pairs to interrogate are the ones lying on the ring itself: a few dozen,
seconds of SAT each, instead of an enumeration that cannot be afforded.
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()


def probe(name, core, pivot, D, k):
    base = rotation_joining(D, F)
    for tag, turns in (("two-fold", (base.about(pivot),)),
                       ("three-fold", (base.about(pivot),
                                       base.inverse().about(pivot)))):
        allp, seen = list(core), set(core)
        for p in core:
            for tr in turns:
                q = tr(p)
                if q not in seen:
                    seen.add(q)
                    allp.append(q)
        n = len(allp)
        E = sorted((min(a, b), max(a, b))
                   for a, b in build_graph(allp).edges())
        ring = [i for i, p in enumerate(allp)
                if (p.x - pivot.x) * (p.x - pivot.x)
                + (p.y - pivot.y) * (p.y - pivot.y) == F.rational(D)]
        cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        if not sv.solve():
            print(f"  {name} D={D} {tag}: {n} points, {len(E)} edges -- NOT "
                  f"{k}-COLOURABLE  [{time.time()-t0:.0f}s]", flush=True)
            with open(SC + f"ringpair_six_{k}.pkl", "wb") as fh:
                pickle.dump((name, str(D), tag,
                             [(str(p.x), str(p.y)) for p in allp], E), fh)
            sv.delete()
            return True
        anti, hits = [], []
        for a in range(len(ring)):
            for b in range(a + 1, len(ring)):
                i, j = ring[a], ring[b]
                pa, pb = allp[i], allp[j]
                if pa.x + pb.x == pivot.x + pivot.x and \
                   pa.y + pb.y == pivot.y + pivot.y:
                    anti.append((i, j))
        # Antipodal pairs only.  Proving one pair forced costs six seconds
        # on an 800-point graph at four colours and far more at five, so the
        # 66 pairs of a 12-point ring are already eleven minutes and the 2500
        # of G*'s are unaffordable.  The criterion does not ask for them: the
        # pair it aims at is the antipodal one, which is what de Grey forces.
        for nth, (i, j) in enumerate(anti):
            t1 = time.time()
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)]):
                hits.append((i, j))
            print(f"      antipodal {nth+1}/{len(anti)}: "
                  f"{'FORCED' if (i, j) in hits else 'free'} "
                  f"[{time.time()-t1:.0f}s]", flush=True)
        sv.delete()
        star = " *** FORCED ***" if hits else ""
        print(f"  {name} D={D} {tag}: {n} points, {len(E)} edges, "
              f"{len(ring)} on the ring, {len(anti)} antipodal, "
              f"{len(hits)} FORCED{star}  [{time.time()-t0:.0f}s]", flush=True)
        if hits:
            with open(SC + f"ringpair_forced_{k}.pkl", "wb") as fh:
                pickle.dump((name, str(D), tag,
                             [(str(p.x), str(p.y)) for p in allp], E, hits,
                             anti), fh)
            return True
    return False


# The control first: the criterion must reproduce de Grey's own forced pair.
ORIG = Point(F.zero(), F.zero())
print("control -- Sa at four colours, the ring the criterion picks:",
      flush=True)
probe("Sa", build_Sa(F), ORIG, Fr(4), 4)

PIV = Point(F.rational(-2), F.zero())
rot = _rot60(F).about(PIV)
G = build_G(F, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
print(f"G*: {len(Gs)} points  [{time.time()-t0:.0f}s]", flush=True)
for D in (Fr(4), Fr(17, 2)):
    if probe("G*", Gs, PIV, D, 5):
        break
