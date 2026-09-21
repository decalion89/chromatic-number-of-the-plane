"""Build the hypothesis instead of hunting it.

Neither Sa nor G realises the chords the width-two hub spindle needs, because
neither was built for them.  But rho^m(w) is constructible wherever rho is, so
closing a seed under the ring's own rotation puts every chord into the graph
by fiat.

Stack: S_j = union of rho^t(seed) for t = 0..j, all sharing the hub, since rho
fixes it.  At each level ask two things.  Is S_j still k-colourable -- because
if it ever is not, the stack alone has done the job.  And does S_j carry the
hub disjunction at some even m, over any two points of any ring about the hub
at the m-step chord?  If it does, the lemma finishes it with m more levels.

Both halves of every disjunction are also tested alone, because a disjunction
that only holds when one half already does is the ordinary spindle in disguise.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD as K, rotation_joining, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import (closable_distance, width_two_spindle_separation,
                       forced_same, doubly_usable_ring)
from pysat.solvers import Solver

which = sys.argv[1] if len(sys.argv) > 1 else "Sa"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
LEVELS = int(sys.argv[3]) if len(sys.argv) > 3 else 5
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 6000
MS = (2, 4, 6)
t0 = time.time()
seed = build_Sa(K) if which == "Sa" else build_G(K, as_graph=False)

# the hub: highest degree, which is where the ring structure is richest
b0 = IntBasis.covering(seed)
r0 = b0.rows(seed)
deg = defaultdict(int)
for a, b in fast_edges_complete(b0, r0):
    deg[a] += 1
    deg[b] += 1
hub_i = max(range(len(seed)), key=lambda i: deg[i])
HUB = seed[hub_i]
print(f"{which}: {len(seed)} pts, hub index {hub_i} degree {deg[hub_i]}, k={k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# rings about the hub that carry a rotation at all
dim0, D20 = b0.dim, b0.D * b0.D
d = r0 - r0[hub_i]
sq = b0._field_square(d[:, :dim0]) + b0._field_square(d[:, dim0:])
rat = np.ones(len(sq), dtype=bool)
for m in range(1, dim0):
    rat &= sq[:, m] == 0
ringpop = defaultdict(int)
for off in np.nonzero(rat)[0]:
    v = Fr(int(sq[off, 0]), D20)
    if v:
        ringpop[v] += 1
rings = [D for D, c in sorted(ringpop.items(), key=lambda t: -t[1])
         if c >= 2 and D != 1 and closable_distance(D)]
print(f"rings about the hub with a rotation: "
      f"{[(str(D), ringpop[D]) for D in rings]}  [{time.time()-t0:.0f}s]",
      flush=True)


def analyse(U, D):
    b = IntBasis.covering(U)
    r = b.rows(U)
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(U)
    hub = U.index(HUB)
    dd = r - r[hub]
    s = b._field_square(dd[:, :dm]) + b._field_square(dd[:, dm:])
    ok = np.ones(len(s), dtype=bool)
    for m in range(1, dm):
        ok &= s[:, m] == 0
    ring = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(s[off, 0]), d2)
        if v:
            ring[v].append(int(off))
    Eset = set(E)
    # pairwise distances inside each ring, for the chord test
    cands = []
    for Dr, mem in ring.items():
        if Dr == 1 or len(mem) < 2 or not closable_distance(Dr):
            continue
        want = {}
        for m in MS:
            v = width_two_spindle_separation(Dr, m)
            if v and v > 0:
                want.setdefault(v, m)
        if not want:
            continue
        sub = r[mem]
        for a in range(len(mem) - 1):
            dv = sub[a + 1:] - sub[a]
            sv2 = b._field_square(dv[:, :dm]) + b._field_square(dv[:, dm:])
            g = np.ones(len(sv2), dtype=bool)
            for m in range(1, dm):
                g &= sv2[:, m] == 0
            for off in np.nonzero(g)[0]:
                bb = mem[a + 1 + int(off)]
                val = Fr(int(sv2[off, 0]), d2)
                if val in want:
                    w, wp = mem[a], bb
                    if (min(hub, w), max(hub, w)) not in Eset and \
                       (min(hub, wp), max(hub, wp)) not in Eset:
                        cands.append((Dr, want[val], w, wp))
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sel0 = 1 + n * k
    for idx, (Dr, m, w, wp) in enumerate(cands):
        s_ = sel0 + idx
        for col in range(k):
            cls.append([-s_, -(1 + hub * k + col), -(1 + w * k + col)])
            cls.append([-s_, -(1 + hub * k + col), -(1 + wp * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    colourable = sv.solve()
    hits = []
    if colourable:
        for idx, (Dr, m, w, wp) in enumerate(cands):
            if not sv.solve(assumptions=[sel0 + idx]):
                halves = (forced_same(sv, hub, w, k),
                          forced_same(sv, hub, wp, k))
                hits.append((Dr, m, w, wp, halves))
    sv.delete()
    return n, len(E), colourable, len(cands), hits


for D in rings:
    rot = rotation_joining(D, K)
    f = rot.about(HUB)
    U, seen = list(seed), set(seed)
    cur = list(seed)
    kind = "doubly usable" if doubly_usable_ring(D) else "spindle only"
    print(f"\n--- ring D={D} ({kind}) ---", flush=True)
    for lev in range(1, LEVELS + 1):
        cur = [f(p) for p in cur]
        fresh = [p for p in cur if p not in seen]
        seen.update(fresh)
        U.extend(fresh)
        if len(U) > CAP:
            print(f"  level {lev}: {len(U)} pts exceeds cap, stop", flush=True)
            break
        n, m, colourable, nc, hits = analyse(U, D)
        tag = "COLOURS" if colourable else f"*** NOT {k}-COLOURABLE ***"
        print(f"  level {lev}: {n} pts, {m} edges, {tag}, {nc} chord pairs, "
              f"{len(hits)} disjunctions  [{time.time()-t0:.0f}s]", flush=True)
        if not colourable:
            print(f"  *** chi >= {k+1} from a plain rho-stack, ring {D}, "
                  f"level {lev} ***", flush=True)
            break
        for Dr, mm, w, wp, halves in hits[:5]:
            tg = "GENUINE" if not any(halves) else "one half alone"
            print(f"      disjunction ring {Dr}, m={mm}: {tg} halves={halves}",
                  flush=True)
        if hits:
            break
print(f"\nDONE  [{time.time()-t0:.0f}s]", flush=True)
