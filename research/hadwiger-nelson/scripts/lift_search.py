"""The lift, with the right criterion: list-colourability, not chromatic number.

de Grey reaches pressure 3 at four colours while his confined set is only
2-chromatic in the worst orientation, so "make the confined subgraph
4-chromatic" was never the condition -- it is a sufficient special case. The
condition is that the auxiliary graph, with the lists the squeezed circle
hands it (k-2 where it sees both circle colours, k-1 where it sees one),
admits no proper colouring, in EVERY orientation.

Angles come from de Grey's own family, 0, +-theta/2, +-(60 - theta) with
cos theta = 5/6 the Moser rotation, plus the solved within-pair values.
Auxiliaries are Minkowski sums of the hexagons, which is all of them.
"""
import sys, time, math, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver

TWO = 2 * math.pi
TH = math.acos(5 / 6)
ANGLES = [0.0, TH / 2, TH, math.radians(60) - TH, math.radians(60) - TH / 2,
          math.radians(28.95502), math.radians(31.04498)]
H = int(sys.argv[1]) if len(sys.argv) > 1 else 3
K = int(sys.argv[2]) if len(sys.argv) > 2 else 5


def build(phis):
    hexs = [[(math.cos(p + TWO * k / 6), math.sin(p + TWO * k / 6))
             for k in range(6)] for p in phis]
    aux = []
    for a, b in itertools.combinations(range(len(phis)), 2):
        for i, u in enumerate(hexs[a]):
            for j, v in enumerate(hexs[b]):
                q = (u[0] + v[0], u[1] + v[1])
                if q[0] ** 2 + q[1] ** 2 > 1e-12:
                    aux.append((q, a, b, i % 2, j % 2))
    return hexs, aux


def refuses(aux, bits, K):
    """No colouring of the auxiliaries respects the induced lists."""
    lists, P = [], []
    for q, a, b, pi, pj in aux:
        ca = pi ^ ((bits >> a) & 1)
        cb = pj ^ ((bits >> b) & 1)
        used = {ca, cb}
        lists.append([c for c in range(K) if c not in used])
        P.append(q)
    n = len(P)
    var, cnt = {}, 0
    for i in range(n):
        for c in lists[i]:
            cnt += 1
            var[(i, c)] = cnt
    cls = [[var[(i, c)] for c in lists[i]] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = (P[i][0] - P[j][0]) ** 2 + (P[i][1] - P[j][1]) ** 2
            if abs(d - 1.0) < 1e-7:
                for c in set(lists[i]) & set(lists[j]):
                    cls.append([-var[(i, c)], -var[(j, c)]])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return not s.solve()


t0, best = time.time(), None
# Coincident hexagons are the SAME six points, so they cannot carry
# independent orientation bits; only strictly increasing angle choices are
# genuine configurations.
for combo in itertools.combinations(range(1, len(ANGLES)), H - 1):
    phis = (0.0,) + tuple(ANGLES[i] for i in combo)
    _hexs, aux = build(phis)
    got = sum(1 for bits in range(1 << H) if refuses(aux, bits, K))
    if best is None or got > best[0]:
        best = (got, phis, len(aux))
        print(f"  angles {[round(math.degrees(p),4) for p in phis]}: refuses "
              f"{got}/{1 << H} orientations, {len(aux)} auxiliaries  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if got == (1 << H):
            print(f"  *** ALL orientations refused at k={K} -- pressure 3 ***",
                  flush=True)
            break
print(f"H={H} k={K}: best {best[0]}/{1 << H} orientations refused, angles "
      f"{[round(math.degrees(p),4) for p in best[1]]}  [{time.time()-t0:.0f}s]",
      flush=True)
