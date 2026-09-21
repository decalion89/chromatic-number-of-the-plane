"""Strip the frontier with a k-core, which is exactly what it is for.

Density climbs with the cap -- 5.65 at 20000 points, 6.14 at 45000, 6.44 at
90000 -- but slowly, because a BFS ball always carries a constant fraction of
its points on the frontier, and those have most of their neighbours missing.
The limit |U| is approached logarithmically and 90000 points is already the
largest graph built here.

The k-core removes precisely those points: delete every vertex of degree below
k, repeat until none remain.  What survives has minimum degree at least k, and
it is a SUBGRAPH, so it inherits the only property that matters -- if the core
needs six colours then so does the graph, and so does the plane.

So build once, then peel.  Each core is denser than the ball it came from and
smaller, which makes it both a better candidate and a cheaper SAT call.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import deque
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from hn.fast import IntBasis
from pysat.solvers import Solver

k5 = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())
W = _rot60(K)
ROT = {D: rotation_joining(Fr(D), K) for D in (2, 3, 4, 7)}


def units(rots, emax):
    out, seen = [], set()
    for a in range(6):
        for es in itertools.product(*[range(-emax, emax + 1) for _ in rots]):
            p = ONE
            for _ in range(a):
                p = W(p)
            for D, e in zip(rots, es):
                r = ROT[D]
                use = r if e >= 0 else Rotation(r.cos, -r.sin)
                for _ in range(abs(e)):
                    p = use(p)
            if p not in seen:
                seen.add(p)
                out.append(p)
    return out


def walk(U, R, cap):
    basis = IntBasis.covering(U + [Point(K.zero(), K.zero())])
    rows = basis.rows(U)
    dim, D = basis.dim, basis.D
    steps = np.concatenate([rows, -rows], axis=0)
    lim = R * R * D * D

    def norms(V):
        return (basis._field_square(V[:, :dim])
                + basis._field_square(V[:, dim:])) @ basis.sqrts

    z = np.zeros(2 * dim, dtype=np.int64)
    idx = {z.tobytes(): 0}
    out = [z]
    q = deque([0])
    while q and len(out) < cap:
        i = q.popleft()
        cand = steps + out[i]
        for w in cand[norms(cand) <= lim]:
            b = w.tobytes()
            if b not in idx:
                if len(out) >= cap:
                    continue
                idx[b] = len(out)
                out.append(w)
                q.append(idx[b])
    arr = np.array(out, dtype=np.int64)
    E = set()
    for u in steps:
        for i, w in enumerate(arr + u):
            j = idx.get(w.tobytes())
            if j is not None and i != j:
                E.add((min(i, j), max(i, j)))
    return len(out), sorted(E)


def kcore(n, E, kk):
    adj = [set() for _ in range(n)]
    for a, b in E:
        adj[a].add(b)
        adj[b].add(a)
    alive = [True] * n
    q = deque(v for v in range(n) if len(adj[v]) < kk)
    while q:
        v = q.popleft()
        if not alive[v]:
            continue
        alive[v] = False
        for u in adj[v]:
            if alive[u]:
                adj[u].discard(v)
                if len(adj[u]) < kk:
                    q.append(u)
    keep = [v for v in range(n) if alive[v]]
    ren = {v: i for i, v in enumerate(keep)}
    EE = [(ren[a], ren[b]) for a, b in E if alive[a] and alive[b]]
    return len(keep), EE


def solve(n, E):
    cls = [[1 + v * k5 + c for c in range(k5)] for v in range(n)]
    for a, b in E:
        for c in range(k5):
            cls.append([-(1 + a * k5 + c), -(1 + b * k5 + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


for rots, emax, R, cap in (((4,), 2, 2.2, 45000), ((4,), 1, 2.2, 45000)):
    U = units(rots, emax)
    n, E = walk(U, R, cap)
    print(f"\nrotations {rots} e<={emax}, {len(U)} units, R={R}: {n} pts, "
          f"{len(E)} edges, {len(E)/n:.2f} per vertex"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for kk in (8, 12, 16, 20, 24, 28):
        m, EE = kcore(n, E, kk)
        if m < 50:
            print(f"   {kk}-core: empty  [{time.time()-t0:.0f}s]", flush=True)
            break
        ok = solve(m, EE)
        print(f"   {kk}-core: {m} pts, {len(EE)} edges, {len(EE)/m:.2f} per "
              f"vertex -> "
              f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if not ok:
            print("*** a subgraph of the plane needs six colours ***",
                  flush=True)
            sys.exit(0)
print("\nDONE", flush=True)
