"""minimize.py -- shrink a non-4-colourable point set: core iteration, then deletion to vertex-criticality.

usage: python3 minimize.py PTS.npy TAG [tlimit_per_call] [tri_a,tri_b,tri_c]
One selector per vertex (ALO clause guarded by -s_v); the colours of a triangle are fixed (valid symmetry
breaking); solve under assumptions {s_v}; get_core -> repeat; then try deleting each vertex (keep the deletion
if still UNSAT, and jump to the new core). Vertices of degree < 4 in the current subgraph are dropped
(4-core), which preserves non-colourability. Every improvement is saved to min_TAG_pts.npy.
"""
import sys, time, json
import numpy as np
from flat import *
from satutil import solve_limited
from pysat.solvers import Solver

D, CD, U = directions()
P = np.load(sys.argv[1])
tag = sys.argv[2]
TL = float(sys.argv[3]) if len(sys.argv) > 3 else 60
logf = open(f"min_{tag}.log", "a")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); logf.write(s + "\n"); logf.flush()


n = len(P)
E, J = build_edges(P, U)
adj = [set() for _ in range(n)]
for a, b in E:
    adj[a].add(int(b)); adj[b].add(int(a))
if len(sys.argv) > 4:
    tri = tuple(int(x) for x in sys.argv[4].split(","))
else:
    S0 = PointSet(P)
    t0p = ([0] * 12, [7] + [0] * 11, list(D[14]))
    idx = [int(S0.index([t])[0]) for t in t0p]
    tri = tuple(idx) if min(idx) >= 0 else find_triangle(n, E.tolist())
assert tri[1] in adj[tri[0]] and tri[2] in adj[tri[0]] and tri[2] in adj[tri[1]]
triset = set(tri)


def var(v, c):
    return 1 + 4 * v + c


def sel(v):
    return 1 + 4 * n + v


s = Solver(name="cadical195")
for v in range(n):
    s.add_clause([-sel(v)] + [var(v, c) for c in range(4)])
for a, b in E:
    for c in range(4):
        s.add_clause([-var(int(a), c), -var(int(b), c)])
for i, v in enumerate(tri):
    s.add_clause([var(v, i)])


def core4(Sset):
    """4-core of the induced subgraph on Sset (triangle kept)"""
    Sset = set(Sset) | triset
    changed = True
    while changed:
        changed = False
        for v in list(Sset):
            if v in triset:
                continue
            if len(adj[v] & Sset) < 4:
                Sset.discard(v); changed = True
    return Sset


def test(Sset, tl):
    r = solve_limited(s, [sel(v) for v in Sset], tlimit=tl)
    if r is False:
        core = set(-l - 0 for l in []) or set(l - 1 - 4 * n for l in s.get_core())
        return False, core4(core)
    return r, None


def save(Sset, stage):
    L = sorted(Sset)
    np.save(f"min_{tag}_pts.npy", P[L])
    json.dump({"n": len(L), "stage": stage, "triangle_in_saved": [L.index(v) for v in tri]},
              open(f"min_{tag}_state.json", "w"))


T0 = time.time()
cur = core4(range(n))
log(f"=== minimize {sys.argv[1]} tag={tag}: {n} vertices, {len(E)} edges; 4-core {len(cur)}; triangle {tri}; tlimit/call {TL}s")
t = time.time()
r, core = test(cur, None)
if r is not False:
    log(f"full graph is not UNSAT (result {r})"); sys.exit(1)
log(f"stage 0: UNSAT in {time.time() - t:.1f}s; core (4-cored) {len(core)}")
cur = core
save(cur, "core0")
# stage 1: core iteration
it = 0
while True:
    it += 1
    t = time.time()
    r, core = test(cur, None)
    assert r is False
    log(f"stage 1 iter {it}: UNSAT in {time.time() - t:.1f}s; core {len(core)}  [total {time.time() - T0:.0f}s]")
    if len(core) >= len(cur):
        break
    cur = core
    save(cur, f"core{it}")
log(f"after core iteration: {len(cur)} vertices")
# stage 2: deletion-based minimisation
tried = set()
passes = 0
while True:
    passes += 1
    progress = False
    order = sorted((v for v in cur if v not in triset and v not in tried), key=lambda v: (len(adj[v] & cur), v))
    for v in order:
        if v not in cur:
            continue
        cand = core4(cur - {v})
        t = time.time()
        r, core = test(cand, TL)
        dt = time.time() - t
        if r is False:
            cur = core
            progress = True
            save(cur, "deletion")
            log(f"  del {v}: UNSAT in {dt:.1f}s -> {len(cur)}  [total {time.time() - T0:.0f}s]")
        else:
            tried.add(v)
            if r is None:
                log(f"  del {v}: time-out ({dt:.0f}s), kept")
    log(f"stage 2 pass {passes}: {len(cur)} vertices  [total {time.time() - T0:.0f}s]")
    if not progress:
        break
    tried = set(v for v in tried if v in cur)
    # after progress, previously-kept vertices might now be removable only if the set changed; recheck all once
    tried = set()
save(cur, "critical")
log(f"final: {len(cur)} vertices (vertex-critical w.r.t. the tests done; time-outs counted as kept)")
