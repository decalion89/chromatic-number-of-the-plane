"""Colouring-guided growth with a numpy candidate pool (grow_kw.py, about 10x less memory).

Same loop as grow_kw.py:
1. repair the colouring with tabu search, and fall back to kissat;
2. insert the R richest rainbow pool points, or conflicting tight pairs when too
   few are rainbow.

The difference is the pool. grow_kw.py keeps every candidate point as a Python
dict entry holding a set of neighbours; with 666 directions and 30 000 points
that is 10+ GB. Here candidates live in sorted numpy arrays of 64-bit hashed
coordinates, with the base vertex, unit index and neighbour count per candidate.
Neighbour colours are recomputed on demand, and only for candidates with at
least 4 neighbours.

Every unit-distance pair is found, not only those along U: a KD-tree finds pairs at
float distance 1, one pair per unknown direction is checked exactly, and each new
direction joins U with its negative and conjugates.  (Without this the colourings
exploited unseen edges: 6 788 monochromatic true edges at 19 851 points in F8.)

MODE plain | apart (c(A) = c(B) imposed) | same (c(A) != c(B) imposed).
Any UNSAT is a claim until verify6 / verify_gadget re-check it exactly.

usage: grow_lean.py <in.json> <out.json> [R] [near_weight]
env: MODE, KISSAT, KTIME, LSIT, IT2, TABU2
"""
import sys, time, json, os, subprocess
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
t0 = time.time(); K = 5
IN, OUT = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 40
WN = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
MODE = os.environ.get("MODE", "plain")
KB = os.environ.get("KISSAT", "kissat"); KTIME = int(os.environ.get("KTIME", "7200"))
_here = os.path.dirname(os.path.abspath(__file__))
TABU2 = os.environ.get("TABU2", os.path.join(_here, "tabu2"))
if not os.path.exists(TABU2):
    subprocess.run(["gcc", "-O2", "-o", TABU2, os.path.join(_here, "tabu2.c")], check=True)
LSIT = int(os.environ.get("LSIT", "3000000")); IT2 = int(os.environ.get("IT2", "0"))

d = json.load(open(IN))
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]
U = [mk(xy) for xy in d["units"]]
A, Bi = d.get("A"), d.get("B")
nU = len(U)
Ux = np.array([u.fx for u in U]); Uy = np.array([u.fy for u in U])
Vx = [p.fx for p in V]; Vy = [p.fy for p in V]

M1 = np.int64(0x9E3779B97F4A7C15 - (1 << 64)); M2 = np.int64(0x632BE59BD9B4E019)


def hkey(x, y):
    """64-bit hash of coordinates rounded to 1e-8 (vectorised)."""
    kx = np.rint(np.asarray(x) * 1e8).astype(np.int64)
    ky = np.rint(np.asarray(y) * 1e8).astype(np.int64)
    with np.errstate(over="ignore"):
        return kx * M1 + ky * M2


# vertex keys: dict for exact bookkeeping, sorted array for vectorised lookups
vkeys = hkey(np.array(Vx), np.array(Vy))
assert len(np.unique(vkeys)) == len(vkeys), "duplicate points or a hash collision among vertices"
vk_sorted_idx = np.argsort(vkeys); vk_sorted = vkeys[vk_sorted_idx]


def vlookup(keys):
    """vertex index for each key, or -1."""
    pos = np.searchsorted(vk_sorted, keys)
    pos[pos >= len(vk_sorted)] = 0
    hit = vk_sorted[pos] == keys
    return np.where(hit, vk_sorted_idx[pos], -1)


def add_vertex_keys(newkeys, newidx):
    global vk_sorted, vk_sorted_idx
    allk = np.concatenate([vk_sorted, newkeys]); alli = np.concatenate([vk_sorted_idx, newidx])
    o = np.argsort(allk, kind="stable"); vk_sorted = allk[o]; vk_sorted_idx = alli[o]


# edges
adj = [set() for _ in V]; EA, EB = [], []
for i in range(len(V)):
    nb = vlookup(hkey(Vx[i] + Ux, Vy[i] + Uy))
    for j in nb[nb > i]:
        j = int(j)
        if j not in adj[i]:
            adj[i].add(j); adj[j].add(i); EA.append(i); EB.append(j)

# candidate pool: sorted keys, base vertex, unit index, neighbour count
ck = np.zeros(0, dtype=np.int64); cbase = np.zeros(0, dtype=np.int32); cunit = np.zeros(0, dtype=np.int16)
ccnt = np.zeros(0, dtype=np.int16)


def feed(idxs, units=None):
    """add the neighbour candidates of the vertices idxs (along the given unit indices, default all) to the pool"""
    global ck, cbase, cunit, ccnt
    if len(idxs) == 0:
        return
    units = np.arange(nU) if units is None else np.asarray(units)
    bx = np.array([Vx[i] for i in idxs]); by = np.array([Vy[i] for i in idxs])
    keys = hkey(bx[:, None] + Ux[None, units], by[:, None] + Uy[None, units]).ravel()
    base = np.repeat(np.array(idxs, dtype=np.int32), len(units)); unit = np.tile(units.astype(np.int16), len(idxs))
    isv = vlookup(keys) >= 0
    keys, base, unit = keys[~isv], base[~isv], unit[~isv]
    uk, first, cnt = np.unique(keys, return_index=True, return_counts=True)
    pos = np.searchsorted(ck, uk); pos2 = np.minimum(pos, max(len(ck) - 1, 0))
    found = (len(ck) > 0) & (ck[pos2] == uk) if len(ck) else np.zeros(len(uk), dtype=bool)
    if found.any():
        np.add.at(ccnt, pos2[found], cnt[found].astype(np.int16))
    nk = uk[~found]
    if len(nk):
        # nk is sorted and disjoint from ck: a linear-time merge
        at = np.searchsorted(ck, nk)
        ck = np.insert(ck, at, nk); cbase = np.insert(cbase, at, base[first[~found]])
        cunit = np.insert(cunit, at, unit[first[~found]]); ccnt = np.insert(ccnt, at, cnt[~found].astype(np.int16))


def remove_candidates(keys):
    global ck, cbase, cunit, ccnt
    if len(ck) == 0: return
    pos = np.searchsorted(ck, keys)
    valid = pos < len(ck)
    pos_v, keys_v = pos[valid], keys[valid]
    drop = pos_v[ck[pos_v] == keys_v]
    keep = np.ones(len(ck), dtype=bool); keep[drop] = False
    ck, cbase, cunit, ccnt = ck[keep], cbase[keep], cunit[keep], ccnt[keep]


from scipy.spatial import cKDTree
ukey = lambda q: (tuple(q.x.c), tuple(q.y.c))
UK = set(ukey(u) for u in U)
ONE = F.rational(1)


def discover(only=None):
    """Find unit-distance pairs whose direction is not in U (float KD-tree, then an exact check of one
    pair per direction), add those directions with their negatives and conjugates to U, add every edge
    along them, and extend the pool.  Returns the number of new directions."""
    global Ux, Uy, nU
    pts = np.column_stack([np.array(Vx), np.array(Vy)]); tree = cKDTree(pts)
    if only is None:
        pr = tree.query_pairs(1 + 1e-7, output_type="ndarray")
    else:
        lst = tree.query_ball_point(pts[only], 1 + 1e-7)
        pr = np.array([(i, j) for i, js in zip(only, lst) for j in js if j != i], dtype=np.int64).reshape(-1, 2)
    if len(pr) == 0: return 0
    dx = pts[pr[:, 1], 0] - pts[pr[:, 0], 0]; dy = pts[pr[:, 1], 1] - pts[pr[:, 0], 1]
    keep = np.abs(np.hypot(dx, dy) - 1) < 1e-7
    pr, dx, dy = pr[keep], dx[keep], dy[keep]
    dk = hkey(dx, dy); known = np.isin(dk, np.array(sorted(set(hkey(Ux, Uy).tolist())), dtype=np.int64))
    pr, dk = pr[~known], dk[~known]
    fresh = []
    seenk = set()
    for (i, j), k in zip(pr, dk):
        if int(k) in seenk: continue
        seenk.add(int(k)); w = V[int(j)] - V[int(i)]
        if w.norm2() != ONE: continue                           # a float coincidence, not a unit
        for z in (w, -w, Point(w.x, -w.y), Point(-w.x, w.y)):
            if ukey(z) not in UK: UK.add(ukey(z)); U.append(z); fresh.append(len(U) - 1)
    if not fresh: return 0
    Ux = np.array([u.fx for u in U]); Uy = np.array([u.fy for u in U]); nU = len(U)
    for v in range(len(V)):                                    # every edge along the new directions
        nb = vlookup(hkey(Vx[v] + Ux[fresh], Vy[v] + Uy[fresh]))
        for j in nb[nb > v]:
            j = int(j)
            if j not in adj[v]:
                adj[v].add(j); adj[j].add(v); EA.append(v); EB.append(j)
    for _s in range(0, len(V), 4000):
        feed(list(range(_s, min(len(V), _s + 4000))), fresh)
    return len(fresh)


for _s in range(0, len(V), 4000):
    feed(list(range(_s, min(len(V), _s + 4000))))
_m0 = len(EA); _nd = discover()
print(f"discovery: {_nd} new unit directions, {len(EA) - _m0} edges along them", flush=True)

ser = lambda q: [[[t.numerator, t.denominator] for t in q.x.c], [[t.numerator, t.denominator] for t in q.y.c]]
col = list(d["colouring"]) if d.get("colouring") and len(d["colouring"]) == len(V) else []


def save(tag):
    json.dump({"field_generators": list(F.gens), "A": A, "B": Bi, "status": tag, "mode": MODE,
               "colouring": col[:len(V)] if len(col) == len(V) else None, "units": [ser(u) for u in U],
               "two_edges": [], "points": [ser(q) for q in V]}, open(OUT, "w"))


tri = next(((a, b, c) for a, b in zip(EA, EB) for c in sorted(adj[a] & adj[b])), None)


def proper(c):
    if any(c[a] == c[b] for a, b in zip(EA, EB)): return False
    if MODE == "apart" and c[A] != c[Bi]: return False
    if MODE == "same" and c[A] == c[Bi]: return False
    return True


def tabucol(n, init, seed):
    ea, eb = EA, EB
    if MODE == "apart":
        ea = [A if x == Bi else x for x in EA]; eb = [A if x == Bi else x for x in EB]
    elif MODE == "same":
        ea = EA + [A]; eb = EB + [Bi]
    inp = (f"{n} {len(ea)} 0 {K} {LSIT} {IT2} {seed}\n" + "\n".join(f"{a} {b}" for a, b in zip(ea, eb)) + "\n"
           + "\n".join(map(str, init)) + "\n")
    out = subprocess.run([TABU2], input=inp, capture_output=True, text=True).stdout.split("\n")
    ok = out[0].startswith("OK"); f = int(out[0].split()[1]); c = [int(x) for x in out[1:1 + n]]
    if MODE == "apart": c[Bi] = c[A]
    return ok, f, c


def kissat_solve(n):
    cnf = OUT + ".cnf"
    X = lambda v, c: 1 + v * K + c
    with open(cnf, "w") as fh:
        cls = [[X(v, c) for c in range(K)] for v in range(n)]
        cls += [[-X(a, c), -X(b, c)] for a, b in zip(EA, EB) for c in range(K)]
        if MODE == "apart": cls += [[-X(A, c), X(Bi, c)] for c in range(K)] + [[X(A, c), -X(Bi, c)] for c in range(K)]
        elif MODE == "same": cls += [[-X(A, c), -X(Bi, c)] for c in range(K)]
        if tri: cls += [[X(v, k)] for k, v in enumerate(tri)]
        fh.write(f"p cnf {n * K} {len(cls)}\n"); fh.write("".join(" ".join(map(str, c)) + " 0\n" for c in cls))
    kout = subprocess.run([KB, f"--time={KTIME}", cnf], capture_output=True, text=True).stdout
    status = next((l for l in kout.split("\n") if l.startswith("s ")), "s UNKNOWN")
    if "UNSATISFIABLE" in status: return False, None, status
    if "SATISFIABLE" in status:
        pos = {int(x) for l in kout.split("\n") if l.startswith("v") for x in l.split()[1:] if int(x) > 0}
        return True, [next((c for c in range(K) if X(v, c) in pos), 0) for v in range(n)], status
    return None, None, status


def neighbour_colours(sel):
    """for candidates sel (indices into the pool): bitmask of colours among graph neighbours, and a
    count of neighbours (in chunks, to bound memory)"""
    colarr = np.array(col + [0], dtype=np.int64)
    Vxa = np.array(Vx); Vya = np.array(Vy)
    masks, degs = [], []
    for s0 in range(0, len(sel), 8000):
        ss = sel[s0:s0 + 8000]
        px = Vxa[cbase[ss]] + Ux[cunit[ss]]; py = Vya[cbase[ss]] + Uy[cunit[ss]]
        nb = vlookup(hkey(px[:, None] + Ux[None, :], py[:, None] + Uy[None, :]))      # m x nU
        bits = np.where(nb >= 0, np.left_shift(1, colarr[nb]), 0)
        masks.append(np.bitwise_or.reduce(bits, axis=1)); degs.append((nb >= 0).sum(axis=1))
    return np.concatenate(masks), np.concatenate(degs), None, None


if A is not None:
    ax, ay, bx_, by_ = Vx[A], Vy[A], Vx[Bi], Vy[Bi]
    near = lambda x, y: np.minimum(np.hypot(x - ax, y - ay), np.hypot(x - bx_, y - by_))
else:
    cx, cy = float(np.mean(Vx)), float(np.mean(Vy))
    near = lambda x, y: np.hypot(x - cx, y - cy)

print(f"{IN}: n={len(V)} m={len(EA)}, {nU} directions, pool {len(ck)}; mode {MODE}"
      + (f", pair {A},{Bi} at d^2 = {V[A].dist2(V[Bi])}" if A is not None else "") + f"   [{time.time()-t0:.0f}s]", flush=True)
lsfail = 0
full = (1 << K) - 1
for it in range(1, 10 ** 7):
    n = len(V)
    init = col + [-1] * (n - len(col))
    ok, info, lcol = tabucol(n, init, it)
    if ok and proper(lcol):
        r = True; col = lcol; how = f"LS {info}"
    else:
        lsfail += 1
        save("pre_cdcl")
        r, kcol, status = kissat_solve(n)
        how = f"KISSAT[{status}]"
        if r: col = kcol
        if r and not proper(col): r = None
    if r is not True:
        save(("UNSAT_" + MODE) if r is False else "hard")
        print(f"  iter {it}: n={n} m={len(EA)}: {'UNSAT (' + MODE + ') -- a claim until verified' if r is False else 'unknown'}   [{time.time()-t0:.0f}s]", flush=True)
        break
    sel = np.nonzero(ccnt >= 4)[0]
    mask, deg, px, py = neighbour_colours(sel) if len(sel) else (np.zeros(0, dtype=np.int64),) * 4
    rb = sel[(mask == full) & (deg >= 5)] if len(sel) else sel
    score_rb = 0.25 * ccnt[rb] - WN * near(np.array([Vx[b] for b in cbase[rb]]) + Ux[cunit[rb]],
                                         np.array([Vy[b] for b in cbase[rb]]) + Uy[cunit[rb]]) if len(rb) else np.zeros(0)
    if it % 10 == 1 or len(rb) < R or how.startswith("KISSAT"):
        print(f"  iter {it}: n={n} m={len(EA)} {how}; {len(rb)} rainbow, pool {len(ck)}   [{time.time()-t0:.0f}s]", flush=True)
    if it % 10 == 0: save("checkpoint")
    pick = list(rb[np.argsort(-score_rb)][:R])
    if len(pick) < R and len(sel):
        # tight: exactly one colour missing among >= 4 neighbours; pair adjacent tight candidates
        # that miss the same colour (the colouring cannot extend to both)
        tmask = mask ^ full
        tight_sel = sel[(np.bitwise_and(tmask, tmask - 1) == 0) & (tmask != 0)]
        tmiss = tmask[(np.bitwise_and(tmask, tmask - 1) == 0) & (tmask != 0)]
        if len(tight_sel):
            tx = np.array([Vx[b] for b in cbase[tight_sel]]) + Ux[cunit[tight_sel]]
            ty = np.array([Vy[b] for b in cbase[tight_sel]]) + Uy[cunit[tight_sel]]
            tk = hkey(tx, ty); o = np.argsort(tk); tks = tk[o]
            nbk = hkey(tx[:, None] + Ux[None, :], ty[:, None] + Uy[None, :])
            posn = np.searchsorted(tks, nbk); posn[posn >= len(tks)] = 0
            hit = tks[posn] == nbk
            ii, jj = np.nonzero(hit)
            jj = o[posn[ii, jj]]
            same = tmiss[ii] == tmiss[jj]
            ii, jj = ii[same], jj[same]
            have = set(int(x) for x in pick)
            for a_, b_ in zip(ii, jj):
                if len(pick) >= R: break
                for x_ in (tight_sel[a_], tight_sel[b_]):
                    if int(x_) not in have: pick.append(x_); have.add(int(x_))
            if it % 10 == 1: print(f"    {len(tight_sel)} tight candidates, {len(ii)} conflicting tight pairs", flush=True)
    if not pick:
        o = np.argsort(-ccnt)[:R]; pick = list(o)
        print(f"    no rainbow, no tight pair: densifying with the {len(pick)} richest candidates", flush=True)
    if not pick:
        save("pool_exhausted"); print("  pool exhausted", flush=True); break
    newidx = []
    for c_ in pick:
        b, k = int(cbase[c_]), int(cunit[c_])
        p = V[b] + U[k]; v = len(V)
        V.append(p); Vx.append(p.fx); Vy.append(p.fy); adj.append(set()); newidx.append(v)
    newkeys = hkey(np.array(Vx[n:]), np.array(Vy[n:]))
    add_vertex_keys(newkeys, np.array(newidx, dtype=np.int64))
    remove_candidates(np.sort(newkeys))
    for v in newidx:
        nb = vlookup(hkey(Vx[v] + Ux, Vy[v] + Uy))
        for j in nb[nb >= 0]:
            j = int(j)
            if j != v and j not in adj[v]:
                adj[v].add(j); adj[j].add(v); EA.append(v); EB.append(j)
    _m0 = len(EA); _nd = discover(newidx)
    if _nd: print(f"    discovery: {_nd} new unit directions ({nU} now), {len(EA) - _m0} edges along them", flush=True)
    feed(newidx)
