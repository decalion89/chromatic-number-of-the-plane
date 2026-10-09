"""growforce4.py -- colouring-guided growth at K colours aimed at a forced pair (env KCOL, default 4).

usage: FORCE="a,b,c,e" python3 growforce4.py TAG [maxn] [kissat_tl] [D] [radius] [per_round]

The plane over Q(sqrt d), d = QD (env, read by kd.py; default 47); points are integer 4-vectors (a, b, c, e) / D
standing for (a + b sqrt d)/D + i (c + e sqrt d)/D, and the directions are all unit vectors with denominator dividing D. FORCE is
a vector m of the module (in units 1/D). The graph always contains 0 and m, and every colouring used must give them
different colours (an extra constraint, not an edge of the plane). Each round K-colours the graph (tabu from the
previous colouring, else kissat), then adds candidates p + u whose neighbours see all K colours, in that colouring or
in MULTI fresh tabu colourings (env, default 3); with none, it recolours from scratch (RECOL) and then grows softly
(at most SOFT rounds that add candidates seeing K - 1 colours). When kissat finds no colouring with c(0) != c(m),
the run checks the graph again without the pair: if it is K-colourable, 0 and m have the same colour in every
K-colouring ('forced'; then a rotation about 0 that moves m by exactly 1 gives a spindle, spin4.py); if not, the graph
is a witness by itself ('plain_unsat'). Neither answer is a proof: certify it separately (DRAT).
Environment: KISSAT (path), SEED (offset for the random seeds), OUTDIR (state and logs; default this directory),
MULTI, RECOL, SOFT. Points are looked up by 64-bit hash keys, every match checked exactly. State saved every round
(OUTDIR/growforce4_TAG_state.json); a run with the same TAG resumes from it."""
import sys, os, json, time, subprocess, itertools
from fractions import Fraction as Fr
from math import gcd
import numpy as np
from kd import mul, conj, unit_from, norm1, R
from tabu import tabucol

HERE = os.path.dirname(os.path.abspath(__file__))
SC = os.path.dirname(HERE)
KISSAT = os.environ.get("KISSAT", f"{SC}/kissat/build/kissat")
OUT = os.environ.get("OUTDIR", HERE)
SEED = int(os.environ.get("SEED", "0"))
tag = sys.argv[1]
MAXN = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
TL = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
NGEN = int(sys.argv[4]) if len(sys.argv) > 4 else 2
RAD = int(sys.argv[5]) if len(sys.argv) > 5 else 2
PER = int(sys.argv[6]) if len(sys.argv) > 6 else 200
logf = open(f"{OUT}/growforce4_{tag}.log", "a")
K = int(os.environ.get("KCOL", "4"))
FORCE = np.array([int(x) for x in os.environ["FORCE"].split(",")], dtype=np.int64)


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); logf.write(s + "\n"); logf.flush()


# ---- directions: every unit vector (a + b r + c i + d i r)/D with integers a, b, c, d (D = NGEN here)
from math import isqrt


def units_den(Dd):
    out = []
    D2 = Dd * Dd
    B = isqrt(D2 // R)
    for b in range(-B, B + 1):
        for d in range(-B, B + 1):
            rest = D2 - R * (b * b + d * d)
            if rest < 0:
                continue
            A = isqrt(rest)
            for a_ in range(-A, A + 1):
                c2 = rest - a_ * a_
                c = isqrt(c2)
                if c * c != c2:
                    continue
                for cc in ({c, -c} if c else {0}):
                    if a_ * b + cc * d == 0:
                        out.append((a_, b, cc, d))
    return sorted(set(out))


D = NGEN
from units_fast import units_fast
U = np.array(units_fast(R, D), dtype=np.int64)
if D <= 300:
    assert sorted(map(tuple, U.tolist())) == units_den(D)
for u in U.tolist():
    t = tuple(Fr(v, D) for v in u)
    assert norm1(t)
T = U
log(f"=== grow3 {tag}: {len(U)} unit vectors with denominator dividing D = {D}")

HW = np.array([1000003, 998244353, 19260817, 1610612741], dtype=np.int64)


OFF = 1 << 15                       # exact keys: four 16-bit fields (every coordinate must be below 2^15 - 1)
SENTINEL = np.uint64(0xFFFFFFFFFFFFFFFF)


_HC = np.array([0x9E3779B97F4A7C15, 0xC2B2AE3D27D4EB4F, 0x165667B19E3779F9, 0xD6E8FEB86659FD93], dtype=np.uint64)


def keys(A):
    """64-bit hash keys of integer 4-vectors of any size (PointSet.index checks every match exactly, and the
    constructor refuses two points with the same key, so a collision can only make it stop, never err)"""
    A = np.asarray(A, dtype=np.int64).reshape(-1, 4)
    with np.errstate(over="ignore"):
        k = (A.astype(np.uint64) * _HC).sum(axis=1, dtype=np.uint64)
        k ^= k >> np.uint64(29); k *= np.uint64(0xBF58476D1CE4E5B9); k ^= k >> np.uint64(32)
    k[k == SENTINEL] = SENTINEL - np.uint64(1)
    return k


class PointSet:
    def __init__(self, P):
        self.P = np.asarray(P, dtype=np.int64).reshape(-1, 4)
        k = keys(self.P)
        if np.any(k == SENTINEL):
            raise ValueError("a coordinate is too large for the exact keys")
        self.order = np.argsort(k, kind='stable')
        self.sk = k[self.order]
        if len(k) > 1 and np.any(self.sk[1:] == self.sk[:-1]):
            raise ValueError("duplicate keys")

    def index(self, Q):
        Q = np.asarray(Q, dtype=np.int64).reshape(-1, 4)
        k = keys(Q)
        pos = np.searchsorted(self.sk, k)
        pos[pos >= len(self.sk)] = 0
        hit = self.sk[pos] == k
        idx = np.where(hit, self.order[pos], -1)
        ok = idx >= 0
        if ok.any():
            same = np.all(self.P[idx[ok]] == Q[ok], axis=1)
            tmp = idx[ok]; tmp[~same] = -1; idx[ok] = tmp
        return idx


def build_edges(P):
    S = PointSet(P)
    E = []
    for u in U:
        idx = S.index(P + u)
        a = np.nonzero(idx >= 0)[0]
        b = idx[a]
        keep = a < b
        E.append(np.stack([a[keep], b[keep]], axis=1))
    return np.concatenate(E)


def kissat3(n, E, tl):
    path = f"/dev/shm/gf4_{tag}_{os.getpid()}.cnf"
    cl = [[K * v + c + 1 for c in range(K)] for v in range(n)]
    cl += [[-(K * a + c + 1), -(K * b + c + 1)] for a, b in E.tolist() for c in range(K)]
    a0, b0 = E[0]
    cl += [[K * int(a0) + 1], [K * int(b0) + 2]]          # one edge (the forcing pair, E[0]) fixed to colours 0, 1
    with open(path, "w") as f:
        f.write(f"p cnf {K * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    out = subprocess.run([KISSAT, f"--time={tl}", path], capture_output=True, text=True).stdout
    st = next((l for l in out.splitlines() if l.startswith("s ")), "s UNKNOWN")
    if "UNSATISFIABLE" in st:
        return False, None
    os.unlink(path)
    if "SATISFIABLE" in st:
        lits = set(int(x) for l in out.splitlines() if l.startswith("v ") for x in l.split()[1:] if int(x) > 0)
        return True, np.array([next(c for c in range(K) if K * v + c + 1 in lits) for v in range(n)])
    return None, None


# ---- seed: sums of at most RAD unit vectors
st_file = f"{OUT}/growforce4_{tag}_state.json"
if os.path.exists(st_file):
    stt = json.load(open(st_file))
    P = np.array(stt["P"], dtype=np.int64)
    col = np.array(stt["col"]) if stt.get("col") else None
    log(f"resumed: {len(P)} points")
else:
    P = np.zeros((1, 4), dtype=np.int64)
    for _ in range(RAD):
        P = np.unique(np.concatenate([P] + [P + u for u in U]), axis=0)
    P = np.unique(np.concatenate([P, P + FORCE]), axis=0)          # balls around 0 and around m
    col = None
S0 = PointSet(P)
i0, im = (int(x) for x in S0.index(np.stack([np.zeros(4, dtype=np.int64), FORCE])))
assert i0 >= 0 and im >= 0
log(f"forcing target m = {FORCE.tolist()} / {D}; vertices {i0} (origin) and {im} (m)")
T0 = time.time()
rnd = 0
soft = 0
RECOL = int(os.environ.get("RECOL", "4"))
SOFT = int(os.environ.get("SOFT", "6"))
while True:
    rnd += 1
    n = len(P)
    E = np.concatenate([np.array([[i0, im]], dtype=np.int64), build_edges(P)])   # E[0] is the forcing pair
    t1 = time.time()
    init = -np.ones(n, dtype=np.int64)
    if col is not None:
        init[:len(col)] = col
    c2 = None
    for s_ in range(6):
        c2 = tabucol(n, E, init=init, k=K, maxiter=2_000_000, seed=100 * rnd + s_ + 1 + 1000003 * SEED)
        if c2 is not None:
            how = f"tabu/{s_ + 1}"; break
    if c2 is None:
        how = "kissat"
        r, c2 = kissat3(n, E, TL)
        if r is False:
            r0, _ = kissat3(n, E[1:], TL)          # the same graph without the forcing pair
            if r0 is not True:
                what = f"has no {K}-colouring by itself" if r0 is False else "is undecided by itself (time-out)"
                log(f"round {rnd}: n={n} e={len(E) - 1}: kissat UNSAT, but the graph {what}: NOT a forced pair  "
                    f"[total {time.time() - T0:.0f}s]")
                np.save(f"{OUT}/growforce4_{tag}_unsat_pts.npy", P)
                json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "m": FORCE.tolist(), "i0": i0, "im": im,
                           "status": "plain_unsat" if r0 is False else "plain_timeout"}, open(st_file, "w"))
                break
            log(f"round {rnd}: n={n} e={len(E) - 1}: kissat UNSAT -- FORCED: 0 and m have the same colour in every "
                f"{K}-colouring, and the graph itself is {K}-colourable  [total {time.time() - T0:.0f}s]")
            np.save(f"{OUT}/growforce4_{tag}_forced_pts.npy", P)
            json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "m": FORCE.tolist(), "i0": i0, "im": im,
                       "status": "forced"}, open(st_file, "w"))
            break
        if r is None:
            log(f"round {rnd}: n={n} e={len(E)}: kissat TIME-OUT  [total {time.time() - T0:.0f}s]")
            json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "col": None, "status": "timeout"}, open(st_file, "w"))
            break
    col = np.asarray(c2)
    assert np.all(col[E[:, 0]] != col[E[:, 1]])
    assert col[i0] != col[im]
    json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "col": col.tolist(), "m": FORCE.tolist(), "i0": i0, "im": im,
               "status": "running"}, open(st_file, "w"))
    if n >= MAXN:
        log(f"round {rnd}: n={n}: {K}-colourable, reached maxn"); break
    S = PointSet(P)
    cand = np.unique(np.concatenate([P + u for u in U]), axis=0)
    cand = cand[S.index(cand) < 0]
    C = PointSet(cand)
    mask = np.zeros(len(cand), dtype=np.int64)
    deg = np.zeros(len(cand), dtype=np.int64)
    for u in U:
        idx = C.index(P + u)                       # candidate p + u is a neighbour of p
        ok = idx >= 0
        np.bitwise_or.at(mask, idx[ok], 1 << col[ok])
        np.add.at(deg, idx[ok], 1)
    blocked = np.nonzero(mask == (1 << K) - 1)[0]
    MULTI = int(os.environ.get("MULTI", "3"))
    extra = []
    for r_ in range(MULTI):                       # more colourings: kill several at once
        c3 = tabucol(n, E, init=None, k=K, maxiter=2_000_000, seed=104729 * rnd + r_ + 1 + 1000003 * SEED)
        if c3 is None:
            continue
        m3 = np.zeros(len(cand), dtype=np.int64)
        for u in U:
            idx = C.index(P + u); ok = idx >= 0
            np.bitwise_or.at(m3, idx[ok], 1 << c3[ok])
        extra.append(np.nonzero(m3 == (1 << K) - 1)[0])
    if extra:
        nb0 = len(blocked)
        blocked = np.unique(np.concatenate([blocked] + extra))
        log(f"round {rnd}: blocked {nb0} by the current colouring, {len(blocked)} by it or {len(extra)} fresh ones")
    if len(blocked) == 0:
        # stuck: recolour from scratch (random tabu seeds, then kissat) and look again; then grow softly
        found = False
        for r_ in range(RECOL):
            c3 = tabucol(n, E, init=None, k=K, maxiter=2_000_000, seed=7919 * rnd + r_ + 1 + 1000003 * SEED)
            if c3 is None:
                continue
            m3 = np.zeros(len(cand), dtype=np.int64)
            for u in U:
                idx = C.index(P + u); ok = idx >= 0
                np.bitwise_or.at(m3, idx[ok], 1 << c3[ok])
            b3 = np.nonzero(m3 == (1 << K) - 1)[0]
            if len(b3):
                col, mask, blocked, found = np.asarray(c3), m3, b3, True
                log(f"round {rnd}: recoloured from scratch (try {r_ + 1}): {len(b3)} blocked")
                break
        if not found:
            soft += 1
            pc = np.array([bin(int(x)).count('1') for x in mask]) if len(mask) else np.zeros(0, dtype=np.int64)
            two = np.nonzero((pc == K - 1) & (deg >= K - 1))[0]
            if len(two) == 0:                      # early on (a sparse seed): candidates seeing at least two colours
                two = np.nonzero((pc >= 2) & (deg >= 2))[0]
            if soft > SOFT or len(two) == 0:
                log(f"round {rnd}: n={n}: no blocked candidate after {RECOL} recolourings and {soft - 1} soft rounds; stop"); break
            blocked = two
            log(f"round {rnd}: soft round {soft}: adding candidates whose neighbours see two colours")
    pick = blocked[np.argsort(-deg[blocked], kind='stable')[:PER]]
    log(f"round {rnd}: n={n} e={len(E)} {K}-col by {how} in {time.time() - t1:.1f}s; cand {len(cand)}, blocked "
        f"{len(blocked)}; add {len(pick)} (deg {deg[pick].min()}-{deg[pick].max()})  [total {time.time() - T0:.0f}s]")
    P = np.concatenate([P, cand[pick]])
    col = np.concatenate([col, -np.ones(len(pick), dtype=np.int64)])
