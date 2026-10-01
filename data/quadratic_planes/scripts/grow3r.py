"""grow3.py -- colouring-guided growth at 3 colours for unit-distance graphs over Q(sqrt47).

usage: python3 grow3d.py TAG [maxn] [kissat_tl] [D] [radius] [per_round]
Points are elements of K = Q(sqrt47)(i) written as integer 4-vectors (a, b, c, d) / D, meaning a + b r + c i + d i r.
U: unit vectors t = s / conj(s) for small s, closed under i (90 degrees), conj, negation and r -> -r.
Round: 3-colour the graph (tabucol from the previous colouring, else kissat); candidates p + u not yet in the graph;
blocked = neighbours see all 3 colours; add the blocked ones with most neighbours. Stop at UNSAT (then the graph needs
4 colours: certify with DRAT separately), maxn, or a kissat time-out. State saved every round."""
import sys, os, json, time, subprocess, itertools
from fractions import Fraction as Fr
from math import gcd
import numpy as np
from kd import mul, conj, unit_from, norm1, R
from tabu import tabucol

HERE = os.path.dirname(os.path.abspath(__file__))
SC = os.path.dirname(HERE)
KISSAT = f"{SC}/kissat/build/kissat"
tag = sys.argv[1]
MAXN = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
TL = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
NGEN = int(sys.argv[4]) if len(sys.argv) > 4 else 2
RAD = int(sys.argv[5]) if len(sys.argv) > 5 else 2
PER = int(sys.argv[6]) if len(sys.argv) > 6 else 200
logf = open(f"{HERE}/grow3r_{tag}.log", "a")


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


def keys(A):
    with np.errstate(over='ignore'):
        return (np.asarray(A, dtype=np.int64) * HW).sum(axis=1)


class PointSet:
    def __init__(self, P):
        self.P = np.asarray(P, dtype=np.int64).reshape(-1, 4)
        k = keys(self.P)
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
    path = f"/dev/shm/g3_{tag}.cnf"
    cl = [[3 * v + c + 1 for c in range(3)] for v in range(n)]
    cl += [[-(3 * a + c + 1), -(3 * b + c + 1)] for a, b in E.tolist() for c in range(3)]
    a0, b0 = E[0]
    cl += [[3 * int(a0) + 1], [3 * int(b0) + 2]]          # one edge fixed to colours 0, 1
    with open(path, "w") as f:
        f.write(f"p cnf {3 * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    out = subprocess.run(["nice", "-n", "19", KISSAT, f"--time={tl}", path], capture_output=True, text=True).stdout
    st = next((l for l in out.splitlines() if l.startswith("s ")), "s UNKNOWN")
    if "UNSATISFIABLE" in st:
        return False, None
    os.unlink(path)
    if "SATISFIABLE" in st:
        lits = set(int(x) for l in out.splitlines() if l.startswith("v ") for x in l.split()[1:] if int(x) > 0)
        return True, np.array([next(c for c in range(3) if 3 * v + c + 1 in lits) for v in range(n)])
    return None, None


# ---- seed: sums of at most RAD unit vectors
st_file = f"{HERE}/grow3_{tag}_state.json"  # same names as grow3q, so pipeline_q.sh works
if os.path.exists(st_file):
    stt = json.load(open(st_file))
    P = np.array(stt["P"], dtype=np.int64)
    col = np.array(stt["col"]) if stt.get("col") else None
    log(f"resumed: {len(P)} points")
else:
    P = np.zeros((1, 4), dtype=np.int64)
    for _ in range(RAD):
        P = np.unique(np.concatenate([P] + [P + u for u in U]), axis=0)
    col = None
T0 = time.time()
rnd = 0
soft = 0
RECOL = int(os.environ.get("RECOL", "4"))
SOFT = int(os.environ.get("SOFT", "6"))
while True:
    rnd += 1
    n = len(P)
    E = build_edges(P)
    t1 = time.time()
    init = -np.ones(n, dtype=np.int64)
    if col is not None:
        init[:len(col)] = col
    c2 = None
    for s_ in range(6):
        c2 = tabucol(n, E, init=init, k=3, maxiter=2_000_000, seed=100 * rnd + s_ + 1)
        if c2 is not None:
            how = f"tabu/{s_ + 1}"; break
    if c2 is None:
        how = "kissat"
        r, c2 = kissat3(n, E, TL)
        if r is False:
            log(f"round {rnd}: n={n} e={len(E)}: kissat UNSAT -- NOT 3-COLOURABLE  [total {time.time() - T0:.0f}s]")
            np.save(f"{HERE}/grow3_{tag}_unsat_pts.npy", P)
            json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "status": "unsat"}, open(st_file, "w"))
            break
        if r is None:
            log(f"round {rnd}: n={n} e={len(E)}: kissat TIME-OUT  [total {time.time() - T0:.0f}s]")
            json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "col": None, "status": "timeout"}, open(st_file, "w"))
            break
    col = np.asarray(c2)
    assert np.all(col[E[:, 0]] != col[E[:, 1]])
    json.dump({"D": D, "U": U.tolist(), "P": P.tolist(), "col": col.tolist(), "status": "running"}, open(st_file, "w"))
    if n >= MAXN:
        log(f"round {rnd}: n={n}: 3-colourable, reached maxn"); break
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
    blocked = np.nonzero(mask == 7)[0]
    if len(blocked) == 0:
        # stuck: recolour from scratch (random tabu seeds, then kissat) and look again; then grow softly
        found = False
        for r_ in range(RECOL):
            c3 = tabucol(n, E, init=None, k=3, maxiter=2_000_000, seed=7919 * rnd + r_ + 1)
            if c3 is None:
                continue
            m3 = np.zeros(len(cand), dtype=np.int64)
            for u in U:
                idx = C.index(P + u); ok = idx >= 0
                np.bitwise_or.at(m3, idx[ok], 1 << c3[ok])
            b3 = np.nonzero(m3 == 7)[0]
            if len(b3):
                col, mask, blocked, found = np.asarray(c3), m3, b3, True
                log(f"round {rnd}: recoloured from scratch (try {r_ + 1}): {len(b3)} blocked")
                break
        if not found:
            soft += 1
            two = np.nonzero(np.isin(mask, [3, 5, 6]) & (deg >= 2))[0]
            if soft > SOFT or len(two) == 0:
                log(f"round {rnd}: n={n}: no blocked candidate after {RECOL} recolourings and {soft - 1} soft rounds; stop"); break
            blocked = two
            log(f"round {rnd}: soft round {soft}: adding candidates whose neighbours see two colours")
    pick = blocked[np.argsort(-deg[blocked], kind='stable')[:PER]]
    log(f"round {rnd}: n={n} e={len(E)} 3-col by {how} in {time.time() - t1:.1f}s; cand {len(cand)}, blocked "
        f"{len(blocked)}; add {len(pick)} (deg {deg[pick].min()}-{deg[pick].max()})  [total {time.time() - T0:.0f}s]")
    P = np.concatenate([P, cand[pick]])
    col = np.concatenate([col, -np.ones(len(pick), dtype=np.int64)])
