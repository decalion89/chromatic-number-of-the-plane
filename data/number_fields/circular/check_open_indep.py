"""Independent exact checker for the open-interval certificates (cert_sqrt{d}_open_*.json.gz).

Written by an internal referee from scratch (no code shared with check_open.py / cert_open.py / nf.py / gen_kappa.py),
then generalised from d = 11 to the d stored in the certificate.

Usage: check_open_indep.py CERT [d,N,P,Q,n1,n2,...]   (the optional list is the claimed configuration:
U = G_N {1, u_n, conj u_n : n in the list}, u_n = (n + i sqrt d)^2/(n^2 + d), one unit per +- pair); without it
only the generic checks A, C, E, F are made.

Statement checked: there is no f in the OPEN box (Q/P, 1 - Q/P)^U with <r, f> in Z for every relation r of the
certificate.  Since for a character xi of ZU the lifts f_u of xi(u) to (Q/P, 1 - Q/P) satisfy <r, f> = xi(sum r_u u)
= xi(0) = 0 mod 1, this shows that no character of ZU maps U into the open interval.

Checks (all in integer arithmetic; rationals are scaled to integers before any comparison):
  A. format: exact key sets, Python ints only (no floats/bools), canonical decimal strings for kids keys and leaf entries;
  B. (d, N, P, Q) equal to the claimed values;
  C. every unit [x0, x1, y0, y1] has x0^2 + d x1^2 + y0^2 + d y1^2 = D^2 and x0 x1 + y0 y1 = 0
     (so ((x0 + x1 s)/D)^2 + ((y0 + y1 s)/D)^2 = 1 with s^2 = d, d not a square);
  D. the units are exactly G_N * {1, u_n, cu_n} up to sign, one per +- pair, pairwise distinct up to sign,
     rebuilt here in an arithmetic of Q(sqrt d)(i) written for this file;
  E. every relation is a NONZERO integer vector with sum_u r_u u = 0 (all four integer coordinates);
  F. tree: every branch node on relation j (not yet fixed on the path) has children exactly for the integers z with
     P*min < P*z < P*max, where min/max are the extremes of <r_j, f> on the CLOSED box (computed exactly as integers
     after multiplying by P); every leaf has a rational vector y (one entry per fixed equation on the path) with
     c = sum_k y_k r_{j_k}: if c = 0 then <y, z> != 0, else <y, z> >= max_closed <c, f> or <= min_closed <c, f>.
"""
import gzip, json, re, sys
from fractions import Fraction as Fr
from math import gcd, isqrt

class Reject(Exception):
    pass

def req(cond, msg):
    if not cond:
        raise Reject(msg)

def is_int(x):
    return type(x) is int          # excludes bool and float

INT_RE = re.compile(r"-?(?:0|[1-9][0-9]*)")
RAT_RE = re.compile(r"(-?(?:0|[1-9][0-9]*))(?:/([1-9][0-9]*))?")

def parse_rat(v):
    if is_int(v):
        return Fr(v)
    req(type(v) is str, f"leaf entry of type {type(v)}")
    m = RAT_RE.fullmatch(v)
    req(m is not None, f"leaf entry not a canonical rational: {v!r}")
    num = int(m.group(1)); den = int(m.group(2)) if m.group(2) else 1
    req(den > 0, "zero denominator")
    return Fr(num, den)

# ---------- arithmetic in L = Q(s)(i), s^2 = DD: element = (a, b, c, e) = (a + b s) + i (c + e s) ----------
DD = None    # set from the certificate
def qmul(x, y):   # (x0 + x1 s)(y0 + y1 s)
    return (x[0] * y[0] + DD * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def qadd(x, y):
    return (x[0] + y[0], x[1] + y[1])
def qneg(x):
    return (-x[0], -x[1])
def lmul(z, w):   # (A + iB)(C + iE) = (AC - BE) + i(AE + BC)
    A, B = (z[0], z[1]), (z[2], z[3]); C, E = (w[0], w[1]), (w[2], w[3])
    re_ = qadd(qmul(A, C), qneg(qmul(B, E))); im_ = qadd(qmul(A, E), qmul(B, C))
    return (re_[0], re_[1], im_[0], im_[1])
def lnorm(z):     # x^2 + y^2 with x = a + b s, y = c + e s
    A, B = (z[0], z[1]), (z[2], z[3])
    return qadd(qmul(A, A), qmul(B, B))

def expected_units(D, N, ns):
    """G_N * V up to sign, V = {1, u_n, conj u_n : n in ns}, as a set of frozensets {v, -v} of integer 4-tuples
    (scaled by D)."""
    one = (Fr(1), Fr(0), Fr(0), Fr(0))
    V = [one]
    for n in ns:
        den = n * n + DD
        u = (Fr(n * n - DD, den), Fr(0), Fr(0), Fr(2 * n, den))     # (n + i s)^2 / (n^2 + d)
        cu = (u[0], u[1], -u[2], -u[3])                              # complex conjugate
        V += [u, cu]
    for v in V:
        assert lnorm(v) == (1, 0)
    G = [(Fr(x, N), Fr(0), Fr(y, N), Fr(0)) for x in range(-N, N + 1) for y in range(-N, N + 1) if x * x + y * y == N * N]
    out = set()
    for g in G:
        for v in V:
            w = lmul(g, v)
            assert lnorm(w) == (1, 0)
            W = tuple(x * D for x in w)
            if not all(x.denominator == 1 for x in W):
                raise Reject(f"D = {D} is not a common denominator of G_N V")
            W = tuple(int(x) for x in W)
            out.add(frozenset({W, tuple(-x for x in W)}))
    return out, len(G) * len(V)

def check(path, config=None, verbose=True):
    global DD
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        C = json.load(fh)
    # A/B: format and parameters
    req(type(C) is dict and set(C.keys()) == {"d", "N", "D", "P", "Q", "units", "relations", "tree"}, "top-level keys")
    d, N, D, P, Q = C["d"], C["N"], C["D"], C["P"], C["Q"]
    req(all(is_int(x) for x in (N, D, P, Q)), "parameters must be ints")
    BIQ = d == "3,11"         # the plane over Q(sqrt3, sqrt11): 8 coordinates over (1, r3, r11, r33)
    req(BIQ or (is_int(d) and d >= 2 and isqrt(d) ** 2 != d), "d not a non-square")
    req(not (BIQ and config is not None), "no configuration check over Q(sqrt3, sqrt11)")
    req(D >= 1, "D must be positive")
    req(Q >= 1 and 2 * Q <= P, "need P/Q >= 2")
    DD = None if BIQ else d
    if config is not None:
        req((d, N, P, Q) == tuple(config[:4]), f"parameters (d, N, P, Q) = {(d, N, P, Q)} differ from the claimed {tuple(config[:4])}")
    U, R = C["units"], C["relations"]
    req(type(U) is list and len(U) > 0 and type(R) is list, "units/relations must be lists")
    n = len(U)
    # C: unit vectors
    def bmul(a, b):    # Q(sqrt3, sqrt11): r3 r11 = r33, r3 r33 = 3 r11, r11 r33 = 11 r3
        return (a[0]*b[0] + 3*a[1]*b[1] + 11*a[2]*b[2] + 33*a[3]*b[3],
                a[0]*b[1] + a[1]*b[0] + 11*a[2]*b[3] + 11*a[3]*b[2],
                a[0]*b[2] + a[2]*b[0] + 3*a[1]*b[3] + 3*a[3]*b[1],
                a[0]*b[3] + a[3]*b[0] + a[1]*b[2] + a[2]*b[1])
    DIM = 8 if BIQ else 4
    for u in U:
        req(type(u) is list and len(u) == DIM and all(is_int(x) for x in u), f"bad unit {u}")
        if BIQ:
            xx, yy = bmul(u[:4], u[:4]), bmul(u[4:], u[4:])
            req(tuple(p + q for p, q in zip(xx, yy)) == (D * D, 0, 0, 0), f"length != 1: {u}")
            continue
        x0, x1, y0, y1 = u
        req(x0 * x0 + d * x1 * x1 + y0 * y0 + d * y1 * y1 == D * D, f"length != 1: {u}")
        req(x0 * x1 + y0 * y1 == 0, f"irrational part of length nonzero: {u}")
    # D: the configuration
    pairs = [frozenset({tuple(u), tuple(-x for x in u)}) for u in U]
    req(len(set(pairs)) == n, "two units are equal or opposite")
    if config is not None:
        exp, count = expected_units(D, N, config[4:])
        req(len(exp) == n and 2 * n == count, f"expected {count // 2} pairs (from {count} products), got {len(exp)} distinct / {n} given")
        req(set(pairs) == exp, "units are not G_N * V up to sign")
    # E: relations
    for r in R:
        req(type(r) is list and len(r) == n and all(is_int(x) for x in r), "bad relation format")
        req(any(x != 0 for x in r), "zero relation")
        for k in range(DIM):
            req(sum(r[i] * U[i][k] for i in range(n) if r[i]) == 0, "not a relation")
    m = len(R)
    # F: the tree, iteratively
    nodes = leaves = 0
    maxdepth = 0
    tight_leaves = 0      # leaves where <y,z> equals the max or min exactly (expected: the closed box is feasible)
    stack = [(C["tree"], ())]
    while stack:
        t, fixed = stack.pop()
        nodes += 1
        maxdepth = max(maxdepth, len(fixed))
        req(type(t) is dict, "node is not an object")
        keys = set(t.keys())
        if keys == {"leaf"}:
            ys = t["leaf"]
            req(type(ys) is list and len(ys) == len(fixed) and len(fixed) >= 1, "leaf vector length")
            y = [parse_rat(v) for v in ys]
            L = 1
            for v in y:
                L = L * v.denominator // gcd(L, v.denominator)
            Y = [int(v * L) for v in y]            # exact: L is a multiple of every denominator
            cvec = [0] * n
            tz = 0
            for yk, (j, z) in zip(Y, fixed):
                if yk:
                    rj = R[j]
                    for i in range(n):
                        if rj[i]:
                            cvec[i] += yk * rj[i]
                    tz += yk * z
            if all(ci == 0 for ci in cvec):
                req(tz != 0, f"leaf with c = 0 and <y,z> = 0 at {fixed}")
            else:
                Pmax = sum((P - Q) * ci if ci > 0 else Q * ci for ci in cvec)    # P * max over closed box
                Pmin = sum(Q * ci if ci > 0 else (P - Q) * ci for ci in cvec)    # P * min over closed box
                req(P * tz >= Pmax or P * tz <= Pmin, f"leaf inequality fails at {fixed}")
                if P * tz == Pmax or P * tz == Pmin:
                    tight_leaves += 1
            leaves += 1
            continue
        req(keys == {"branch", "lo", "hi", "kids"}, f"bad node keys {keys}")
        j = t["branch"]
        req(is_int(j) and 0 <= j < m, "branch index out of range")
        req(all(j != a for a, _ in fixed), "relation branched twice on a path")
        r = R[j]
        A = sum(Q * ci if ci > 0 else (P - Q) * ci for ci in r)     # P * min
        B = sum((P - Q) * ci if ci > 0 else Q * ci for ci in r)     # P * max
        req(A < B, "degenerate range (zero relation)")
        zlo = A // P + 1                 # least z with P z > A
        zhi = -((-B) // P) - 1           # greatest z with P z < B
        need_set = set(range(zlo, zhi + 1))
        req(is_int(t["lo"]) and is_int(t["hi"]) and (t["lo"], t["hi"]) == (zlo, zhi), f"lo/hi mismatch at {fixed}")
        kids = t["kids"]
        req(type(kids) is dict, "kids not an object")
        got = []
        for k in kids:
            req(INT_RE.fullmatch(k) is not None, f"non-canonical kid key {k!r}")
            got.append(int(k))
        req(len(got) == len(set(got)) and set(got) == need_set,
            f"children {sorted(got)} != required integers {sorted(need_set)} at {fixed}")
        for k, kid in kids.items():
            stack.append((kid, fixed + ((j, int(k)),)))
    res = dict(units=n, relations=m, nodes=nodes, leaves=leaves, maxdepth=maxdepth, tight_leaves=tight_leaves, D=D)
    if verbose:
        conf = f"(= G_{N} V up to sign, V from n = {config[4:]})" if config is not None else "(configuration not checked)"
        print(f"ACCEPTED {path}: d = {d}, N = {N}, D = {D}, P/Q = {P}/{Q}: {n} units {conf}, "
              f"{m} nonzero relations, {nodes} nodes, {leaves} leaves, max depth {maxdepth}, "
              f"{tight_leaves} leaves tight (equality)")
    return res

if __name__ == "__main__":
    try:
        check(sys.argv[1], config=[int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else None)
    except Reject as e:
        print("REJECTED:", e)
        sys.exit(1)
