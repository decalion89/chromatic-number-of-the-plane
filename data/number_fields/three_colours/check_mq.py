"""check_mq.py CERT.json[.gz]: independent exact check of a winding certificate over a real multiquadratic field
F = Q(sqrt a_1, ..., sqrt a_k) (integers and fractions only, no solver, no field library).
A unit is an integer vector (X_b)_b, (Y_b)_b over the basis e_b = prod_j sqrt(a_j)^(b_j), b in {0,1}^k, meaning the
vector (X/D, Y/D) of F^2.  Checks:
1. every unit has X^2 + Y^2 = D^2 in F (e_b e_c = prod_j a_j^((b_j + c_j) // 2) e_(b + c mod 2));
2. every relation r satisfies sum_u r_u u = 0;
3. the tree covers every integer value z_j in the box range of <r_j, f> (f in [1/3, 2/3]^U) at each branch, and
   every leaf's vector y proves that its fixed equations <r_j, f> = z_j have no solution f in the box.
Conclusion: no character of ZU sends every u into [1/3, 2/3]; by Theorem W, Cay(ZU, U) is not 3-colourable, and
chi(F^2) >= 4."""
import gzip, json, math, re, sys
from fractions import Fraction as Fr

LO, HI = Fr(1, 3), Fr(2, 3)
_INT = re.compile(r"-?(0|[1-9][0-9]*)")
_FRAC = re.compile(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?")


class CertificateError(ValueError):
    pass


def _need(cond, msg):
    if not cond:
        raise CertificateError(msg)


def _is_int(x):
    return type(x) is int


def _squarefree(a):
    return a >= 2 and all(a % (p * p) for p in range(2, math.isqrt(a) + 1))


def square(gens, basis, X):
    """X^2 in the multiquadratic basis, as a dict basis-tuple -> integer"""
    out = {}
    for b, xb in zip(basis, X):
        if not xb: continue
        for c, xc in zip(basis, X):
            if not xc: continue
            coef = xb * xc
            for j, a in enumerate(gens):
                if b[j] + c[j] == 2:
                    coef *= a
            key = tuple((bj + cj) % 2 for bj, cj in zip(b, c))
            out[key] = out.get(key, 0) + coef
    return out


def check(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        C = json.load(fh)
    _need(type(C) is dict and set(C) == {"gens", "basis", "N", "D", "units", "relations", "tree"}, "bad keys")
    gens, basis, D, U, R = C["gens"], C["basis"], C["D"], C["units"], C["relations"]
    _need(type(gens) is list and gens and all(_is_int(a) and _squarefree(a) for a in gens), "bad generators")
    k = len(gens)
    _need(type(basis) is list and sorted(tuple(b) for b in basis) == sorted(
        tuple((m >> j) & 1 for j in range(k)) for m in range(2 ** k)), "basis is not {0,1}^k")
    basis = [tuple(b) for b in basis]
    # the generators must be independent modulo squares, so that F has degree 2^k
    for m in range(1, 2 ** k):
        prod = 1
        for j in range(k):
            if (m >> j) & 1: prod *= gens[j]
        _need(math.isqrt(prod) ** 2 != prod, f"generators dependent modulo squares (mask {m})")
    _need(_is_int(D) and D >= 1, "bad denominator")
    nb = len(basis)
    _need(type(U) is list and len(U) >= 1, "no units")
    n = len(U)
    one = tuple([0] * k)
    for u in U:
        _need(type(u) is list and len(u) == 2 * nb and all(_is_int(x) for x in u), f"bad unit {u!r}")
        sx, sy = square(gens, basis, u[:nb]), square(gens, basis, u[nb:])
        tot = {key: sx.get(key, 0) + sy.get(key, 0) for key in set(sx) | set(sy)}
        _need(all(v == (D * D if key == one else 0) for key, v in tot.items()) and tot.get(one, 0) == D * D,
              f"not a unit vector: {u}")
    _need(type(R) is list and len(R) >= 1, "no relations")
    for r in R:
        _need(type(r) is list and len(r) == n and all(_is_int(x) for x in r), "bad relation")
        _need(all(sum(r[i] * U[i][c] for i in range(n)) == 0 for c in range(2 * nb)), f"not a relation: {r}")

    def rng(r):
        lo = sum(c * (LO if c > 0 else HI) for c in r)
        hi = sum(c * (HI if c > 0 else LO) for c in r)
        return math.ceil(lo), math.floor(hi)

    nodes = leaves = 0
    stack = [(C["tree"], [])]
    while stack:
        t, fixed = stack.pop()
        nodes += 1
        _need(type(t) is dict, "bad tree node")
        if "leaf" in t:
            _need(set(t) == {"leaf"} and type(t["leaf"]) is list, "bad leaf")
            y = []
            for v in t["leaf"]:
                _need(_is_int(v) or (type(v) is str and _FRAC.fullmatch(v)), f"bad leaf entry {v!r}")
                y.append(Fr(v))
            _need(len(y) == len(fixed) and len(fixed) >= 1, "leaf vector has the wrong length")
            c = [sum(y[q] * R[j][i] for q, (j, _) in enumerate(fixed)) for i in range(n)]
            mx = sum(ci * (HI if ci > 0 else LO) for ci in c)
            mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
            val = sum(y[q] * z for q, (_, z) in enumerate(fixed))
            _need(val > mx or val < mn, f"bad leaf at {fixed}")
            leaves += 1
            continue
        _need(set(t) == {"branch", "lo", "hi", "kids"} and type(t["kids"]) is dict, "bad branch node")
        j = t["branch"]
        _need(_is_int(j) and 0 <= j < len(R), "bad branch index")
        _need(j not in [jj for jj, _ in fixed], f"relation {j} branched twice")
        lo, hi = rng(R[j])
        _need((t["lo"], t["hi"]) == (lo, hi), f"bad range at {fixed}")
        kids = []
        for z, kid in t["kids"].items():
            _need(type(z) is str and _INT.fullmatch(z), f"bad branch value {z!r}")
            kids.append((int(z), kid))
        _need(sorted(z for z, _ in kids) == list(range(lo, hi + 1)), f"branch does not cover its range at {fixed}")
        for z, kid in kids:
            stack.append((kid, fixed + [(j, z)]))
    return (f"VERIFIED: F = Q(sqrt {', sqrt '.join(map(str, gens))}), D = {D}: {n} unit vectors, {len(R)} relations, "
            f"{nodes} nodes, {leaves} leaves; no character of ZU maps U into [1/3, 2/3]")


if __name__ == "__main__":
    try:
        print(check(sys.argv[1]))
    except CertificateError as err:
        print(f"REJECTED: {err}")
        sys.exit(1)
