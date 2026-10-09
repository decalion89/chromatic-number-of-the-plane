"""check_w.py CERT.json[.gz]: independent exact check of a winding certificate (integers and fractions only, no solver).
1. every unit ((a + b sqrt d)/D, (c + e sqrt d)/D) has length 1: a^2 + d b^2 + c^2 + d e^2 = D^2 and ab + ce = 0;
2. every relation r satisfies sum_u r_u u = 0;
3. the tree covers every integer value z_j in the box range of <r_j, f> (f in [1/3, 2/3]^U) at each branch, and
   every leaf's vector y proves that its fixed equations <r_j, f> = z_j have no solution f in the box.
Conclusion: no character of ZU sends every u into [1/3, 2/3], so by Theorem W (notes/winding_lemma.md)
Cay(ZU, U) is not 3-colourable, and chi(Q(sqrt d)^2) >= 4.

Every check raises CertificateError (it does not use assert, so python -O checks the same), and only integers
(JSON numbers without a fraction part, or decimal strings for the branch values) and exact fractions written as
strings "p" or "p/q" are accepted: a float anywhere makes the certificate invalid."""
import gzip
import json
import math
import re
import sys
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
    return type(x) is int  # excludes bool and float


def _key_int(z):
    _need(type(z) is str and _INT.fullmatch(z), f"bad branch value {z!r}")
    return int(z)


def _frac(v):
    if _is_int(v):
        return Fr(v)
    _need(type(v) is str and _FRAC.fullmatch(v), f"bad leaf entry {v!r}")
    return Fr(v)


def check(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        C = json.load(fh)
    _need(type(C) is dict and set(C) == {"d", "D", "units", "relations", "tree"}, "bad top-level keys")
    d, D, U, R = C["d"], C["D"], C["units"], C["relations"]
    _need(_is_int(d) and d >= 2 and math.isqrt(d) ** 2 != d, f"d = {d!r} is not a non-square integer >= 2")
    _need(_is_int(D) and D >= 1, f"bad denominator {D!r}")
    _need(type(U) is list and len(U) >= 1, "no unit vectors")
    n = len(U)
    for u in U:
        _need(type(u) is list and len(u) == 4 and all(_is_int(x) for x in u), f"bad unit vector {u!r}")
        a, b, c, e = u
        _need(a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0, f"not a unit vector: {u}")
    _need(type(R) is list and len(R) >= 1, "no relations")
    for r in R:
        _need(type(r) is list and len(r) == n and all(_is_int(x) for x in r), f"bad relation {r!r}")
        _need(all(sum(r[i] * U[i][k] for i in range(n)) == 0 for k in range(4)), f"not a relation: {r}")

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
            _need(set(t) == {"leaf"} and type(t["leaf"]) is list, "bad leaf node")
            y = [_frac(v) for v in t["leaf"]]
            _need(len(y) == len(fixed) and len(fixed) >= 1, "leaf vector has the wrong length")
            c = [sum(y[k] * R[j][i] for k, (j, _) in enumerate(fixed)) for i in range(n)]
            mx = sum(ci * (HI if ci > 0 else LO) for ci in c)
            mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
            val = sum(y[k] * z for k, (_, z) in enumerate(fixed))
            _need(val > mx or val < mn, f"bad leaf at {fixed}")
            leaves += 1
            continue
        _need(set(t) == {"branch", "lo", "hi", "kids"} and type(t["kids"]) is dict, "bad branch node")
        j = t["branch"]
        _need(_is_int(j) and 0 <= j < len(R), f"bad branch index {j!r}")
        _need(j not in [jj for jj, _ in fixed], f"relation {j} branched twice")
        lo, hi = rng(R[j])
        _need(_is_int(t["lo"]) and _is_int(t["hi"]) and (t["lo"], t["hi"]) == (lo, hi), f"bad range at {fixed}")
        kids = [(_key_int(z), kid) for z, kid in t["kids"].items()]
        _need(sorted(z for z, _ in kids) == list(range(lo, hi + 1)), f"branch does not cover its range at {fixed}")
        for z, kid in kids:
            stack.append((kid, fixed + [(j, z)]))
    return (f"VERIFIED: d={d} D={D}: {n} unit vectors, {len(R)} relations, {nodes} nodes, {leaves} leaves; "
            f"no character of ZU maps U into [1/3, 2/3]")


if __name__ == "__main__":
    try:
        print(check(sys.argv[1]))
    except CertificateError as err:
        print(f"REJECTED: {err}")
        sys.exit(1)
