"""check_open.py CERT.json: independent exact check (integers and fractions only) that no character of ZU maps every
unit vector of U into the OPEN interval (Q/P, 1 - Q/P), for U in the plane over Q(sqrt d) or over Q(sqrt3, sqrt11).
A unit is an integer vector [x0, x1, y0, y1] meaning ((x0 + x1 sqrt d)/D, (y0 + y1 sqrt d)/D), or, for d = "3,11",
[x0, x1, x2, x3, y0, y1, y2, y3] meaning ((x0 + x1 s3 + x2 s11 + x3 s33)/D, (y0 + ...)/D), s_k = sqrt k.  Checks:
1. x0^2 + d x1^2 + y0^2 + d y1^2 = D^2 and x0 x1 + y0 y1 = 0 (length 1), d not a square (for d = "3,11": the
   square of the length, computed in Q(sqrt3, sqrt11), is 1);
2. every relation r is nonzero and has sum_u r_u u = 0 (a zero relation would have an empty open range);
3. each branch on relation j covers every integer strictly between the minimum and the maximum of <r_j, f> over the
   closed box [Q/P, 1 - Q/P]^U (the values <r_j, f> takes on the open box), and each leaf's vector y satisfies:
   c = sum_k y_k r_{j_k} is zero and <y, z> != 0, or <y, z> >= max_{closed box} <c, f> or <= min_{closed box} <c, f>;
   then no f in the open box satisfies the fixed equations.
Conclusion: every character xi of ZU has some u with ||xi(u)|| <= Q/P; by Theorem W+, Cay(ZU, U) has no homomorphism
to K_{p/q} for p/q < P/Q, so chi_c >= P/Q (P/Q <= 4; for P/Q = 4, kappa(U) <= 1/4)."""
import gzip, json, math, sys, re
from fractions import Fraction as Fr

class CertificateError(ValueError):
    pass

def need(c, m):
    if not c: raise CertificateError(m)

def isint(x): return type(x) is int

def check(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        C = json.load(fh)
    need(set(C) == {"d", "N", "D", "P", "Q", "units", "relations", "tree"}, "bad keys")
    d, D, P, Q, U, R = C["d"], C["D"], C["P"], C["Q"], C["units"], C["relations"]
    need(isint(P) and isint(Q) and 0 < 4 * Q and 2 * Q < P <= 4 * Q, "need 2 < P/Q <= 4")
    need(isint(D) and D >= 1, "bad D")
    LO, HI = Fr(Q, P), 1 - Fr(Q, P)
    n = len(U)
    if d == "3,11":
        dim = 8

        def sq(a0, a1, a2, a3):     # square in Q(sqrt3, sqrt11): s3 s11 = s33, s3 s33 = 3 s11, s11 s33 = 11 s3
            return (a0 * a0 + 3 * a1 * a1 + 11 * a2 * a2 + 33 * a3 * a3, 2 * a0 * a1 + 22 * a2 * a3,
                    2 * a0 * a2 + 6 * a1 * a3, 2 * a0 * a3 + 2 * a1 * a2)
        for u in U:
            need(type(u) is list and len(u) == 8 and all(isint(x) for x in u), "bad unit")
            a, b = sq(*u[:4]), sq(*u[4:])
            need(tuple(p + q for p, q in zip(a, b)) == (D * D, 0, 0, 0), f"not a unit: {u}")
    else:
        dim = 4
        need(isint(d) and d >= 2 and math.isqrt(d) ** 2 != d, "d must be a non-square")
        for u in U:
            need(type(u) is list and len(u) == 4 and all(isint(x) for x in u), "bad unit")
            x0, x1, y0, y1 = u
            need(x0 * x0 + d * x1 * x1 + y0 * y0 + d * y1 * y1 == D * D and x0 * x1 + y0 * y1 == 0, f"not a unit: {u}")
    for r in R:
        need(type(r) is list and len(r) == n and all(isint(x) for x in r), "bad relation")
        need(any(r), "zero relation")   # a zero relation has an empty open range and would close any node
        need(all(sum(r[i] * U[i][c] for i in range(n)) == 0 for c in range(dim)), "not a relation")
    nodes = leaves = 0
    stack = [(C["tree"], [])]
    while stack:
        t, fixed = stack.pop(); nodes += 1
        if "leaf" in t:
            need(set(t) == {"leaf"}, "bad leaf")
            y = [Fr(v) for v in t["leaf"]]
            for v in t["leaf"]:
                need(isint(v) or (type(v) is str and re.fullmatch(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?", v)), "bad entry")
            need(len(y) == len(fixed) >= 1, "leaf length")
            c = [sum(y[k] * R[j][i] for k, (j, _) in enumerate(fixed)) for i in range(n)]
            tv = sum(y[k] * z for k, (_, z) in enumerate(fixed))
            if all(ci == 0 for ci in c):
                need(tv != 0, f"bad leaf at {fixed}")
            else:
                mx = sum(ci * (HI if ci > 0 else LO) for ci in c); mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
                need(tv >= mx or tv <= mn, f"bad leaf at {fixed}")
            leaves += 1
            continue
        need(set(t) == {"branch", "lo", "hi", "kids"}, "bad branch node")
        j = t["branch"]
        need(isint(j) and 0 <= j < len(R) and j not in [a for a, _ in fixed], "bad branch index")
        r = R[j]
        mn = sum(ci * (LO if ci > 0 else HI) for ci in r); mx = sum(ci * (HI if ci > 0 else LO) for ci in r)
        lo, hi = math.floor(mn) + 1, math.ceil(mx) - 1
        need((t["lo"], t["hi"]) == (lo, hi), f"bad range at {fixed}")
        keys = sorted(int(z) for z in t["kids"])
        need(keys == list(range(lo, hi + 1)), f"branch does not cover its range at {fixed}")
        for z, kid in t["kids"].items():
            stack.append((kid, fixed + [(j, int(z))]))
    return (f"VERIFIED: d = {d}, D = {D}: {n} units, {len(R)} relations, {nodes} nodes, {leaves} leaves; no character "
            f"maps U into the open interval ({Q}/{P}, {P - Q}/{P}); chi_c(Cay(ZU, U)) >= {P}/{Q}"
            + ("" if 4 * Q != P else "; kappa(U) <= 1/4"))

if __name__ == "__main__":
    try:
        print(check(sys.argv[1]))
    except CertificateError as e:
        print("REJECTED:", e); sys.exit(1)
