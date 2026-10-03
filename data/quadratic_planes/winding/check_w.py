"""check_w.py CERT.json[.gz]: independent exact check of a winding certificate (integers and fractions only, no solver).
1. every unit ((a + b sqrt d)/D, (c + e sqrt d)/D) has length 1: a^2 + d b^2 + c^2 + d e^2 = D^2 and ab + ce = 0;
2. every relation r satisfies sum_u r_u u = 0;
3. the tree covers every integer value z_j in the box range of <r_j, f> (f in [1/3, 2/3]^U) at each branch, and
   every leaf's vector y proves that its fixed equations <r_j, f> = z_j have no solution f in the box.
Conclusion: no character of ZU sends every u into [1/3, 2/3], so by Theorem W (notes/winding_lemma.md)
Cay(ZU, U) is not 3-colourable, and chi(Q(sqrt d)^2) >= 4."""
import sys, json, math, gzip
from fractions import Fraction as Fr

LO, HI = Fr(1, 3), Fr(2, 3)


def check(path):
    C = json.load(gzip.open(path, "rt") if path.endswith(".gz") else open(path))
    d, D = C["d"], C["D"]; U = C["units"]; R = C["relations"]; n = len(U)
    for a, b, c, e in U:
        assert a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0, (a, b, c, e)
    for r in R:
        assert len(r) == n and all(sum(r[i] * U[i][k] for i in range(n)) == 0 for k in range(4)), r

    def rng(r):
        lo = sum(c * (LO if c > 0 else HI) for c in r); hi = sum(c * (HI if c > 0 else LO) for c in r)
        return math.ceil(lo), math.floor(hi)

    cnt = [0, 0]
    stack = [(C["tree"], [])]
    while stack:
        t, fixed = stack.pop()
        cnt[0] += 1
        if "leaf" in t:
            y = [Fr(v) for v in t["leaf"]]
            assert len(y) == len(fixed) and fixed
            c = [sum(y[k] * R[j][i] for k, (j, _) in enumerate(fixed)) for i in range(n)]
            mx = sum(ci * (HI if ci > 0 else LO) for ci in c); mn = sum(ci * (LO if ci > 0 else HI) for ci in c)
            val = sum(y[k] * z for k, (_, z) in enumerate(fixed))
            assert val > mx or val < mn, "bad leaf"
            cnt[1] += 1
            continue
        j = t["branch"]; lo, hi = rng(R[j])
        assert 0 <= j < len(R) and j not in [jj for jj, _ in fixed] and (t["lo"], t["hi"]) == (lo, hi)
        assert sorted(int(z) for z in t["kids"]) == list(range(lo, hi + 1))
        for z, kid in t["kids"].items():
            stack.append((kid, fixed + [(j, int(z))]))
    return (f"VERIFIED: d={d} D={D}: {n} unit vectors, {len(R)} relations, {cnt[0]} nodes, {cnt[1]} leaves; "
            f"no character of ZU maps U into [1/3, 2/3]")


if __name__ == "__main__":
    print(check(sys.argv[1]))
