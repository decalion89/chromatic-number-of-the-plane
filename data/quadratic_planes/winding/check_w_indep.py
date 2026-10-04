"""Second checker for the winding certificates (cert_*.json.gz here and in family/), written by a referee of the
winding paper from the certificate format alone; it shares no code with check_w.py. Every failed check raises
ValueError (no assert), so python -O checks the same. usage: python3 check_w_indep.py CERT.json.gz ...

For each certificate:
  1. every unit (a,b,c,e) satisfies a^2+d b^2+c^2+d e^2 = D^2 and ab+ce = 0 (so it is a unit vector of Q(sqrt d)^2
     with denominator D); no two units are equal or opposite;
  2. every relation is an integer relation (sum r_u u = 0 coordinatewise);
  3. tree: at each branch the children are exactly the integers in [ceil(min), floor(max)] of <r_j, f> over the box
     [1/3,2/3]^U; at each leaf, with c = sum_k y_k r_{j_k}, the value sum_k y_k z_k lies outside [min, max] of <c,f>.
     Arithmetic is done in integers after scaling by 3 (f = g/3 with g in [1,2]) and clearing denominators of y.
  4. rank of the unit matrix (over Q) and rank of the relation matrix; whether the relation lattice is saturated
     (all elementary divisors 1), i.e. the listed relations generate ALL integer relations among the units.
  5. completeness: the set of units equals the set of ALL unit vectors with denominator D, up to sign.
"""
import gzip, json, math, sys
from fractions import Fraction
import flint


def _req(cond, msg):
    if not cond:
        raise ValueError(msg)


def load(path):
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        return json.load(fh)


def parse_frac(v):
    if isinstance(v, bool):
        raise ValueError("bool")
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str):
        return Fraction(v)
    raise ValueError(f"non-exact entry {v!r}")


def all_units(d, D):
    """All (a,b,c,e) in Z^4 with a^2+d b^2+c^2+d e^2 = D^2, ab+ce = 0."""
    out = set()
    B = math.isqrt(D * D // d)
    for b in range(-B, B + 1):
        for e in range(-B, B + 1):
            M = D * D - d * (b * b + e * e)
            if M < 0:
                continue
            if b == 0 and e == 0:
                # a^2 + c^2 = D^2
                for a in range(-D, D + 1):
                    c2 = M - a * a
                    if c2 < 0:
                        continue
                    c = math.isqrt(c2)
                    if c * c == c2:
                        out.add((a, 0, c, 0))
                        out.add((a, 0, -c, 0))
                continue
            g = math.gcd(b, e)
            bp, ep = b // g, e // g
            nn = bp * bp + ep * ep
            if M % nn:
                continue
            k2 = M // nn
            k = math.isqrt(k2)
            if k * k != k2:
                continue
            for kk in {k, -k}:
                a, c = -kk * ep, kk * bp
                _req(a * b + c * e == 0, "check failed")
                out.add((a, b, c, e))
    return out


def canon(u):
    u = tuple(u)
    m = tuple(-x for x in u)
    return max(u, m)


def check(path):
    C = load(path)
    d, D, U, R, tree = C["d"], C["D"], C["units"], C["relations"], C["tree"]
    isint = lambda x: isinstance(x, int) and not isinstance(x, bool)
    _req(isint(d) and d > 0 and isint(D) and D > 0, "d and D must be positive integers")
    _req(all(len(u) == 4 and all(isint(x) for x in u) for u in U), "units must be integer 4-tuples")
    _req(all(all(isint(x) for x in r) for r in R), "relations must be integer vectors")
    n = len(U)
    for (a, b, c, e) in U:
        _req(a * a + d * b * b + c * c + d * e * e == D * D, "length")
        _req(a * b + c * e == 0, "orthogonality")
    cs = [canon(u) for u in U]
    _req(len(set(cs)) == n, "duplicate or opposite units")
    for r in R:
        _req(len(r) == n, "check failed")
        for k in range(4):
            _req(sum(r[i] * U[i][k] for i in range(n)) == 0, "relation fails")
    # ranks and saturation
    Um = flint.fmpz_mat([[int(x) for x in u] for u in U])
    rankU = Um.rank()
    Rm = flint.fmpz_mat([[int(x) for x in r] for r in R])
    rankR = Rm.rank()
    # saturation: Smith normal form of R has all nonzero invariant factors equal to 1
    snf = Rm.snf()
    diag = [snf[i, i] for i in range(min(snf.nrows(), snf.ncols()))]
    nonzero = [x for x in diag if x != 0]
    saturated = all(abs(int(x)) == 1 for x in nonzero)
    # completeness of U_D (informational; skipped, reported as None, when D/sqrt(d) > 3000)
    full = {canon(u) for u in all_units(d, D)} if math.isqrt(D * D // d) <= 3000 else None
    complete = (full == set(cs)) if full is not None else None
    # tree check, scaled: f_i = g_i / 3, g_i in [1, 2].  <r, f> = <r, g>/3.
    def rng3(r):  # min and max of <r, g> over g in [1,2]^n
        lo = sum(x * (1 if x > 0 else 2) for x in r)
        hi = sum(x * (2 if x > 0 else 1) for x in r)
        return lo, hi
    nodes = leaves = 0
    maxdepth = 0
    stack = [(tree, ())]
    while stack:
        t, fixed = stack.pop()
        nodes += 1
        maxdepth = max(maxdepth, len(fixed))
        if "leaf" in t:
            y = [parse_frac(v) for v in t["leaf"]]
            _req(len(y) == len(fixed) and len(fixed) > 0, "check failed")
            L = 1
            for q in y:
                L = L * q.denominator // math.gcd(L, q.denominator)
            yi = [int(q * L) for q in y]
            cvec = [0] * n
            rhs = 0
            for k, (j, z) in enumerate(fixed):
                if yi[k]:
                    rj = R[j]
                    for i in range(n):
                        cvec[i] += yi[k] * rj[i]
                    rhs += yi[k] * z
            lo, hi = rng3(cvec)  # range of 3 <cvec, f>
            _req(3 * rhs < lo or 3 * rhs > hi, f"leaf not infeasible at {fixed}")
            leaves += 1
            continue
        j = t["branch"]
        _req(isinstance(j, int) and 0 <= j < len(R), "check failed")
        _req(j not in [jj for jj, _ in fixed], "check failed")
        lo, hi = rng3(R[j])
        zlo = -((-lo) // 3)  # ceil(lo/3)
        zhi = hi // 3  # floor(hi/3)
        keys = sorted(int(k) for k in t["kids"])
        _req(all(str(int(k)) == k for k in t["kids"]), "noncanonical key")
        _req(keys == list(range(zlo, zhi + 1)), f"coverage fails at {fixed}")
        for k, kid in t["kids"].items():
            stack.append((kid, fixed + ((j, int(k)),)))
    return dict(d=d, D=D, units=n, relations=len(R), nodes=nodes, leaves=leaves, maxdepth=maxdepth,
                rankU=rankU, rankR=rankR, saturated=saturated, complete_UD=complete, all_UD=(len(full) if full is not None else None))


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print(p, check(p), flush=True)
