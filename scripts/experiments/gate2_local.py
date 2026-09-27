"""The 2-adic gate for any field Q(sqrt d_1, ..., sqrt d_n) with odd d_j, at k colours.

When every d_j is odd, a place of L above 2 has completion inside K = Q_2(sqrt3, sqrt5), and it splits
in L(i) as soon as i lies in it (the classes of the d_j modulo squares include 7, or both 3 and 5).
There z = x + i y goes to the pair (iota(x) + i iota(y), iota(x) - i iota(y)) in K x K, a unit vector to
(t, 1/t), and a finite unit set, scaled by a power of 2 to become integral, to a finite Cayley graph on
O_K / pi^A x O_K / pi^B (notes/local_colourings.md, sections 9 and 12).  A proper k-colouring of it
k-colours every unit-distance graph whose edge vectors are in the set.  module_gate.py skips 2; this
script is the missing case.

Model: O_K = Z_2[w][pi] with w^2 = -1 - w (unramified, residue field F_4) and pi^2 = 2 - 2 pi, so that
pi = sqrt3 - 1 is a uniformiser; Z_2-basis 1, w, pi, w pi; arithmetic modulo 2^N.
  sqrt3 = 1 + pi,  sqrt(-3) = 1 + 2w,  i = sqrt3 sqrt(-3) / 3,
and sqrt d = sqrt c * sqrt(d / c) with c in {1, 3, -1, -3} of the class of d mod 8 and d / c = 1 mod 8.

usage: gate2_local.py units.json k A,B [A,B ...] [--kissat SECONDS] [--places N]
Levels A, B are absolute: the quotient O_K/pi^A x O_K/pi^B, after scaling the units by 2^E so that
all their coordinates are 2-integral (E is printed).
"""
import os, sys, json, time, itertools, subprocess
from fractions import Fraction as Fr
import numpy as np

N = 64                     # precision: arithmetic modulo 2^N
MOD = 1 << N
KISSAT = os.environ.get("KISSAT", "kissat")


def mul(a, b):
    """product in Z_2[w, pi]/(w^2 + w + 1, pi^2 + 2 pi - 2), basis 1, w, pi, w pi (tuples mod 2^N)"""
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    # (a0 + a1 w)(b0 + b1 w) with w^2 = -1 - w
    def m2(x0, x1, y0, y1):
        p0 = x0 * y0 - x1 * y1
        p1 = x0 * y1 + x1 * y0 - x1 * y1
        return p0, p1
    c0, c1 = m2(a0, a1, b0, b1)                   # pi^0 part
    d0, d1 = m2(a0, a1, b2, b3)
    e0, e1 = m2(a2, a3, b0, b1)                   # pi^1 parts
    f0, f1 = m2(a2, a3, b2, b3)                   # pi^2 part = f (2 - 2 pi)
    r0 = c0 + 2 * f0
    r1 = c1 + 2 * f1
    r2 = d0 + e0 - 2 * f0
    r3 = d1 + e1 - 2 * f1
    return tuple(x % MOD for x in (r0, r1, r2, r3))


def add(a, b):
    return tuple((x + y) % MOD for x, y in zip(a, b))


def scal(q, a):
    """a rational q with odd denominator times a"""
    q = Fr(q)
    assert q.denominator % 2 == 1
    s = q.numerator * pow(q.denominator, -1, MOD)
    return tuple(x * s % MOD for x in a)


ONE = (1, 0, 0, 0)
SQRT3 = (1, 0, 1, 0)
SQRTM3 = (1, 2, 0, 0)
I = scal(Fr(1, 3), mul(SQRT3, SQRTM3))
assert mul(SQRT3, SQRT3) == scal(3, ONE) and mul(SQRTM3, SQRTM3) == scal(-3, ONE) and mul(I, I) == scal(-1, ONE)


def sqrt2adic(u):
    """a square root in Z_2 of a 2-adic integer u = 1 mod 8 (given mod 2^(N+3)), modulo 2^N"""
    assert u % 8 == 1
    x = 1
    for k in range(3, N + 2):                     # x^2 = u mod 2^k, lift to 2^(k+1)
        if (x * x - u) % (1 << (k + 1)):
            x += 1 << (k - 1)
    x %= MOD
    assert (x * x - u) % MOD == 0
    return x


def sqrt_of(d):
    """sqrt d in the model, d odd squarefree"""
    c = {1: 1, 3: 3, 7: -1, 5: -3}[d % 8]
    base = {1: ONE, 3: SQRT3, -1: I, -3: SQRTM3}[c]
    u = Fr(d, c)                                   # a 2-adic unit = 1 mod 8
    M = 1 << (N + 3)
    ui = u.numerator * pow(u.denominator, -1, M) % M
    assert ui % 8 == 1
    return scal(sqrt2adic(ui), base)


def embed_basis(gens, signs):
    """images of the basis products prod_{j in mask} sqrt d_j, indexed by bitmask"""
    roots = [scal(s, sqrt_of(d)) for d, s in zip(gens, signs)]
    out = []
    for m in range(1 << len(gens)):
        v = ONE
        for j in range(len(gens)):
            if m >> j & 1:
                v = mul(v, roots[j])
        out.append(v)
    return out


PI = (0, 0, 1, 0)
ONE_MINUS_PI_INV = None


def inv_unit(a):
    """inverse of a unit of O_K modulo 2^N, by Newton iteration from the inverse modulo pi"""
    # the residue field is F_4 = {0, 1, w, w^2}; find b0 with a b0 = 1 mod pi by trying 1, w, w^2 = -1 - w
    for b0 in [(1, 0, 0, 0), (0, 1, 0, 0), (MOD - 1, MOD - 1, 0, 0)]:
        r = mul(a, b0)
        if (r[0] - 1) % 2 == 0 and r[1] % 2 == 0:        # a b0 = 1 mod (2, pi)?  check mod pi below
            b = b0
            break
    for _ in range(8):                                   # b <- b (2 - a b), quadratic convergence
        b = mul(b, add(scal(2, ONE), scal(-1, mul(a, b))))
    assert mul(a, b) == ONE, "not a unit"
    return b


def div_pi(v, j):
    """v / pi^j, exact"""
    global ONE_MINUS_PI_INV
    if ONE_MINUS_PI_INV is None:
        ONE_MINUS_PI_INV = inv_unit(add(ONE, scal(-1, PI)))
    for _ in range(j):
        w = mul(mul(v, PI), ONE_MINUS_PI_INV)             # v pi / (1 - pi) = 2 v / pi
        assert all(x % 2 == 0 for x in w), "not divisible by pi"
        v = tuple((x // 2) % MOD for x in w)            # exact halving (N bits of precision drop by one)
    return v


def val(v, cap=40):
    """pi-adic valuation of v (capped)"""
    j = 0
    while j < cap:
        w = mul(mul(v, PI), inv_unit(add(ONE, scal(-1, PI))) if ONE_MINUS_PI_INV is None else ONE_MINUS_PI_INV)
        if any(x % 2 for x in w):
            return j
        v = tuple((x // 2) % MOD for x in w); j += 1
    return j


def row_hnf(rows):
    A = [list(map(int, r)) for r in rows]
    n = len(A[0]); out = []
    for c in range(n):
        cand = [r for r in A if r[c] != 0]
        rest = [r for r in A if r[c] == 0]
        while len(cand) > 1:
            cand.sort(key=lambda r: abs(r[c]))
            p = cand[0]; new = [p]
            for r in cand[1:]:
                q = r[c] // p[c]
                r2 = [a - q * b for a, b in zip(r, p)]
                (new if r2[c] != 0 else rest).append(r2)
            cand = new
        p = cand[0] if cand[0][c] > 0 else [-a for a in cand[0]]
        out.append(p); A = rest
    return out


def ideal_rows(A):
    """Z-lattice of pi^A O_K + 2^N-stuff, in the basis 1, w, pi, w pi (exact integer HNF)"""
    pa = ONE
    for _ in range(A):
        pa = mul(pa, PI)
    rows = []
    for e in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]:
        v = mul(pa, e)
        rows.append([x if x < MOD // 2 else x - MOD for x in v])
    # 2^A kills everything below pi^A? not needed: pi^A O_K has index 4^A and the rows span it exactly
    return row_hnf(rows)


def gate(gens, units, k, levels, signs, ktime=None):
    D2 = 1                                          # common 2-power of the denominators
    for u in units:
        for c in list(u[0]) + list(u[1]):
            while (c * D2).denominator % 2 == 0:
                D2 *= 2
    B = embed_basis(gens, signs)

    def image(coords):
        v = (0, 0, 0, 0)
        for c, b in zip(coords, B):
            if c:
                v = add(v, scal(c * D2, b))
        return v
    pairs = []
    for x, y in units:
        ix, iy = image(x), image(y)
        t = add(ix, mul(I, iy)); tb = add(ix, scal(-1, mul(I, iy)))
        pairs.append((t, tb))
    # divide every component by the largest common power of pi, so that the scaling is pi^m with
    # m the largest depth |v(t)| in the set (as in split_gate_q35.py)
    s = min(val(c) for pr in pairs for c in pr)
    pairs = [(div_pi(t, s), div_pi(tb, s)) for t, tb in pairs]
    E2 = (D2.bit_length() - 1, s)
    results = []
    for A, Bv in levels:
        HA = np.array(ideal_rows(A), dtype=object); HB = np.array(ideal_rows(Bv), dtype=object)
        dA = [int(HA[t][t]) for t in range(4)]; dB = [int(HB[t][t]) for t in range(4)]
        dims = dA + dB
        n = int(np.prod(dims))

        def red(v, H, d):
            v = list(v)
            for t in range(4):
                q = v[t] // d[t]
                if q:
                    v = [a - q * int(b) for a, b in zip(v, H[t])]
            return v
        gensQ = set()
        for t, tb in pairs:
            a = red([x % MOD for x in t], HA, dA)
            b = red([x % MOD for x in tb], HB, dB)
            gensQ.add(tuple(a + b))
        if any(all(x == 0 for x in g) for g in gensQ):
            print(f"  level ({A},{Bv}): a unit reduces to 0"); results.append(None); continue
        G = np.array(sorted(gensQ), dtype=np.int64)
        radix = np.cumprod([1] + dims[:-1]).astype(np.int64)
        HAi = np.array(HA.tolist(), dtype=np.int64); HBi = np.array(HB.tolist(), dtype=np.int64)
        dAi = np.array(dA, dtype=np.int64); dBi = np.array(dB, dtype=np.int64)

        def redv(V):
            V = V.copy()
            for t in range(4):
                V[:, :4] -= np.floor_divide(V[:, t], dAi[t])[:, None] * HAi[t][None, :]
            for t in range(4):
                V[:, 4:] -= np.floor_divide(V[:, 4 + t], dBi[t])[:, None] * HBi[t][None, :]
            return V
        allK = np.arange(n, dtype=np.int64)
        X = np.stack([(allK // radix[t]) % dims[t] for t in range(8)], 1)
        E = []
        for g in G:
            Y = (redv(X + g[None, :]) * radix[None, :]).sum(1)
            m = allK < Y
            E.append(np.stack([allK[m], Y[m]], 1))
        E = np.unique(np.concatenate(E), axis=0)
        del X
        print(f"  level ({A},{Bv}): {n} points, {len(G)} unit residues, {len(E)} edges", flush=True)
        if len(E) < 12_000_000:
            from pysat.solvers import Solver
            s = Solver(name="cd19")
            for v in range(n):
                s.add_clause([1 + v * k + c for c in range(k)])
            for u_, v_ in E.tolist():
                for c in range(k):
                    s.add_clause([-(1 + u_ * k + c), -(1 + v_ * k + c)])
            s.add_clause([1])
            res = s.solve(); s.delete()
        elif ktime:
            fn = f"gate2_{os.getpid()}.cnf"
            with open(fn, "w") as f:
                f.write(f"p cnf {n * k} {n + len(E) * k + 1}\n1 0\n")
                V = np.arange(n, dtype=np.int64)
                np.savetxt(f, np.stack([1 + V * k + c for c in range(k)] + [np.zeros(n, dtype=np.int64)], 1), fmt="%d")
                for c in range(k):
                    np.savetxt(f, np.stack([-(1 + E[:, 0] * k + c), -(1 + E[:, 1] * k + c),
                                            np.zeros(len(E), dtype=np.int64)], 1), fmt="%d")
            r = subprocess.run(["timeout", str(ktime), KISSAT, fn], capture_output=True, text=True)
            os.remove(fn)
            st = [l for l in r.stdout.splitlines() if l.startswith("s ")]
            res = None if not st else ("UNSAT" not in st[0])
        else:
            print("  too large for pysat; pass --kissat SECONDS"); results.append(None); continue
        print(f"    {k}-colourable: {res}", flush=True)
        results.append(res)
    return E2, results


if __name__ == "__main__":
    argv = sys.argv[1:]
    ktime = None; nplaces = None
    if "--kissat" in argv:
        i = argv.index("--kissat"); ktime = int(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    if "--places" in argv:
        i = argv.index("--places"); nplaces = int(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    d = json.load(open(argv[0]))
    gens = tuple(d["field_generators"])
    assert all(g % 2 for g in gens), "every generator must be odd"
    rd = lambda xy: ([Fr(a, b) for a, b in xy[0]], [Fr(a, b) for a, b in xy[1]])
    units = [rd(u) for u in d["units"]]
    k = int(argv[1])
    levels = [tuple(map(int, s.split(","))) for s in argv[2:]]
    classes = {g % 8 for g in gens}
    split = 7 in classes or {3, 5} <= classes or ({3, 7} <= classes) or ({5, 7} <= classes)
    print(f"{len(units)} units over Q(sqrt {gens}); classes mod 8: {sorted(classes)}; "
          f"{'split' if split else 'NOT split'} at 2", flush=True)
    seen = 0
    for signs in itertools.product((1, -1), repeat=len(gens)):
        if signs[0] == -1 and len(gens) > 0:
            continue                               # -sqrt d_1 gives a conjugate embedding of the same place
        E2, res = gate(gens, units, k, levels, signs, ktime)
        print(f"place {signs}: coordinates scaled by 2^{E2[0]}, then divided by pi^{E2[1]} "
              f"(largest depth m = {2 * E2[0] - E2[1]}); {res}", flush=True)
        seen += 1
        if nplaces and seen >= nplaces:
            break
