"""Level-m 2-adic gate for d = 3 (mod 8). The place P above 2 in K(i), K = Q(sqrt d), is unique (residue F_4, e = 2),
so every unit vector is P-integral, and reduction modulo 2^m is a graph homomorphism from any graph built with a
direction set U to Cay(O/2^m, U mod 2^m), O = Z_4[i], Z_4 = Z_2[w], w^2 + w + 1 = 0.
In L = Q_2(sqrt d, i): sqrt(-d) = x (1 + 2w) with x^2 = d/3 in Z_2, and a point ((a + b r) + i (c + e r))/D,
r = sqrt d = -i sqrt(-d), is z = (A + i B)/D with A = a + t e, B = c - t b, t = sqrt(-d).
If the Cayley graph is 3-colourable, every graph built from U is 3-colourable."""
import sys, subprocess, os, itertools
from units_fast import units_fast
SC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KISSAT = f"{SC}/kissat/build/kissat"


def sqrt2adic(a, k):                # a = 1 mod 8
    x = 1
    for j in range(3, k + 1):
        if (x * x - a) % (2 ** j):
            x += 2 ** (j - 2)
    return x % 2 ** k


def mul(p, q, M):                   # elements (a0, a1, b0, b1) = (a0 + a1 w) + i (b0 + b1 w)
    def m4(x, y):                   # Z_4 product: (x0 + x1 w)(y0 + y1 w), w^2 = -1 - w
        return ((x[0] * y[0] - x[1] * y[1]) % M, (x[0] * y[1] + x[1] * y[0] - x[1] * y[1]) % M)
    A, B, C, D = p[:2], p[2:], q[:2], q[2:]
    ac, bd, ad, bc = m4(A, C), m4(B, D), m4(A, D), m4(B, C)
    return ((ac[0] - bd[0]) % M, (ac[1] - bd[1]) % M, (ad[0] + bc[0]) % M, (ad[1] + bc[1]) % M)


def residues(d, D, m):
    assert d % 8 == 3
    s = 0
    Dp = D
    while Dp % 2 == 0:
        Dp //= 2; s += 1
    k = m + s + 4
    M = 2 ** k
    x = sqrt2adic((d * pow(3, -1, M)) % M, k)
    t = (x % M, (2 * x) % M)        # t = x + 2x w
    S = set()
    for a, b, c, e in units_fast(d, D):
        # A = a + t e, B = c - t b  (elements of Z_4)
        A = ((a + t[0] * e) % M, (t[1] * e) % M)
        B = ((c - t[0] * b) % M, (-t[1] * b) % M)
        z = (A[0], A[1], B[0], B[1])
        assert all(v % (2 ** s) == 0 for v in z) or s == 0, "not divisible by 2^s"
        if s:
            # (A + iB) is divisible by (1+i)^(2s) = (2i)^s, i.e. by 2^s in Z_4[i]: check and divide
            if not all(v % 2 ** s == 0 for v in z):
                raise ValueError("unexpected 2-adic valuation")
            z = tuple(v // 2 ** s for v in z)
        inv = pow(Dp, -1, 2 ** m)
        z = tuple((v * inv) % 2 ** m for v in z)
        S.add(z)
    return S


def colourable(m, S, k=3):
    Mm = 2 ** m
    elems = list(itertools.product(range(Mm), repeat=4))
    idx = {e: j for j, e in enumerate(elems)}
    n = len(elems)
    cl = [[k * v + c + 1 for c in range(k)] for v in range(n)]
    for e in elems:
        for z in S:
            f = tuple((e[j] + z[j]) % Mm for j in range(4))
            a, b = idx[e], idx[f]
            if a < b:
                cl += [[-(k * a + c + 1), -(k * b + c + 1)] for c in range(k)]
    cl.append([1])                  # vertex 0 gets colour 0
    path = f"/dev/shm/gate2_{os.getpid()}.cnf"
    with open(path, "w") as fh:
        fh.write(f"p cnf {k * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    r = subprocess.run([KISSAT, "--time=600", path], capture_output=True, text=True)
    os.unlink(path)
    return {10: "yes (SAT)", 20: "no (UNSAT)"}.get(r.returncode, f"unknown ({r.returncode})")


if __name__ == "__main__":
    m = int(sys.argv[1])
    for pair in sys.argv[2:]:
        d, D = map(int, pair.split(":"))
        S = residues(d, D, m)
        print(f"d={d} D={D} m={m}: {len(units_fast(d, D))} units -> {len(S)} residues mod 2^{m}; "
              f"Cay 3-colourable? {colourable(m, S)}", flush=True)
