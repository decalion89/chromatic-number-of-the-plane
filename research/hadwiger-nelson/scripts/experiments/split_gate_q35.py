"""Periodic 4-colourings of a unit module over L = Q(sqrt3, sqrt5) at the split primes above 2, 3, 5.

K = L(i) = Q(zeta12, sqrt5) has ring of integers O = Z[zeta12, phi], phi = (1 + sqrt5)/2, with basis
zeta^j phi^f (j < 4, f < 2).  The primes 2, 3, 5 each have one prime p of L above them (uniformisers
sqrt3 - 1, sqrt3, sqrt5; residue fields F_4, F_9, F_25), and p splits in K as P Pbar.  A unit vector u
has v_P(u) = -v_Pbar(u), and it can be large: phi2 = (1 + i sqrt15)/4 has valuations (2, -2) at 2.

Scale the unit set U by pi^m so that every unit is integral at P and Pbar.  For A, B >= m + 1 the map
    z -> z mod P^A Pbar^B,   O -> O / P^A Pbar^B   (a group of order q^(A+B))
is a homomorphism of groups, and it sends pi^m Z[U] to the Cayley graph Q = Cay(O / P^A Pbar^B, pi^m U).
No scaled unit maps to 0, since its valuations are (m + v, m - v) with |v| <= m < A, B.  So a proper
4-colouring of Q 4-colours every unit-distance graph whose edge vectors lie in U.  Colourability is
monotone: a colouring at (A, B) pulls back to every finer level.

Levels are given relative to m: "a,b" means A = m + a, B = m + b.  Small quotients are solved with
pysat; larger ones are written as DIMACS and handed to kissat (--kissat SECONDS).

usage: split_gate_q35.py units.json p a,b [a,b ...] [--kissat SECONDS]      (p in 2, 3, 5)
Results (notes/local_colourings.md, section 12):
  zeta, tau                      p = 2: 4-colourable at a,b = 1,1 (the image is K_4)
  zeta, tau, phi2 (108 units)    p = 2: not at 1,1; 4-colourable at 2,2
  zeta, o1 (36 units)            p = 2: 4-colourable at 1,1
  zeta, o1, phi2 (84 units)      p = 2: not 4-colourable at 1,1, 2,1, 2,2 (1 048 576 points; kissat)
"""
import os, sys, json, time, subprocess
from fractions import Fraction as Fr
import numpy as np

KISSAT = os.environ.get("KISSAT", "kissat")


# ---- the ring O = Z[zeta12, phi] -------------------------------------------------------------------
def omul(a, b):
    out = [Fr(0)] * 8
    for i in range(8):
        if not a[i]:
            continue
        for j in range(8):
            if not b[j]:
                continue
            z = [Fr(0)] * 7
            z[(i % 4) + (j % 4)] = Fr(1)
            for t in range(6, 3, -1):                      # zeta^4 = zeta^2 - 1
                if z[t]:
                    c = z[t]; z[t] = Fr(0); z[t - 2] += c; z[t - 4] -= c
            ph = [Fr(0)] * 3
            ph[(i // 4) + (j // 4)] = Fr(1)
            if ph[2]:                                      # phi^2 = phi + 1
                c = ph[2]; ph[2] = Fr(0); ph[1] += c; ph[0] += c
            for zz in range(4):
                for f in range(2):
                    if z[zz] and ph[f]:
                        out[zz + 4 * f] += a[i] * b[j] * z[zz] * ph[f]
    return out


def const(q):
    v = [Fr(0)] * 8
    v[0] = Fr(q)
    return v


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sc(q, a):
    return [Fr(q) * x for x in a]


def powe(x, n):
    out = const(1)
    for _ in range(n):
        out = omul(out, x)
    return out


ZETA = [Fr(int(t == 1)) for t in range(8)]
PHI = [Fr(int(t == 4)) for t in range(8)]
I = powe(ZETA, 3)
SQRT3 = add(sc(2, ZETA), sc(-1, I))
SQRT5 = add(sc(2, PHI), const(-1))
assert omul(SQRT3, SQRT3) == const(3) and omul(SQRT5, SQRT5) == const(5) and omul(I, I) == const(-1)
BASIS_L = [const(1), SQRT3, SQRT5, omul(SQRT3, SQRT5)]     # the order of hn.field.Field((3, 5))


def from_xy(xc, yc):
    """z = x + i y for coordinates x, y on the basis 1, sqrt3, sqrt5, sqrt15"""
    z = [Fr(0)] * 8
    for c, b in zip(xc, BASIS_L):
        z = add(z, sc(c, b))
    for c, b in zip(yc, BASIS_L):
        z = add(z, sc(c, omul(I, b)))
    return z


def row_hnf(rows):
    """upper-triangular basis of the full-rank Z-lattice spanned by integer rows"""
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
        if not cand:
            raise ValueError("not full rank")
        p = cand[0] if cand[0][c] > 0 else [-a for a in cand[0]]
        out.append(p); A = rest
    return out


# uniformiser pi of L at p, and gamma with valuations (4, 0) at (P, Pbar): gamma O = P^4
PLACE = {
    2: (add(SQRT3, const(-1)),
        omul(from_xy([Fr(1, 4), 0, 0, 0], [0, 0, 0, Fr(1, 4)]), powe(add(SQRT3, const(-1)), 2)),
        omul(from_xy([Fr(1, 4), 0, 0, 0], [0, 0, 0, Fr(-1, 4)]), powe(add(SQRT3, const(-1)), 2))),
    3: (SQRT3, from_xy([2, 0, 0, 0], [0, 0, 1, 0]), from_xy([2, 0, 0, 0], [0, 0, -1, 0])),
    5: (SQRT5, from_xy([3, 0, 0, 0], [4, 0, 0, 0]), from_xy([3, 0, 0, 0], [-4, 0, 0, 0])),
}


def ideal(p, A, B):
    """HNF rows of P^A Pbar^B"""
    pi, gam, gamb = PLACE[p]
    lo = min(A, B); g = gam if A >= B else gamb
    q = (abs(A - B) + 3) // 4
    gens = [omul(powe(pi, lo), powe(g, q)), powe(pi, max(A, B))]
    rows = []
    for x in gens:
        for s in range(8):
            rows.append([int(c) for c in omul(x, [Fr(int(t == s)) for t in range(8)])])
    return row_hnf(rows)


def gate(units, p, levels, ktime=None, tag="gate"):
    pi = PLACE[p][0]
    m = 0
    while not all(c.denominator % p for u in units for c in omul(powe(pi, m), u)):
        m += 1
    su = [omul(powe(pi, m), u) for u in units]
    print(f"{len(units)} units, scaled by pi^{m} at p = {p}", flush=True)
    results = []
    for a, b in levels:
        A, B = m + a, m + b
        H = np.array(ideal(p, A, B), dtype=np.int64)
        dg = np.array([H[t, t] for t in range(8)], dtype=np.int64)
        n = int(np.prod(dg))
        radix = np.cumprod(np.concatenate([[1], dg[:-1]])).astype(np.int64)

        def red(V):
            V = V.copy()
            for t in range(8):
                V -= np.floor_divide(V[:, t], dg[t])[:, None] * H[t][None, :]
            return V
        key = lambda V: (V * radix[None, :]).sum(1)
        Nm = p ** 24
        G = np.array([[c.numerator * pow(c.denominator, -1, Nm) % Nm for c in u] for u in su], dtype=np.int64)
        G = np.unique(red(G), axis=0)
        if np.any(key(G) == 0):
            print(f"level ({A},{B}): a unit reduces to 0"); results.append(None); continue
        allK = np.arange(n, dtype=np.int64)
        X = np.stack([(allK // radix[t]) % dg[t] for t in range(8)], 1)
        E = []
        for g in G:
            Y = key(red(X + g[None, :]))
            msk = allK < Y
            E.append(np.stack([allK[msk], Y[msk]], 1))
        E = np.unique(np.concatenate(E), axis=0)
        del X
        print(f"level ({A},{B}): {n} points, {len(G)} unit residues, {len(E)} edges", flush=True)
        K = 4
        if len(E) < 12_000_000:
            from pysat.solvers import Solver
            s = Solver(name="cd19")
            for v in range(n):
                s.add_clause([1 + v * K + c for c in range(K)])
            for u, v in E.tolist():
                for c in range(K):
                    s.add_clause([-(1 + u * K + c), -(1 + v * K + c)])
            s.add_clause([1])
            res = s.solve(); s.delete()
        elif ktime:
            fn = f"{tag}_p{p}_{A}_{B}.cnf"
            with open(fn, "w") as f:
                f.write(f"p cnf {n * K} {n + len(E) * K + 1}\n1 0\n")
                V = np.arange(n, dtype=np.int64)
                np.savetxt(f, np.stack([1 + V * K + c for c in range(K)] + [np.zeros(n, dtype=np.int64)], 1), fmt="%d")
                for c in range(K):
                    np.savetxt(f, np.stack([-(1 + E[:, 0] * K + c), -(1 + E[:, 1] * K + c),
                                            np.zeros(len(E), dtype=np.int64)], 1), fmt="%d")
            r = subprocess.run(["timeout", str(ktime), KISSAT, fn], capture_output=True, text=True)
            os.remove(fn)
            st = [l for l in r.stdout.splitlines() if l.startswith("s ")]
            res = None if not st else ("UNSAT" not in st[0])
        else:
            print("  too large for pysat; pass --kissat SECONDS"); results.append(None); continue
        print(f"  4-colourable: {res}", flush=True)
        results.append(res)
    return results


if __name__ == "__main__":
    argv = sys.argv[1:]
    ktime = None
    if "--kissat" in argv:
        i = argv.index("--kissat"); ktime = int(argv[i + 1]); argv = argv[:i] + argv[i + 2:]
    d = json.load(open(argv[0]))
    rd = lambda xy: ([Fr(a, b) for a, b in xy[0]], [Fr(a, b) for a, b in xy[1]])
    U = [from_xy(*rd(u)) for u in d["units"]]
    p = int(argv[1])
    levels = [tuple(map(int, s.split(","))) for s in argv[2:]]
    gate(U, p, levels, ktime, tag=os.path.basename(argv[0]).replace(".json", ""))
