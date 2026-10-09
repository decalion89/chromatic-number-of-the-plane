#!/usr/bin/env python3
"""Tiny sanity check of the remark 'Lemma F17 makes Corollary F14 effective' on distance graphs.

U = D, a set of positive integers (vectors of Z, so U and -U are disjoint), with kappa(D) <= 1/4 checked
exactly (no t with ||t d|| > 1/4 for all d in D; the characters of Z are x -> t x mod 1).
Steps: Z-basis r_j of Rel; Pi; for each p in Pi an integer y != 0 in the cone
C_p = {y : <y,p> >= max over the box [1/4,3/4]^U of <yR,.>} (checked on all box vertices), rho_p = yR;
infeasibility of the integer system on S; chains from 0 and bubble sorts in Cay(Z, +-D); H = vertices used;
G = induced subgraph of the distance graph on H.  Then, independently of Lemma P, a SAT check that G has
no proper 4-colouring with an acyclic tight digraph (levels encoding), i.e. chi_c(G) >= 4 (Guichard).
Control: the same for the graph built from the chains of the r_j alone.
Usage: python3 remark_check.py d1 d2 ... (default 2 3 5 6); kissat from the environment variable KISSAT, else PATH.
"Lemma P" is Lemma F17 of the note (Lemma 24 of the paper), "the remark" the remark after it.
"""
import atexit
import itertools
import os
import subprocess
import sys
import tempfile
from fractions import Fraction as Fr

KISSAT = os.environ.get("KISSAT", "kissat")
OUT = tempfile.mkdtemp(prefix="remark_check_")
atexit.register(os.rmdir, OUT)

D = [int(a) for a in sys.argv[1:]] or [2, 3, 5, 6]
m = len(D)


def kappa_gt_quarter(D):
    I = [(Fr(0), Fr(1))]
    for d in D:
        J = [(Fr(4 * k + 1, 4 * d), Fr(4 * k + 3, 4 * d)) for k in range(d)]
        I = [(max(a0, b0), min(a1, b1)) for a0, a1 in I for b0, b1 in J if max(a0, b0) < min(a1, b1)]
    return bool(I)


def kernel_basis(row):
    """Z-basis of {n in Z^m : sum row_i n_i = 0} by unimodular column operations."""
    m = len(row)
    V = [[int(i == j) for j in range(m)] for i in range(m)]   # columns of V
    a = list(row)
    while sum(1 for x in a if x) > 1:
        i = min((abs(x), i) for i, x in enumerate(a) if x)[1]
        for j in range(m):
            if j != i and a[j]:
                q = a[j] // a[i]
                a[j] -= q * a[i]
                for t in range(m):
                    V[t][j] -= q * V[t][i]
    return [tuple(V[t][j] for t in range(m)) for j in range(m) if a[j] == 0]


def PN(rho):
    return sum(x for x in rho if x > 0), -sum(x for x in rho if x < 0)


def in_range(v, rho):
    P, N = PN(rho)
    return Fr(P - 3 * N, 4) < v < Fr(3 * P - N, 4)


def box_max(w):
    return sum(max(Fr(x, 4), Fr(3 * x, 4)) for x in w)


assert not kappa_gt_quarter(D), "kappa(D) > 1/4"
R = kernel_basis(D)
k = len(R)
assert all(sum(r[i] * D[i] for i in range(m)) == 0 for r in R)
print("D =", D, " kappa(D) <= 1/4 (exact); Z-basis of Rel:", R)

ranges = []
for r in R:
    P, N = PN(r)
    lo, hi = Fr(P - 3 * N, 4), Fr(3 * P - N, 4)
    ranges.append([z for z in range(int(lo) - 2, int(hi) + 3) if lo < z < hi])
Pi = list(itertools.product(*ranges))
print("|Pi| =", len(Pi))

rho_p = {}
vertices = list(itertools.product([Fr(1, 4), Fr(3, 4)], repeat=m))
for p in Pi:
    found = None
    for B in range(1, 16):
        for y in itertools.product(range(-B, B + 1), repeat=k):
            if max(abs(c) for c in y) != B:
                continue
            w = tuple(sum(y[j] * R[j][i] for j in range(k)) for i in range(m))
            if sum(yj * pj for yj, pj in zip(y, p)) >= box_max(w):
                found = (y, w)
                break
        if found:
            break
    assert found, ("no rho_p found for", p)
    y, w = found
    # the cone description: <y, p - R b> >= 0 for every vertex b of the box
    assert all(sum(y[j] * (p[j] - sum(R[j][i] * b[i] for i in range(m))) for j in range(k)) >= 0
               for b in vertices)
    v = sum(yj * pj for yj, pj in zip(y, p))
    assert any(w) and not in_range(v, w)
    rho_p[p] = (y, w)
print("rho_p found for every p in Pi; max |y_j| =", max(max(abs(c) for c in y) for y, w in rho_p.values()))

S = [(tuple(int(j == i) for j in range(k)), R[i]) for i in range(k)] + list(rho_p.values())
# (iii): no p in Z^k satisfies all ranges on S (only p in Pi can satisfy the r_j ranges; checked on a wider box)
wide = list(itertools.product(*[range(min(rg) - 3, max(rg) + 4) if rg else range(-3, 4) for rg in ranges]))
bad = [p for p in wide if all(in_range(sum(a * b for a, b in zip(coef, p)), rho) for coef, rho in S)]
print("(iii) p in a box of %d candidates satisfying every range on S: %d" % (len(wide), len(bad)))
assert not bad


def walk_steps(rho, sign=1):
    st = []
    for i, n in enumerate(rho):
        st += [(i, 1 if n > 0 else -1)] * abs(n)
    if sign < 0:
        st = [(i, -s) for (i, s) in reversed(st)]
    return st


def pts_of(st):
    x = 0
    out = [0]
    for i, s in st:
        x += s * D[i]
        out.append(x)
    return out


def build_H(S):
    H = {0}
    for coef, rho in S:
        cur = [0] * k
        seq = []
        for j, a in enumerate(coef):
            seq += [(j, 1 if a > 0 else -1)] * abs(a)
        for j, s in seq:
            r0 = tuple(sum(cur[i] * R[i][u] for i in range(k)) for u in range(m))
            cur[j] += s
            r1 = tuple(sum(cur[i] * R[i][u] for i in range(k)) for u in range(m))
            st = walk_steps(r0) + walk_steps(R[j], s)
            pts = pts_of(st)
            H.update(pts)
            changed = True
            while changed:
                changed = False
                a = 0
                while a < len(st) - 1:
                    x, z = st[a], st[a + 1]
                    if x[0] == z[0] and x[1] == -z[1]:
                        del st[a:a + 2]
                        del pts[a + 1:a + 3]
                        changed = True
                        a = max(a - 1, 0)
                        continue
                    if x[0] > z[0]:
                        new = pts[a] + z[1] * D[z[0]]
                        H.add(new)
                        st[a], st[a + 1] = z, x
                        pts[a + 1] = new
                        changed = True
                    a += 1
            assert st == walk_steps(r1)
            H.update(pts_of(walk_steps(r1)))
        assert tuple(cur) == tuple(coef) and tuple(sum(cur[i] * R[i][u] for i in range(k))
                                                   for u in range(m)) == tuple(rho)
    return sorted(H)


def no_tight_cycle_colourable(H, tag):
    """SAT iff G[H] has a proper 4-colouring whose tight digraph is acyclic (rank function encoding)."""
    n = len(H)
    idx = {x: i for i, x in enumerate(H)}
    E = [(idx[x], idx[x + d]) for x in H for d in D if x + d in idx]
    X = lambda v, c: c * n + v + 1
    arcs = E + [(b, a) for (a, b) in E]
    T = {a: 4 * n + 1 + i for i, a in enumerate(arcs)}
    base = 4 * n + len(arcs)
    O = lambda v, i: base + v * (n - 1) + i   # O(v,i) = [level(v) >= i], i = 1..n-1
    cls = []
    for v in range(n):
        cls.append([X(v, c) for c in range(4)])
        for c in range(4):
            for c2 in range(c + 1, 4):
                cls.append([-X(v, c), -X(v, c2)])
        for i in range(1, n - 1):
            cls.append([-O(v, i + 1), O(v, i)])
    for (a, b) in E:
        for c in range(4):
            cls.append([-X(a, c), -X(b, c)])
    for (a, b) in arcs:
        t = T[(a, b)]
        for c in range(4):
            cls.append([-X(a, c), -X(b, (c + 1) % 4), t])
        # t(a,b) -> level(b) > level(a):  level(a) >= i  ->  level(b) >= i+1, i = 0..n-1
        for i in range(0, n):
            lit_a = [] if i == 0 else [-O(a, i)]
            if i + 1 <= n - 1:
                cls.append([-t] + lit_a + [O(b, i + 1)])
            else:
                cls.append([-t] + lit_a)
    nv = base + n * (n - 1)
    path = os.path.join(OUT, "remark_%s.cnf" % tag)
    with open(path, "w") as fh:
        fh.write("p cnf %d %d\n" % (nv, len(cls)))
        for c in cls:
            fh.write(" ".join(map(str, c)) + " 0\n")
    r = subprocess.run([KISSAT, path], capture_output=True, text=True)
    os.remove(path)
    s = [l for l in r.stdout.splitlines() if l.startswith("s ")]
    print("%s: |H| = %d, |E| = %d, CNF %d vars %d clauses -> %s" % (tag, n, len(E), nv, len(cls), s))
    return s


H = build_H(S)
s_full = no_tight_cycle_colourable(H, "G_S")
H0 = build_H(S[:k])
s_ctrl = no_tight_cycle_colourable(H0, "G_basis_only")
ok = s_full == ["s UNSATISFIABLE"]
print("RESULT:", "every proper 4-colouring of G has a tight cycle (chi_c(G) >= 4)" if ok else "UNEXPECTED")


def circular_colourable(H, p, q, tag):
    """SAT iff G[H] has a (p,q)-colouring: c: V -> Z_p with q <= |c(x)-c(y)| <= p-q on edges."""
    n = len(H)
    idx = {x: i for i, x in enumerate(H)}
    E = [(idx[x], idx[x + d]) for x in H for d in D if x + d in idx]
    X = lambda v, c: c * n + v + 1
    cls = [[X(v, c) for c in range(p)] for v in range(n)]
    for (a, b) in E:
        for c in range(p):
            for c2 in range(p):
                dd = (c - c2) % p
                if min(dd, p - dd) < q:
                    cls.append([-X(a, c), -X(b, c2)])
    path = os.path.join(OUT, "remark_circ_%s.cnf" % tag)
    with open(path, "w") as fh:
        fh.write("p cnf %d %d\n" % (p * n, len(cls)))
        for c in cls:
            fh.write(" ".join(map(str, c)) + " 0\n")
    r = subprocess.run([KISSAT, path], capture_output=True, text=True)
    os.remove(path)
    return [l for l in r.stdout.splitlines() if l.startswith("s ")]


Q = (len(H) + 1) // 4
if len(H) <= 64:
  print("independent check: (%d,%d)-colourability of G (largest fraction < 4 with numerator <= |H|):"
      % (4 * Q - 1, Q), circular_colourable(H, 4 * Q - 1, Q, "G_S"),
      "; 4-colourable:", circular_colourable(H, 4, 1, "G_S4"))
else:
  print("independent (4Q-1,Q) check skipped (|H| > 64); 4-colourable:", circular_colourable(H, 4, 1, "G_S4"))
