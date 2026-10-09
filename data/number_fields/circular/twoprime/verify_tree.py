"""Independent verifier of tree_r7_25.txt (completeness of the enumeration of S^r(1,1) at r = 7/25).

Claim checked: every c in C with Re(conj(c) g) in [r, 1-r] + Z for all g in G(1,1) = {u rho^j sigma^l}
has, after translation by 65 Z[i], the strip-index vector of one of the 13 components listed in
prop7_certificates.txt.

Logic: translate c so that its cell (strip indices of x = Re(conj c) and -y = Im(conj c)) is (a, -b-1) with
0 <= a, b <= 64.  Along the tree, for each functional f_t the strip index of c is one of the candidates (all n with
[n + r, n + 1 - r] meeting the range of f_t on the cell, recomputed here from the corners).  A candidate marked D
carries a Farkas certificate (lam >= 0, sum lam (a, b) = 0, sum lam c < 0 for valid constraints a x + b y <= c on that
branch), so no point lies on that branch; hence c follows S-branches to a leaf, whose index vector is checked to be
one of the 13.  The functionals are rebuilt here from Gaussian integers (no shared code)."""
from fractions import Fraction as Fr
import sys


def gdiv(a, b):
    n = b[0] * b[0] + b[1] * b[1]
    return (Fr(a[0] * b[0] + a[1] * b[1], n), Fr(a[1] * b[0] - a[0] * b[1], n))


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z = gdiv((1, 0), z)
        e = -e
    for _ in range(e):
        w = gmul(w, z)
    return w


RHO = gdiv((2, 1), (2, -1))
SIG = gdiv((3, 2), (3, -2))
FUN = []
for j in (-1, 0, 1):
    for l in (-1, 0, 1):
        g = gmul(gpow(RHO, j), gpow(SIG, l))
        FUN.append((g[0], g[1]))          # Re(conj(c) g) = g0 x + g1 y
        FUN.append((g[1], -g[0]))         # Im(conj(c) g) = g1 x - g0 y
assert FUN[8] == (1, 0) and FUN[9] == (0, -1)


def fl(q):
    return q.numerator // q.denominator


r = Fr(7, 25)
LO, HI = r, 1 - r
targets = set()
for ln in open("prop7_certificates.txt"):
    if ln.startswith("component"):
        targets.add(tuple(int(x) for x in ln.split("|")[0].split("indices")[1].split()))
assert len(targets) == 13

lines = open("tree_r7_25.txt").read().split("\n")
header = lines[0]
assert header.startswith("# r = 7/25")
ORDER = [int(x) for x in header.split("functional order:")[1].split()]
assert sorted(ORDER + [8, 9]) == list(range(18))
pos = 1
cells_seen = set()
nD = nS = nL = 0
leaves = set()


def constraint(ref, cell, idx):
    """the constraint a x + b y <= c named ref on the current branch"""
    a, b = cell
    name, side = ref.split()
    if name == 'x':
        return (Fr(1), Fr(0), a + HI) if side == 'hi' else (Fr(-1), Fr(0), -(a + LO))
    if name == 'y':
        return (Fr(0), Fr(1), b + HI) if side == 'hi' else (Fr(0), Fr(-1), -(b + LO))
    t = int(name)
    n = idx[t]                         # KeyError if the functional is not fixed on this branch
    fa, fb = FUN[t]
    return (fa, fb, n + HI) if side == 'hi' else (-fa, -fb, -(n + LO))


def check_node(depth, cell, idx, corners):
    global pos, nD, nS, nL
    if depth == len(ORDER):
        assert lines[pos] == "L", (pos, lines[pos])
        pos += 1
        iv = tuple(idx[t] for t in range(18))
        assert iv in targets, iv
        leaves.add(iv)
        nL += 1
        return
    t = ORDER[depth]
    fa, fb = FUN[t]
    vals = [fa * p[0] + fb * p[1] for p in corners]
    mn, mx = min(vals), max(vals)
    cands = [n for n in range(fl(mn) - 3, fl(mx) + 3) if n + HI >= mn and n + LO <= mx]
    for n in cands:
        ln = lines[pos]
        pos += 1
        kind, rest = ln[0], ln[2:]
        if kind == 'D':
            nstr, cert = rest.split(" ", 1)
            assert int(nstr) == n, (ln, n)
            idx2 = dict(idx)
            idx2[t] = n
            sa = sb = sc = Fr(0)
            for term in cert.split(" ; "):
                parts = term.split()
                ref, lam = parts[0] + " " + parts[1], Fr(parts[2])
                assert lam >= 0
                a, b, c = constraint(ref, cell, idx2)
                sa += lam * a; sb += lam * b; sc += lam * c
            assert sa == 0 and sb == 0 and sc < 0, ln
            nD += 1
        else:
            assert kind == 'S' and int(rest) == n, (ln, n)
            nS += 1
            idx2 = dict(idx)
            idx2[t] = n
            check_node(depth + 1, cell, idx2, corners)


while pos < len(lines) and lines[pos]:
    ln = lines[pos]
    assert ln.startswith("C "), ln
    a, b = (int(x) for x in ln.split()[1:])
    assert 0 <= a <= 64 and 0 <= b <= 64 and (a, b) not in cells_seen
    cells_seen.add((a, b))
    pos += 1
    corners = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
    check_node(0, (a, b), {8: a, 9: -b - 1}, corners)
assert len(cells_seen) == 65 * 65
assert leaves == targets
print(f"tree verified: {len(cells_seen)} cells, {nS} branch nodes, {nD} Farkas-closed branches, {nL} leaves = the 13 "
      f"index vectors of prop7_certificates.txt")
