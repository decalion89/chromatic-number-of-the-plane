"""Independent verifier for cert_K{K}M{M}.txt and tree_K{K}M{M}.txt (no code shared with certify_window.py,
probe2d.py or lpexact.py).

Checks, with the rotations rho^j sigma^l rebuilt from Gaussian integers ((2+i)/(2-i), (3+2i)/(3-2i)):
 (1) tree: every point c with Re(conj(c) g) in [r0, 1-r0] + Z for all g in the window has, after translation by
     N Z[i], one of the listed index vectors (branch tree with exact Farkas certificates, candidates recomputed
     from the cell corners);
 (2) MAIN lines: the type point N h or (N/3)(a+bi), translated by N m, has exactly the listed indices;
 (3) SEVEN lines: the listed point P satisfies 7P/N in Z[i] (a 7-torsion point), has the listed indices, and all
     its margins are >= 2/7 (so its values on the window, which contains rho mod 7 of order 8 times units, i.e. all
     of mu_8 mod 7 when K >= 1, lie in {2,..,5}/7);
 (4) EXTRA lines: lam >= 0, sum lam = 1, the linear parts cancel and sum lam*const equals the stated kappa; prints the
     largest such kappa.
usage: python3 verify_window.py K M"""
from fractions import Fraction as Fr
import sys

K, M = int(sys.argv[1]), int(sys.argv[2])
N = 5 ** K * 13 ** M


def gdiv(a, b):
    n = b[0] * b[0] + b[1] * b[1]
    return (Fr(a[0] * b[0] + a[1] * b[1], n), Fr(a[1] * b[0] - a[0] * b[1], n))


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z, e = gdiv((1, 0), z), -e
    for _ in range(e):
        w = gmul(w, z)
    return w


RHO, SIG = gdiv((2, 1), (2, -1)), gdiv((3, 2), (3, -2))


def functional(j, l, e):
    g = gmul(gpow(RHO, j), gpow(SIG, l))
    return (g[0], g[1]) if e == 'Re' else (g[1], -g[0])


def fl(q):
    return q.numerator // q.denominator


cert = open(f"cert_K{K}M{M}.txt").read().splitlines()
h0 = cert[0]
r0 = Fr(h0.split("r0=")[1].split(";")[0])
LO, HI = r0, 1 - r0
keys = [tuple(x.strip("()").split(",")) for x in cert[1].split(":", 1)[1].split()]
keys = [(int(j), int(l), e) for (j, l, e) in keys]
assert sorted((j, l) for (j, l, e) in keys if e == 'Re') == sorted((j, l) for j in range(-K, K + 1) for l in range(-M, M + 1))
FUN = [functional(*k) for k in keys]
t_re, t_im = keys.index((0, 0, 'Re')), keys.index((0, 0, 'Im'))
assert FUN[t_re] == (1, 0) and FUN[t_im] == (0, -1)

targets = set()
counts = {}
maxkappa = None
for ln in cert[2:]:
    if not ln.startswith("component"):
        continue
    head, lab = ln.split(" | ", 1)
    iv = tuple(int(x) for x in head.split("indices")[1].split())
    assert len(iv) == len(FUN) and iv not in targets
    targets.add(iv)
    parts = lab.split()

    def same_indices(P):
        for f, n in zip(FUN, iv):
            v = f[0] * P[0] + f[1] * P[1]
            if not (n + LO <= v <= n + HI):
                return False
        return True

    if parts[0] == 'MAIN':
        name, m1, m2 = parts[1], int(parts[2]), int(parts[3])
        T = (Fr(N, 2), Fr(N, 2)) if name == 'C' else (Fr(N * int(name[1]), 3), Fr(N * int(name[2]), 3))
        assert same_indices((T[0] + N * m1, T[1] + N * m2)), ln
        counts[name[0]] = counts.get(name[0], 0) + 1
    elif parts[0] == 'SEVEN':
        P = (Fr(parts[1]), Fr(parts[2]))
        assert (7 * P[0] / N).denominator == 1 and (7 * P[1] / N).denominator == 1
        assert same_indices(P), ln
        mn = min(min(f[0] * P[0] + f[1] * P[1] - n, n + 1 - (f[0] * P[0] + f[1] * P[1])) for f, n in zip(FUN, iv))
        assert mn == Fr(2, 7), mn
        counts['7'] = counts.get('7', 0) + 1
    else:
        assert parts[0] == 'EXTRA' and parts[1] == 'kappa' and parts[2] == '='
        kap = Fr(parts[3])
        terms = lab.split(" ; cert ")[1].split(" ; ")
        sa = sb = sc = tot = Fr(0)
        for term in terms:
            j, l, e, sg, lam = term.split()
            key = (int(j), int(l), e)
            lam = Fr(lam)
            assert lam >= 0
            f = FUN[keys.index(key)]
            n = iv[keys.index(key)]
            if sg == '+':
                sa += lam * f[0]; sb += lam * f[1]; sc += lam * (-n)
            else:
                sa -= lam * f[0]; sb -= lam * f[1]; sc += lam * (n + 1)
            tot += lam
        assert tot == 1 and sa == 0 and sb == 0 and sc == kap, ln
        counts['extra'] = counts.get('extra', 0) + 1
        maxkappa = kap if maxkappa is None or kap > maxkappa else maxkappa

# ---------------- the tree
lines = open(f"tree_K{K}M{M}.txt").read().split("\n")
assert lines[0].startswith(f"# window K={K} M={M} N={N} r0={r0}")
ORDER = [int(x) for x in lines[0].split("functional order:")[1].split()]
assert sorted(ORDER + [t_re, t_im]) == list(range(len(FUN)))
pos = 1
seen = set()
leaves = set()
nS = nD = 0


def constraint(ref, cell, idx):
    a, b = cell
    name, side = ref.split()
    if name == 'x':
        return (Fr(1), Fr(0), a + HI) if side == 'hi' else (Fr(-1), Fr(0), -(a + LO))
    if name == 'y':
        return (Fr(0), Fr(1), b + HI) if side == 'hi' else (Fr(0), Fr(-1), -(b + LO))
    t = int(name)
    n = idx[t]
    fa, fb = FUN[t]
    return (fa, fb, n + HI) if side == 'hi' else (-fa, -fb, -(n + LO))


def node(depth, cell, idx, corners):
    global pos, nS, nD
    if depth == len(ORDER):
        assert lines[pos] == "L"
        pos += 1
        iv = tuple(idx[t] for t in range(len(FUN)))
        assert iv in targets
        leaves.add(iv)
        return
    t = ORDER[depth]
    fa, fb = FUN[t]
    vals = [fa * p[0] + fb * p[1] for p in corners]
    mn, mx = min(vals), max(vals)
    for n in [n for n in range(fl(mn) - 3, fl(mx) + 3) if n + HI >= mn and n + LO <= mx]:
        ln = lines[pos]
        pos += 1
        idx2 = dict(idx)
        idx2[t] = n
        if ln[0] == 'D':
            nstr, rest = ln[2:].split(" ", 1)
            assert int(nstr) == n
            sa = sb = sc = Fr(0)
            for term in rest.split(" ; "):
                p1, p2, lam = term.split()
                lam = Fr(lam)
                assert lam >= 0
                a, b, c = constraint(p1 + " " + p2, cell, idx2)
                sa += lam * a; sb += lam * b; sc += lam * c
            assert sa == 0 and sb == 0 and sc < 0, ln
            nD += 1
        else:
            assert ln == f"S {n}", (ln, n)
            nS += 1
            node(depth + 1, cell, idx2, corners)


while pos < len(lines) and lines[pos]:
    a, b = (int(x) for x in lines[pos].split()[1:])
    assert lines[pos].startswith("C ") and 0 <= a < N and 0 <= b < N and (a, b) not in seen
    seen.add((a, b))
    pos += 1
    corners = [(a + LO, b + LO), (a + HI, b + LO), (a + HI, b + HI), (a + LO, b + HI)]
    node(0, (a, b), {t_re: a, t_im: -b - 1}, corners)
assert len(seen) == N * N and leaves == targets
print(f"window (K,M)=({K},{M}), N={N}, r0={r0}: tree verified ({len(seen)} cells, {nS} branch nodes, {nD} Farkas-closed "
      f"branches, {len(leaves)} leaves); labels {counts}; largest kappa of the EXTRA components = {maxkappa}")
