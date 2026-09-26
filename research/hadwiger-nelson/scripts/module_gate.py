"""Split places matter for modules: reduce a graph's module at a place and colour the finite image.

Take the points of a unit-distance graph and its growth unit set. Suppose all their coordinates are
integral at a place w of F = Q(sqrt d_1, ..., sqrt d_n) above an odd prime p. Here "integral" means
p-adic valuation >= 0 in the completion, which is checked exactly with p-adic arithmetic; a coordinate
may have p in its denominator and still be integral at w.
- Then reduction mod w is a group homomorphism from the module M they generate into F_{p^2}^2.
- Every unit vector goes to the circle a^2 + b^2 = 1, since the reduction is a ring map on
  coordinates.
- So every unit-distance graph with vertices in M maps homomorphically to Cay(A, A ∩ circle), where
  A is the image of M. This covers every unit edge that could ever appear in M, not only the
  present ones.
- If that finite Cayley graph is k-colourable, no growth inside M can pass k colours.

At the level of the whole field only non-split places give such colourings
(notes/local_colourings.md). For a finitely generated module, split places do too, until a unit that
is not integral at w is added.

Primes unramified in F are examined directly. At a prime that ramifies in exactly one generator,
ramified_place() writes each coordinate as A + B sqrt(d_r), with A, B in the unramified part, and
reduces that; any other ramified prime is reported and skipped.

usage: module_gate.py graph.json [pmax] [k] [maxsize]
"""
import sys, json, os, time, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from sympy import isprime
from pysat.solvers import Solver

NPREC = 30


def load(path):
    d = json.load(open(path)); gens = tuple(d["field_generators"])
    rd = lambda xy: ([Fr(a, b) for a, b in xy[0]], [Fr(a, b) for a, b in xy[1]])
    pts = [rd(xy) for xy in d["points"]]
    units = [rd(xy) for xy in d.get("units", [])]
    return gens, pts, units


class Place:
    """An embedding of F = Q(sqrt d_j) into Q_{p^2} = Q_p(tau), tau^2 = nr, to precision p^NPREC."""

    def __init__(self, p, gens, signs):
        self.p = p; self.P = p ** NPREC
        self.nr = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        n = len(gens); roots = []
        for dj, sg in zip(gens, signs):
            if pow(dj % p, (p - 1) // 2, p) == 1:
                r = self._sqrt(dj); roots.append(((sg * r) % self.P, 0))
            else:
                r = self._sqrt(dj * pow(self.nr, -1, self.P)); roots.append((0, (sg * r) % self.P))
        self.bimg = []
        for m in range(1 << n):
            v = (1, 0)
            for j in range(n):
                if m >> j & 1: v = self.mul(v, roots[j])
            self.bimg.append(v)

    def _sqrt(self, a):
        p, P = self.p, self.P
        a %= P
        r = next(t for t in range(1, p) if (t * t - a) % p == 0)
        for _ in range(8):                       # Newton: doubles the precision each step
            r = (r - (r * r - a) * pow(2 * r, -1, P)) % P
        assert (r * r - a) % P == 0
        return r

    def mul(self, x, y):
        P = self.P
        return ((x[0] * y[0] + self.nr * x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)

    def image(self, coeffs):
        """(valuation, residue in F_p^2) of sum_m c_m sqrt(d_m)."""
        p, P = self.p, self.P
        ks = []
        for c in coeffs:
            k, den = 0, c.denominator
            while den % p == 0: den //= p; k += 1
            ks.append(k)
        K = max(ks) if ks else 0
        a0 = a1 = 0
        for c, k, (u0, u1) in zip(coeffs, ks, self.bimg):
            if c == 0: continue
            den = c.denominator // p ** k
            t = c.numerator * p ** (K - k) * pow(den, -1, P)
            a0 += t * u0; a1 += t * u1
        a0 %= P; a1 %= P
        if a0 == 0 and a1 == 0: return NPREC, (0, 0)
        v = 0
        while a0 % p == 0 and a1 % p == 0: a0 //= p; a1 //= p; v += 1
        val = v - K
        if val < 0: return val, None
        if val > 0: return val, (0, 0)
        return 0, (a0 % p, a1 % p)


def span(vecs, p):
    basis = []
    for v in vecs:
        v = [x % p for x in v]
        for piv, b in basis:
            if v[piv]:
                f = v[piv]; v = [(a - f * c) % p for a, c in zip(v, b)]
        nz = [i for i, x in enumerate(v) if x]
        if nz:
            piv = nz[0]; inv = pow(v[piv], -1, p); v = [x * inv % p for x in v]
            basis = [(pv, [(a - bb[piv] * c) % p for a, c in zip(bb, v)]) for pv, bb in basis]
            basis.append((piv, v))
    return basis


def tabu_colour(A, S, p, k, iters=30000000, seed=5):
    import subprocess
    idx = {v: i for i, v in enumerate(A)}
    E = set()
    for v in A:
        for s_ in S:
            w = tuple((a + b) % p for a, b in zip(v, s_))
            i, j = idx[v], idx[w]
            if i < j: E.add((i, j))
    TC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tabucol")
    inp = f"{len(A)} {len(E)} {k} {iters} {seed}\n" + "".join(f"{i} {j}\n" for i, j in E) + " ".join(["-1"] * len(A)) + "\n"
    out = subprocess.run([TC], input=inp, capture_output=True, text=True).stdout.split()
    return out[0] == "OK", (0 if out[0] == "OK" else int(out[1]))


def colourable(A, S, p, k):
    idx = {v: i for i, v in enumerate(A)}
    E = set()
    for v in A:
        for s in S:
            w = tuple((a + b) % p for a, b in zip(v, s))
            i, j = idx[v], idx[w]
            if i < j: E.add((i, j))
    X = lambda v, c: 1 + v * k + c
    cnf = [[X(v, c) for c in range(k)] for v in range(len(A))] + [[-X(i, c), -X(j, c)] for i, j in E for c in range(k)]
    cnf.append([X(0, 0)])
    s = Solver(name="cd19", bootstrap_with=cnf); r = s.solve(); s.delete()
    return r, len(E)


def ramified_place(p, gens, elems, KMAX, MAXSIZE):
    """p divides exactly one generator d_r (all fields used here). Write each coordinate as
    A + B sqrt(d_r) with A, B in the unramified part; sqrt(d_r) has valuation 1/2, so the coordinate
    is integral iff A and B are, and its residue is A's."""
    rs = [j for j, dj in enumerate(gens) if dj % p == 0]
    if len(rs) != 1 or (gens[rs[0]] // p) % p == 0:
        print(f"  p={p}: ramified in a way not handled", flush=True); return
    r = rs[0]; n = len(gens)
    unr = [j for j in range(n) if j != r]
    sub = tuple(gens[j] for j in unr)

    def split(coeffs):
        A = [Fr(0)] * (1 << len(unr)); B = [Fr(0)] * (1 << len(unr))
        for m, c in enumerate(coeffs):
            if c == 0: continue
            mm = sum(1 << t for t, j in enumerate(unr) if m >> j & 1)
            (B if m >> r & 1 else A)[mm] += c
        return A, B

    seen = set(); killed = 0; tot = 0
    for signs in itertools.product((1, -1), repeat=len(unr)):
        tot += 1
        pl = Place(p, sub, signs)
        imgs = []; ok = True
        for x, y in elems:
            row = []
            for co in (x, y):
                A, B = split(list(co))
                va, ra = pl.image(A); vb, _ = pl.image(B)
                if va < 0 or vb < 0: ok = False; break
                row += list(ra)
            if not ok: break
            imgs.append(tuple(row))
        if not ok: killed += 1; continue
        B_ = span(imgs, p)
        key = tuple(sorted(tuple(b) for _, b in B_))
        if key in seen: continue
        seen.add(key)
        report(p, f"ramified, place {signs}", B_, pl.nr, KMAX, MAXSIZE)
    if killed == tot:
        print(f"  p={p}: ramified; every place has a non-integral generator (killed)", flush=True)


def report(p, label, B, nr, KMAX, MAXSIZE):
    size = p ** len(B)
    if size > MAXSIZE:
        print(f"  p={p} {label}: integral; image has {size} points (> {MAXSIZE}), not tested", flush=True)
        return
    A = []
    for co in itertools.product(range(p), repeat=len(B)):
        v = [0, 0, 0, 0]
        for c, (_, b) in zip(co, B): v = [(a + c * bb) % p for a, bb in zip(v, b)]
        A.append(tuple(v))
    sq = lambda a, b: ((a * a + nr * b * b) % p, (2 * a * b) % p)
    S = [v for v in A if tuple((s + t) % p for s, t in zip(sq(v[0], v[1]), sq(v[2], v[3]))) == (1, 0)]
    t0 = time.time(); best = None
    if len(A) <= 400:
        for k in range(2, KMAX + 1):
            r, ne = colourable(A, S, p, k)
            if r: best = k; break
        verdict = f"chi = {best}" if best else f"chi > {KMAX} (SAT)"
    else:
        ok, conf = tabu_colour(A, S, p, KMAX)
        best = KMAX if ok else None
        verdict = f"{KMAX}-colourable (tabu)" if ok else f"tabu found no {KMAX}-colouring (best {conf} conflicts)"
    print(f"  p={p} {label}: integral; image {size} points, {len(S)} unit residues, {verdict}  [{time.time()-t0:.0f}s]", flush=True)
    if best:
        print(f"  ==> every graph in this module is {best}-colourable (place above {p})", flush=True)


def generators(path):
    """Growth units, the distinct unit vectors along edges, and one point per component minus the first
    point: the module every growth step stays in is spanned by these."""
    import numpy as np
    from scipy.spatial import cKDTree
    from hn.field import Field
    from hn.geometry import Point
    gens, pts, units = load(path)
    F = Field(gens)
    P = [Point(F.element(list(x)), F.element(list(y))) for x, y in pts]
    xy = np.array([[q.fx, q.fy] for q in P]); tree = cKDTree(xy)
    pr = tree.query_pairs(1 + 1e-6, output_type="ndarray")
    dd = ((xy[pr[:, 0]] - xy[pr[:, 1]]) ** 2).sum(1); pr = pr[np.abs(dd - 1) < 1e-6]
    one = F.one(); U = {}; parent = list(range(len(P)))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for i, j in pr:
        u = P[j] - P[i]
        if u.x * u.x + u.y * u.y == one:
            U[(tuple(u.x.c), tuple(u.y.c))] = (tuple(u.x.c), tuple(u.y.c))
            parent[find(i)] = find(j)
    reps = {}
    for i in range(len(P)): reps.setdefault(find(i), i)
    comps = [(tuple(a - b for a, b in zip(pts[i][0], pts[0][0])), tuple(a - b for a, b in zip(pts[i][1], pts[0][1])))
             for i in reps.values() if i != 0]
    elems = [(tuple(x), tuple(y)) for x, y in units] + list(U.values()) + comps
    return gens, elems, len(pts), len(units), len(U), len(reps)


def main():
    path = sys.argv[1]
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    MAXSIZE = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    PMIN = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    gens, elems, npts, nun, nedge, ncomp = generators(path)
    print(f"{path}: {npts} points, {nun} growth units, {nedge} edge unit vectors, {ncomp} components; field sqrt{gens}", flush=True)
    for p in range(PMIN, PMAX + 1):
        if not isprime(p): continue
        if any(dj % p == 0 for dj in gens):
            ramified_place(p, gens, elems, KMAX, MAXSIZE); continue
        seen = set(); killed = 0
        for signs in itertools.product((1, -1), repeat=len(gens)):
            pl = Place(p, gens, signs)
            imgs = []; ok = True
            for x, y in elems:
                vx, rx = pl.image(list(x)); vy, ry = pl.image(list(y))
                if vx < 0 or vy < 0: ok = False; break
                imgs.append(rx + ry)
            if not ok: killed += 1; continue
            B = span(imgs, p)
            key = tuple(sorted(tuple(b) for _, b in B))
            if key in seen: continue
            seen.add(key)
            report(p, f"place {signs}", B, pl.nr, KMAX, MAXSIZE)
        if killed == 1 << len(gens):
            print(f"  p={p}: every place has a non-integral generator (killed)", flush=True)


if __name__ == "__main__":
    main()
