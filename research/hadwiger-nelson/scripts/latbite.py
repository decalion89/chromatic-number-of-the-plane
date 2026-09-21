"""Bite a lattice: de Grey's move applied to someone else's point set.

The two-distance lattice sweep failed for a clean reason -- a Cayley graph of
Z^2 is homogeneous, and homogeneity buys periodic colourings.  What Sa has and
a lattice has not is SCALE: Sa is closed under a rotation of infinite order,
so its points sit at radii 1, 1/3, 5/9, 4/3, 5/3, and no translation preserves
it.

So give the lattice that.  Take a triangular patch -- uniform degree six,
3-chromatic on its own, far denser locally than any orbit -- and close it under
rho_D, the rotation that makes the ring of squared radius D bite its own image.
The union is inhomogeneous, multi-scale, and is not de Grey's point set; it is
his MOVE on someone else's.

D must be Loeschian, so the lattice has a ring there, and closable, so the
spindle exists afterwards: D in {3, 4, 7, 9} all qualify over Q(v3,v5,v7,v11),
wanting sqrt(11), sqrt(15), sqrt(27), sqrt(35) respectively.

One SAT call per step asks the only thing that ends it outright: is the union
still 5-colourable.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.field import Field
from hn.geometry import Point, rotation_joining, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
half = K.rational(Fr(1, 2))
rt3 = K.sqrt(3) * K.rational(Fr(1, 2))


def patch(R):
    """Triangular lattice, spacing one, inside radius R."""
    out = []
    n = int(R) + 2
    for a in range(-n, n + 1):
        for b in range(-n, n + 1):
            if a * a + a * b + b * b <= R * R:
                out.append(Point(K.rational(a) + half * K.rational(b),
                                 rt3 * K.rational(b)))
    return out


def solve(P):
    g = build_graph(P)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return g.n, len(E), ok


O = Point(K.zero(), K.zero())
for D in (Fr(3), Fr(4), Fr(7), Fr(9)):
    for R in (5, 8):
        L = patch(R)
        ring = sum(1 for p in L if p.norm2() == D)
        b = rotation_joining(D, K)
        f, iv = b.about(O), Rotation(b.cos, -b.sin).about(O)
        print(f"\nD={D}, patch radius {R}: {len(L)} points, {ring} on the "
              f"ring; bite cos {b.cos} sin {b.sin}  [{time.time()-t0:.0f}s]",
              flush=True)
        seen, U = set(L), list(L)
        front = list(L)
        for m in range(1, 5):
            nxt = []
            for q in front:
                for r in (f, iv):
                    w = r(q)
                    nxt.append(w)
                    if w not in seen:
                        seen.add(w)
                        U.append(w)
            front = nxt[:len(L)]
            n, e, ok = solve(U)
            print(f"   m={m}: {n} pts, {e} edges, {e/n:.2f} per vertex -> "
                  f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            if not ok:
                print("*** chi >= 6 ***", flush=True)
                sys.exit(0)
            if n > 9000:
                break
print("\nDONE", flush=True)
