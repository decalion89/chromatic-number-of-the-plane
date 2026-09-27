"""Rank three: keep the Eisenstein density and block anyway.

A discrete module never blocks -- rank at most 2, three unit directions
against PG(1,5)'s six points.  But rank THREE can: PG(2,5) has 31 points and
its hyperplanes are lines of 6, and six CONCURRENT lines cover all 31.  So the
cheapest possible blocked module is Z[zeta_6] plus one more modulus-one step,
and that keeps the triangular lattice -- triangles, density, every translate
landing on a point -- inside it.

The graph is then layers: a patch of the Eisenstein lattice, the same patch
shifted by u, by 2u, and so on.  Within a layer it is the triangular lattice
and as dense as a unit-distance graph gets; between layers the edges are the
u-steps and any other unit differences that happen to close.

This looks for a u making {1, zeta_6, u} block, then builds and colours it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

F = CycloField(21)
D = F.degree
one = F.rational(1)
z3 = F.zeta(7)
z6 = F.neg(F.mul(z3, z3))
z7 = F.zeta(3)


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def prim(vs):
    den = 1
    for v in vs:
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
    out = set()
    for v in vs:
        w = tuple(int(x * den) for x in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w) if g > 1 else w)
    return sorted(out)


# candidate extra steps: modulus-one elements of Q(zeta_7) inside Q(zeta_21)
cands = []
for coeffs in itertools.product(range(-2, 3), repeat=4):
    if not any(coeffs):
        continue
    a, p = F.zero(), one
    for c in coeffs:
        a = F.add(a, tuple(Fraction(c) * x for x in p))
        p = F.mul(p, z7)
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) == one:
        cands.append(u)
print(f"{len(cands)} candidate extra steps", flush=True)

# the module <1, zeta_6, u>: its unit vectors are the six roots of unity,
# plus u and its Eisenstein multiples, plus whatever else closes up
eis, z = [], one
for _ in range(6):
    z = F.mul(z, z6)
    eis.append(z)

found = None
for u in cands:
    steps = set(eis) | {F.mul(e, u) for e in eis}
    steps |= {F.neg(s) for s in steps}
    phi, why = has_homomorphism(prim(steps), 5)
    if phi is None:
        found = (u, steps)
        print(f"  BLOCKS with {len(steps)} unit steps", flush=True)
        break
if found is None:
    print("  no single extra step blocks; the Eisenstein directions are "
          "three of six and one more orbit does not reach a pencil")
    # try two extra steps
    for u, v in itertools.combinations(cands[:40], 2):
        steps = set(eis)
        for w in (u, v):
            steps |= {F.mul(e, w) for e in eis}
        steps |= {F.neg(s) for s in steps}
        if has_homomorphism(prim(steps), 5)[0] is None:
            found = ((u, v), steps)
            print(f"  BLOCKS with two extra steps, {len(steps)} unit steps",
                  flush=True)
            break

if found is None:
    print("  and no pair of extra steps blocks either")
    sys.exit()

u, steps = found
steps = sorted(steps)
pts = {F.zero()}
for _ in range(3):
    pts = {F.add(p, s) for p in pts for s in steps} | pts
    print(f"    {len(pts)} points", flush=True)
    if len(pts) > 20000:
        break

P = sorted(pts)
idx = {p: i for i, p in enumerate(P)}
edges = sorted({(min(i, j), max(i, j))
                for i, p in enumerate(P) for s in steps
                for j in [idx.get(F.add(p, s))] if j is not None and j != i})
adj = [set() for _ in P]
for a, b in edges:
    adj[a].add(b)
    adj[b].add(a)
tri = sum(len(adj[a] & adj[b]) for a, b in edges) // 3
print(f"  graph: {len(P)} vertices, {len(edges)} edges, {tri} triangles, "
      f"average degree {2*len(edges)/len(P):.2f}", flush=True)

# the densest part: the k-core, which is where any high chromatic number lives
for kc in (6, 8, 10, 12):
    alive = set(range(len(P)))
    deg = {v: len(adj[v]) for v in alive}
    ch = True
    while ch:
        ch = False
        for v in list(alive):
            if deg[v] < kc:
                alive.discard(v)
                for w in adj[v] & alive:
                    deg[w] -= 1
                ch = True
    print(f"  {kc}-core: {len(alive)} vertices", flush=True)
    if not alive:
        break

for k in (3, 4, 5):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"  {k}-colourable: {r}" + ("" if r else f"   *** chi > {k} ***"),
          flush=True)
    if r:
        break
