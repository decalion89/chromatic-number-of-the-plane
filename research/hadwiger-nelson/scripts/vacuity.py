"""Is Sa u rot(Sa) still 4-chromatic?  If not, the rho drop is vacuous.

rho is only meaningful when proper k-colourings exist.  If a union stops being
k-colourable there are none, every set is "forcing" for want of a
counterexample, and the measurement is empty -- which is exactly what the
triangular patch just did at three colours, returning rho = 1, below the floor
rho >= k.

So the same question has to be asked of the result it threatens: the drop from
7 to 5 on Sa u rot(Sa) is only real if that union is still 4-colourable.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))


def colourable(g, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


cur = list(degrey.build_Sa())
for copies in range(1, 8):
    if copies > 1:
        seen, grown = set(cur), list(cur)
        for p in cur:
            q = ROT(p)
            if q not in seen:
                seen.add(q)
                grown.append(q)
        cur = grown
    g = build_graph(cur)
    t = time.time()
    ok = colourable(g, 4)
    print(f"  {copies} copies: n={g.n:5} m={g.m:6}  4-colourable: {ok}"
          + ("" if ok else "   *** rho measurement was VACUOUS ***")
          + f"  [{time.time()-t:.0f}s]", flush=True)
    if not ok:
        break
