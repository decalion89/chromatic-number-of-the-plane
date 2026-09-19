"""Check the rho framework where the method DOES work: four colours.

A core of size r at p means N(p) u T uses all k colours in every colouring,
so rho <= deg(p) + r.  At five colours de Grey's G would need rho <= 63 and
has rho = 1581.  The framework is only worth that much if it comes out right
at four, where the construction actually succeeds: de Grey's Sa reaches
pressure 3 with cores down to 7 at a pivot of degree 30, which predicts

    rho(Sa, 4) <= 37.

That is a prediction with a number in it, and it is either met or the whole
reading of rho is wrong.  Measured here alongside the Moser spindle, whose
seven vertices at four colours are the smallest case there is.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.graph import build_graph
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.forced import ColourRelations, forcing_set, shrink_forcing_set

FLD = Field((3, 11))


def moser_spindle():
    h = FLD.rational(Fraction(1, 2))
    u = Point(h, FLD.sqrt(3) * h)
    one = Point(FLD.one(), FLD.zero())
    rh = [Point(FLD.zero(), FLD.zero()), one, u, one + u]
    rho = Rotation(FLD.rational(Fraction(5, 6)),
                   FLD.sqrt(11) * FLD.rational(Fraction(1, 6)))
    return rh + [rho(p) for p in rh[1:]]


for name, pts, k in (("Moser spindle", moser_spindle(), 4),
                     ("de Grey Sa", degrey.build_Sa(), 4),
                     ("de Grey Sa", degrey.build_Sa(), 5)):
    g = build_graph(pts)
    rel = ColourRelations(g, k)
    t = time.time()
    S, closed = forcing_set(rel, limit=500)
    if closed:
        S = shrink_forcing_set(rel, S)
    deg = max(g.degrees())
    print(f"{name}, k={k}: n={g.n}, max degree {deg}; "
          f"rho {'=' if closed else '>'} {len(S)}"
          f"   (a core of r at the hub needs rho <= {deg} + r)"
          f"  [{time.time()-t:.0f}s]", flush=True)
