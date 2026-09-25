"""The law has to be applied per vertex, not to the mean degree.

free@k counts the vertices with a spare colour, and a vertex of degree d is
free with probability about (k-1)(1 - 1/(k-1))^d if its neighbourhood looked
random.  Summing that over the vertices is not the same as evaluating it at the
mean degree: the function is CONVEX in d, so by Jensen

    mean_v  f(deg v)   >=   f(mean deg)

with equality only when every degree is the same.  A graph whose degrees are
spread -- a dense middle and a thin boundary, which is every union of rotated
copies here -- therefore has more free vertices than the mean-degree law
predicts, and comparing against that law manufactures an apparent "worse than
random".

That is exactly what happened.  The tuned chain at mean degree 13.35 measured
free@5 = 13.86 % against a mean-degree law of 8.60 %, reported as 0.62x --
worse than random, which is not a thing a real graph can be for this quantity
by any structural mechanism.  The per-vertex law is the honest null model, and
it needs no solving at all: only the degree sequence.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, math
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = HN_DIR
K = 5
# measured free@5, from scripts/freecurve.py
MEASURED = {
    "five_247_c.json": 0.15733, "five_247_b.json": 0.13617,
    "five_247.json": 0.12496, "five_tuned_1_1.json": 0.13861,
    "five_tuned_4_1.json": None, "five_tuned_1_3.json": None,
    "five_symmetric.json": None,
}
def f(d):
    return (K - 1) * ((K - 2) / (K - 1)) ** d

print(f"{'graph':24s} {'n':>6s} {'deg':>6s} {'sd':>5s} "
      f"{'mean-law':>9s} {'per-vertex':>11s} {'measured':>9s} {'ratio':>6s}")
for name, meas in MEASURED.items():
    try:
        d = json.load(open(f"{ROOT}/data/{name}"))
    except FileNotFoundError:
        continue
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    degs = [len(a) for a in g.adj]
    mu = sum(degs) / n
    sd = math.sqrt(sum((x - mu) ** 2 for x in degs) / n)
    mean_law = f(mu)
    pv = sum(f(x) for x in degs) / n
    r = f"{meas/pv:6.2f}" if meas else "     -"
    m = f"{100*meas:8.2f}%" if meas else "        -"
    print(f"{name:24s} {n:6d} {mu:6.2f} {sd:5.2f} {100*mean_law:8.3f}% "
          f"{100*pv:10.3f}% {m} {r}", flush=True)
