#!/usr/bin/env python3
"""Search for a forced disjunction, the weakest hypothesis the spindle needs.

For a pivot p and the vertices Q at a given distance from it, one SAT query
asks whether some k-colouring separates p from every element of Q.

  SAT   -> no disjunction over Q is forced.  Dead end, one query.
  UNSAT -> the disjunction over Q is forced; the UNSAT core is shrunk to a
           minimal forcing subset.

  core size 1                -> a plain forced pair: two rotated copies about
                                p give a non-k-colourable graph (any distance
                                with d^2 >= 1/4).
  core size 2 and d^2 = 1/3  -> three rotated copies do, by pigeonhole: the
                                images form a unit equilateral triangle.

Sizes above 2 are not usable, because more than three points pairwise at
distance 1 do not fit on a circle.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from hn.coloring import is_k_colorable
from hn.fast import IntBasis, dist2_rational, fast_graph_lazy, fast_walk
from hn.field import QSQRT3_11 as FIELD
from hn.generate import unit_vectors
from hn.geometry import rotation_joining
from hn.spindle import SeparationTest, spindle_union, triple_spindle_union

K = int(os.environ.get("HN_K", "4"))
CAP = int(os.environ.get("HN_CAP", "60000"))
NPIVOTS = int(os.environ.get("HN_PIVOTS", "6"))
OUT = os.environ.get("HN_OUT", "/tmp/hn")      # working directory for the output
os.makedirs(OUT, exist_ok=True)


def spindleable_distances(basis, rows, pivot, max_d2=40.0):
    d2 = dist2_rational(basis, rows, pivot)
    groups = {}
    D2 = basis.D * basis.D
    dim = basis.dim
    diff = rows - rows[pivot]
    sq = basis._field_square(diff[:, :dim]) + basis._field_square(diff[:, dim:])
    rational = np.ones(len(rows), dtype=bool)
    for m in range(1, dim):
        rational &= sq[:, m] == 0
    for j in np.nonzero(rational)[0]:
        if j == pivot:
            continue
        val = Fraction(int(sq[j, 0]), D2)
        if val < Fraction(1, 4) or val > max_d2:
            continue
        groups.setdefault(val, []).append(int(j))
    usable = {}
    for val, js in groups.items():
        try:
            rotation_joining(val, FIELD)
        except ValueError:
            continue
        usable[val] = js
    return usable


def probe(g, basis, rows, pivot, tag):
    usable = spindleable_distances(basis, rows, pivot)
    allt = sorted({j for js in usable.values() for j in js})
    if not allt:
        return None
    test = SeparationTest(g, K, pivot, allt)
    try:
        for val in sorted(usable):
            js = usable[val]
            sep, core = test.run(subset=js)
            if sep:
                continue
            print(f"  {tag} d2={val}: FORCED disjunction, core size {len(core)}", flush=True)
            if len(core) == 1:
                return ("pair", pivot, core[0], val)
            if len(core) == 2 and val == Fraction(1, 3):
                return ("triple", pivot, tuple(core), val)
            print(f"    core size {len(core)} at d2={val} is not usable", flush=True)
    finally:
        test.close()
    return None


def attempt(m_max, steps, radius):
    U = unit_vectors(m_max=m_max)
    b = IntBasis.covering(U)
    t = time.time()
    rows = fast_walk(b, b.rows(U), steps=steps, radius=radius, cap=CAP)
    if len(rows) > CAP:
        print(f"m={m_max} st={steps} R={radius}: n={len(rows)} over cap", flush=True)
        return None
    g = fast_graph_lazy(b, rows).k_core(K)
    if g.n == 0:
        return None
    rows2 = g.vertices.rows
    fx = b.floats(rows2)
    order = np.argsort(fx[:, 0] ** 2 + fx[:, 1] ** 2)      # pivots nearest the centre
    print(f"m={m_max} st={steps} R={radius}: {g} [{time.time()-t:.0f}s]", flush=True)
    for pivot in order[:NPIVOTS]:
        hit = probe(g, b, rows2, int(pivot), f"pivot#{int(pivot)}")
        if hit:
            return hit, g
    return None


def main():
    configs = []
    for m_max in (1, 2, 3):
        for radius in (2.0, 2.5, 3.0, 3.5):
            for steps in (5, 7, 9):
                configs.append((m_max, steps, radius))
    for cfg in configs:
        try:
            r = attempt(*cfg)
        except Exception as e:
            print(f"{cfg}: {type(e).__name__}: {e}", flush=True)
            continue
        if not r:
            continue
        hit, g = r
        kind = hit[0]
        print(f"\n*** {kind.upper()} SPINDLE AVAILABLE: {hit} ***", flush=True)
        if kind == "pair":
            spun = spindle_union(g, hit[1], hit[2], FIELD)
        else:
            spun = triple_spindle_union(g, hit[1])
        sat, _ = is_k_colorable(spun, K, timeout=3600)
        print(f"spindled: {spun}  {K}-colourable={sat}", flush=True)
        if sat is False:
            import json
            with open(os.path.join(OUT, f"found_k{K}.json"), "w") as f:
                json.dump(spun.to_dict(), f)
            print(f"*** NON-{K}-COLOURABLE GRAPH SAVED ***", flush=True)
            return


if __name__ == "__main__":
    main()
