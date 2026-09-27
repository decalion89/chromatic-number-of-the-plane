#!/usr/bin/env python3
"""Grow unit-distance graphs asymmetrically, because every generator here
symmetrises and symmetry is what floors the core.

The diagnosis. `hex_ball` is the triangular lattice, six-fold by definition.
`walk_ball` sums unit vectors omega^a, and that step set is six-fold. de Grey's
own construction says it outright: Sa is "all points obtained by rotating S
about the origin by multiples of 60 degrees and/or negating their
y-coordinates" -- explicit symmetrisation by a twelve-element group.

And the measurement: the eleven surviving targets sit at 38.55 + 60k and
69.6 + 60k, two six-fold orbits, and orbits are exactly what stops a core
narrowing, since a colouring argument cannot separate points a symmetry
permutes.

So the search space was symmetric by construction and could not contain a
core of 1 or 2, which needs targets that do not form a full orbit. This grows
graphs by accretion instead: start from a seed, repeatedly attach a point at
distance 1 from an existing vertex chosen at random, never symmetrise, and
measure the forced core as it goes.
"""
import os
import random
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.geometry import SPINDLE, eisenstein, origin
from hn.generate import unit_vectors
from hn.graph import build_graph
from hn.spindle import SeparationTest, spindle_union_auto, triple_spindle_union

K = int(os.environ.get("HN_K", "5"))
TRIALS = int(os.environ.get("HN_TRIALS", "40"))
TARGET_N = int(os.environ.get("HN_N", "900"))
SEED = int(os.environ.get("HN_SEED", "1"))
OUT = os.environ.get("HN_OUT", "/tmp/hn")      # working directory for the output
os.makedirs(OUT, exist_ok=True)


def smallest_core(g, scan=10):
    best = None
    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:scan]:
        p = g.vertices[bp]
        groups = {}
        for j in range(g.n):
            if j == bp:
                continue
            d2 = p.dist2(g.vertices[j])
            if not d2.is_rational():
                continue
            v = d2.c[0]
            if Fraction(1, 4) <= v <= Fraction(24):
                groups.setdefault(v, []).append(j)
        if not groups:
            continue
        st = SeparationTest(g, K, bp, sorted({j for js in groups.values() for j in js}))
        try:
            for val in sorted(groups):
                sep, core = st.run(subset=groups[val])
                if sep or not core:
                    continue
                if best is None or len(core) < best[0]:
                    best = (len(core), bp, val, core)
                    if len(core) <= 2:
                        return best
        finally:
            st.close()
    return best


def grow(rnd, n_target, units, pool=24):
    """Accrete points at distance 1, picking the one that creates most edges.

    Attaching at random breaks the symmetry but throws away the density with
    it: 699 vertices came out with 1070 edges, average degree 3.06, against
    10.0 for de Grey's graph -- and a sparse graph forces nothing, since it
    colours with room to spare. Symmetry was traded for density and both are
    needed.

    So each step samples a pool of candidate attachments and keeps whichever
    lands adjacent to the most existing points. The result stays asymmetric,
    nothing here symmetrises, while growing in dense clusters rather than
    branches.
    """
    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    pts = list(rh) + [SPINDLE(p) for p in rh]      # a Moser spindle as seed
    seen = set(pts)
    while len(pts) < n_target:
        best, best_deg = None, -1
        for _ in range(pool):
            base = pts[rnd.randrange(len(pts))]
            u = units[rnd.randrange(len(units))]
            q = base + u
            if q in seen:
                continue
            deg = sum(1 for w in pts if q.is_unit_apart(w))
            if deg > best_deg:
                best, best_deg = q, deg
        if best is None:
            continue
        seen.add(best)
        pts.append(best)
    return build_graph(pts)


def main():
    t0 = time.time()
    units = unit_vectors(m_max=2)
    print(f"asymmetric accretion at k={K}: {TRIALS} trials of ~{TARGET_N} vertices, "
          f"{len(units)} unit steps", flush=True)
    best = None
    for t in range(TRIALS):
        rnd = random.Random(SEED * 10000 + t)
        g = grow(rnd, TARGET_N, units)
        res = smallest_core(g)
        tag = "none" if res is None else f"CORE={res[0]} at d2={res[2]}"
        print(f"  trial {t+1:>3}: {g}  {tag}  [{time.time()-t0:.0f}s]", flush=True)
        if res is None:
            continue
        size, pivot, val, core = res
        if best is None or size < best[0]:
            best = (size, t + 1)
        if size <= 2:
            print(f"\n*** core {size} at k={K} ***", flush=True)
            spun = (spindle_union_auto(g, pivot, core[0])[0] if size == 1
                    else triple_spindle_union(g, pivot))
            ok, _ = is_k_colorable(spun, K, timeout=7200)
            print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
            if ok is False:
                save_certificate(spun, os.path.join(OUT, f"GROWN_k{K}.json"), K,
                                 f"chi(R^2) >= {K+1}: no proper {K}-colouring")
                print(f"\n*** chi(R^2) >= {K+1} -- CERTIFICATE SAVED ***", flush=True)
                return
    print(f"\nbest core over {TRIALS} trials: {best} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
