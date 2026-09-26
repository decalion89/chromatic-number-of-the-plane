#!/usr/bin/env python3
"""Hunt for a forced monochromatic pair under k colours, then spindle it.

If a pair (p, q) at a spindle-able distance is monochromatic in *every*
k-colouring of G, then G union rho_p(G) is not k-colourable at all.  For k=4
that reproduces chi(R^2) >= 5; for k=5 it would settle the open problem.

Forcing is monotone -- a pair forced in a subgraph stays forced in any
supergraph -- so the search wants the largest ball it can afford, not a clever
small one.
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
from hn.spindle import ForcedPairFinder, spindle_union

K = int(os.environ.get("HN_K", "4"))
CAP = int(os.environ.get("HN_CAP", "90000"))
OUT = os.environ.get("HN_OUT", "/tmp/hn")      # working directory for the output
os.makedirs(OUT, exist_ok=True)


def spindleable_pairs(basis, rows, pivot, max_per_distance=3):
    """Candidate (pivot, q) pairs whose spindle rotation stays in the field."""
    d2 = dist2_rational(basis, rows, pivot)
    finite = np.nonzero(~np.isnan(d2))[0]
    counts, out, ok_cache = {}, [], {}
    D2 = basis.D * basis.D
    dim = basis.dim
    diff = rows[finite] - rows[pivot]
    sq = basis._field_square(diff[:, :dim]) + basis._field_square(diff[:, dim:])
    for pos, j in enumerate(finite):
        if j == pivot:
            continue
        val = Fraction(int(sq[pos, 0]), D2)
        if val < Fraction(1, 4):
            continue  # rotating can never separate such a pair by 1
        if counts.get(val, 0) >= max_per_distance:
            continue
        ok = ok_cache.get(val)
        if ok is None:
            try:
                rotation_joining(val, FIELD)
                ok = True
            except ValueError:
                ok = False
            ok_cache[val] = ok
        if not ok:
            continue
        counts[val] = counts.get(val, 0) + 1
        out.append((pivot, int(j), val))
    return out, sorted(ok_cache and [v for v, k in ok_cache.items() if k])


def origin_index(rows):
    z = np.nonzero(~rows.any(axis=1))[0]
    return int(z[0]) if len(z) else 0


def attempt(m_max, steps, radius):
    U = unit_vectors(m_max=m_max)
    b = IntBasis.covering(U)
    t = time.time()
    rows = fast_walk(b, b.rows(U), steps=steps, radius=radius, cap=CAP)
    if len(rows) > CAP:
        print(f"m={m_max} st={steps} R={radius}: n={len(rows)} over cap, skip", flush=True)
        return None
    g = fast_graph_lazy(b, rows)
    g = g.k_core(K)
    if g.n == 0:
        return None
    rows2 = g.vertices.rows
    o = origin_index(rows2)
    cands, dists = spindleable_pairs(b, rows2, o)
    t_setup = time.time() - t
    t = time.time()
    finder = ForcedPairFinder(g, K, [(p, q) for p, q, _ in cands])
    forced = finder.run(verbose=False)
    print(
        f"m={m_max} st={steps} R={radius}: n={g.n} m={g.m} | dists={[str(d) for d in dists][:8]}"
        f" | {len(cands)} cands -> {len(forced)} forced  [setup {t_setup:.0f}s sat {time.time()-t:.0f}s]",
        flush=True,
    )
    if not forced:
        return None
    p, q = forced[0]
    print(f"  *** forced pair found; spindling ***", flush=True)
    sp = spindle_union(g, p, q, FIELD)
    sat, _ = is_k_colorable(sp, K, timeout=3600)
    print(f"  spindled: {sp}  {K}-colourable={sat}", flush=True)
    return sp if sat is False else None


def main():
    configs = []
    for m_max in (1, 2):
        for radius in (3.0, 4.0, 4.5, 5.0, 6.0):
            for steps in (6, 8, 10, 12):
                configs.append((m_max, steps, radius))
    for cfg in sorted(configs, key=lambda c: (c[2], c[1])):
        try:
            r = attempt(*cfg)
        except Exception as e:
            print(f"{cfg}: {type(e).__name__}: {e}", flush=True)
            continue
        if r is not None:
            print(f"\n*** NON-{K}-COLOURABLE UNIT-DISTANCE GRAPH FOUND ***", flush=True)
            import json
            with open(os.path.join(OUT, f"found_k{K}.json"), "w") as f:
                json.dump(r.to_dict(), f)
            return


if __name__ == "__main__":
    main()
