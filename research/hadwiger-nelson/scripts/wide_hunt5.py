#!/usr/bin/env python3
"""Breadth at five colours, where breadth is cheap.

The asymmetry that has shaped this whole search works in our favour here.
A 5-colouring of these graphs is found in milliseconds while the k=4 UNSAT on
the same object takes twenty minutes, so scanning many configurations at k=5
costs about ten seconds each -- and what we need is one configuration born
with a forced core of at most 2.

Includes the tightened graphs left by the asymmetric narrowing runs. Those
were built to be hard at k=4 and have never been looked at with five colours;
they are the most constrained objects this project has produced.
"""
import glob
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import load_certificate, save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G, build_Sa, build_Y
from hn.fast import IntBasis, fast_graph_complete, fast_walk
from hn.generate import hex_ball, minkowski, unit_vectors
from hn.graph import build_graph
from hn.spindle import (SeparationTest, local_ball, spindle_union_auto,
                        triple_spindle_union)

K = int(os.environ.get("HN_K", "5"))
SCAN = int(os.environ.get("HN_SCAN", "14"))
CAP = int(os.environ.get("HN_CAP", "45000"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def smallest_core(g, scan=SCAN):
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
            if Fraction(1, 4) <= v <= Fraction(40):
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


def configurations():
    # the tightened checkpoints first: most constrained, never seen at k=5
    for path in sorted(glob.glob(os.path.join(OUT, "asymp_round*.json"))):
        try:
            pts, _ = load_certificate(path)
            yield os.path.basename(path), build_graph(pts)
        except Exception:
            continue
    for path in (os.path.join(OUT, "f4_core.json"),):
        if os.path.exists(path):
            pts, _ = load_certificate(path)
            yield os.path.basename(path), build_graph(pts)

    G = build_G()
    yield "de Grey G", G
    yield "Sa", build_graph(build_Sa())
    yield "Y", build_graph(build_Y())

    # spindled unions of G, which nothing has examined at k=5
    centre = min(range(G.n), key=lambda v: G.vertices[v].fx ** 2 + G.vertices[v].fy ** 2)
    for d2 in (Fraction(1, 3), Fraction(1, 4), Fraction(1, 2), Fraction(3)):
        p = G.vertices[centre]
        tg = [j for j in range(G.n) if j != centre and p.dist2(G.vertices[j]).is_rational()
              and p.dist2(G.vertices[j]).c[0] == d2]
        if not tg:
            continue
        try:
            u, _ = spindle_union_auto(G, centre, tg[0])
            if u.n <= CAP:
                yield f"G u rho(G) d2={d2}", u
        except (ValueError, AssertionError):
            pass

    U = unit_vectors(m_max=1)
    b = IntBasis.covering(U)
    for steps, rad in ((5, 2.0), (6, 2.0), (7, 2.5), (8, 2.5), (8, 3.0), (9, 3.0)):
        rows = fast_walk(b, b.rows(U), steps=steps, radius=rad, cap=CAP)
        if len(rows) <= CAP:
            yield f"walk s={steps} R={rad}", fast_graph_complete(b, rows)

    for r in (4, 5, 6):
        hb = hex_ball(r)
        yield f"hex r={r}", build_graph(hb)
        if r <= 4:
            yield f"hex r={r} + spindled", build_graph(minkowski(hb, hex_ball(1)))


def main():
    t0 = time.time()
    print(f"wide hunt at k={K}: looking for a forced core <= 2", flush=True)
    seen = 0
    best = None
    for name, g in configurations():
        t = time.time()
        try:
            res = smallest_core(g)
        except Exception as e:
            print(f"  {name:<26} {type(e).__name__}: {e}", flush=True)
            continue
        seen += 1
        if res is None:
            print(f"  {name:<26} n={g.n:>6}  none  [{time.time()-t:.0f}s]", flush=True)
            continue
        size, pivot, val, core = res
        print(f"  {name:<26} n={g.n:>6}  CORE={size} at d2={val}  [{time.time()-t:.0f}s]", flush=True)
        if best is None or size < best[0]:
            best = (size, name)
        if size <= 2:
            print(f"\n*** core {size} at k={K} -- spindleable ***", flush=True)
            spun = (spindle_union_auto(g, pivot, core[0])[0] if size == 1
                    else triple_spindle_union(g, pivot))
            ok, _ = is_k_colorable(spun, K, timeout=7200)
            print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
            if ok is False:
                save_certificate(spun, os.path.join(OUT, f"WIDE_k{K}.json"), K,
                                 f"chi(R^2) >= {K+1}: no proper {K}-colouring")
                print(f"\n*** chi(R^2) >= {K+1} -- CERTIFICATE SAVED ***", flush=True)
                return
    print(f"\n{seen} configurations, best core {best} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
