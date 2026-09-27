#!/usr/bin/env python3
"""Tighten asymmetrically, because symmetry is what holds the core wide.

The forced core on the measured k=4 configuration is 34 of 36 targets, and all
36 form a closed orbit under the 60-degree rotation about the pivot.  A
colouring argument cannot single out targets a symmetry permutes, so that
orbit is a floor.  The rotation is not a full automorphism -- 280 of 359
vertices return -- and that partial asymmetry is exactly what bought the two
exclusions from 36 down to 34.

So the tightening used until now was self-defeating: unioning a graph with
copies rotated about the *same* pivot raises symmetry and pushes the core the
wrong way, which is how the effort score could climb through 88, 2110, 22697
while the core never moved.

This tightens about *other* pivots instead, and picks the union that leaves
fewest vertices fixed by the rotation about the forcing pivot -- the most
symmetry broken per round.  The objective is the core size, which has to reach
1, or 2 at exactly d^2 = 1/3.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import load_certificate, save_certificate
from hn.coloring import is_k_colorable
from hn.geometry import _rot60
from hn.graph import build_graph
from hn.spindle import (SeparationTest, spindle_union_auto, triple_spindle_union)

K = int(os.environ.get("HN_K", "4"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "10"))
FANOUT = int(os.environ.get("HN_FANOUT", "6"))
BUDGET = int(os.environ.get("HN_BUDGET", "30000"))
SCAN = int(os.environ.get("HN_SCAN", "12"))   # pivots examined per candidate
SAVE = os.environ.get("HN_SAVE", "asym")
OUT = os.environ.get("HN_OUT", "/tmp/hn")
SRC = os.environ.get("HN_SRC", os.path.join(OUT, "f4_core.json"))


def symmetry_fraction(g, pivot):
    """How much of the graph the 60-degree rotation about `pivot` preserves.
    Lower is more asymmetric, which is what buys exclusions."""
    turn = _rot60(g.vertices[pivot].field).about(g.vertices[pivot])
    inside = set(g.vertices)
    return sum(1 for v in g.vertices if turn(v) in inside) / g.n


def forced_core(g, pivot, d2):
    p = g.vertices[pivot]
    tg = [j for j in range(g.n) if j != pivot and p.dist2(g.vertices[j]).is_rational()
          and p.dist2(g.vertices[j]).c[0] == d2]
    if not tg:
        return None
    st = SeparationTest(g, K, pivot, tg)
    try:
        sep, core = st.run(subset=tg)
    finally:
        st.close()
    return None if sep else core


def find_pivot(g, d2, scan=None):
    """Locate a forcing pivot.  Scanning every high-degree vertex costs 40 SAT
    calls per candidate, which at 1400 vertices dominates the round; after the
    first tightening the forcing pivot stays in much the same place, so a
    shorter scan buys most of the same answers at a sixth of the price."""
    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:(scan or SCAN)]:
        c = forced_core(g, bp, d2)
        if c:
            return bp, c
    return None, None


def main():
    t0 = time.time()
    pts, _ = load_certificate(SRC)
    g = build_graph(pts)
    third = Fraction(1, 3)
    pivot, core = find_pivot(g, third)
    if pivot is None:
        print("no forced core in the source", flush=True)
        return
    print(f"start: {g}  pivot {pivot}  core {len(core)}  "
          f"symmetry {symmetry_fraction(g, pivot):.3f}", flush=True)

    history = [len(core)]
    for rnd in range(1, ROUNDS + 1):
        best = None
        # tighten about pivots OTHER than the forcing one, to break its symmetry
        cands = [v for v in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:40] if v != pivot]
        for q in cands[:FANOUT]:
            pq = g.vertices[q]
            tg = [j for j in range(g.n) if j != q and pq.dist2(g.vertices[j]).is_rational()
                  and Fraction(1, 4) <= pq.dist2(g.vertices[j]).c[0] <= Fraction(8)]
            if not tg:
                continue
            try:
                cand, _ = spindle_union_auto(g, q, tg[0])
            except (ValueError, AssertionError):
                continue
            if cand.n > BUDGET:
                continue
            np_, ncore = find_pivot(cand, third)
            if np_ is None:
                continue
            sym = symmetry_fraction(cand, np_)
            print(f"    via pivot {q}: n={cand.n} core={len(ncore)} sym={sym:.3f}", flush=True)
            key = (len(ncore), sym)
            if best is None or key < best[0]:
                best = (key, cand, np_, ncore)
        if best is None:
            print("  no asymmetric candidate; stopping", flush=True)
            return
        (size, sym), g, pivot, core = best
        history.append(size)
        # checkpoint: without this an interrupted run loses every round it did
        save_certificate(g, os.path.join(OUT, f"{SAVE}_round{rnd}.json"), K,
                         f"asymmetric narrowing round {rnd}, forced core {size}")
        print(f"round {rnd}: n={g.n} core={size} symmetry={sym:.3f} "
              f"history={history} [{time.time()-t0:.0f}s]", flush=True)

        if size == 1 or (size == 2):
            print(f"  *** core {size} -- spindleable ***", flush=True)
            spun = (spindle_union_auto(g, pivot, core[0])[0] if size == 1
                    else triple_spindle_union(g, pivot))
            ok, _ = is_k_colorable(spun, K, timeout=7200)
            print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
            if ok is False:
                save_certificate(spun, os.path.join(OUT, f"asym_k{K}.json"), K,
                                 f"no proper {K}-colouring, from a core of size {size}")
                print("  *** CERTIFICATE SAVED ***", flush=True)
                return
        if len(history) >= 4 and min(history[-3:]) >= history[-4]:
            print(f"  core has not narrowed in three rounds -- asymmetry is not "
                  f"buying exclusions here either", flush=True)
            return
    print(f"finished: history {history} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
