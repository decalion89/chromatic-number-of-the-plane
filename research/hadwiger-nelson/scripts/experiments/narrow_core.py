#!/usr/bin/env python3
"""Does tightening narrow the forced core, or is there a floor?

The spindle needs a narrow forced core: size 1, or size 2 at exactly
d^2 = 1/3.  On the configuration measured here the core is 34, so the question
that decides everything is whether tightening moves that number.

It can be asked at k=4, where forcing already exists, rather than waiting for
it to appear at k=5.  Each round unions the graph with a rotated copy and
re-measures the smallest forced core over a sample of pivots.  A number that
falls says narrowing is possible and at what rate; a number that sits still
says there is a structural floor, which would itself explain why the published
constructions are built rather than found.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.spindle import (SeparationTest, local_ball, spindle_union_auto,
                        triple_spindle_union)

K = int(os.environ.get("HN_K", "4"))
RADIUS = float(os.environ.get("HN_RADIUS", "3.0"))
SAMPLE = int(os.environ.get("HN_SAMPLE", "12"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "8"))
BUDGET = int(os.environ.get("HN_BUDGET", "40000"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def groups_for(g, pivot):
    p = g.vertices[pivot]
    out = {}
    for j in range(g.n):
        if j == pivot:
            continue
        d2 = p.dist2(g.vertices[j])
        if not d2.is_rational():
            continue
        v = d2.c[0]
        if Fraction(1, 4) <= v <= Fraction(40):
            out.setdefault(v, []).append(j)
    return out


def narrowest(g, sample):
    """(core size, pivot, d2, core) for the narrowest forced core found."""
    best = None
    for pivot in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:sample]:
        grp = groups_for(g, pivot)
        if not grp:
            continue
        allt = sorted({j for js in grp.values() for j in js})
        st = SeparationTest(g, K, pivot, allt)
        try:
            for val in sorted(grp):
                sep, core = st.run(subset=grp[val])
                if sep or not core:
                    continue
                if best is None or len(core) < best[0]:
                    best = (len(core), pivot, val, core)
                    if len(core) <= 2:
                        return best
        finally:
            st.close()
    return best


def main():
    t0 = time.time()
    g, _ = local_ball(build_G(), 
                      min(range(build_G().n),
                          key=lambda v: build_G().vertices[v].fx ** 2 + build_G().vertices[v].fy ** 2),
                      RADIUS)
    print(f"start: {g}  k={K}", flush=True)
    history = []
    for rnd in range(1, ROUNDS + 1):
        best = narrowest(g, SAMPLE)
        if best is None:
            print(f"round {rnd}: n={g.n} no forced core at all [{time.time()-t0:.0f}s]", flush=True)
            return
        size, pivot, val, core = best
        history.append(size)
        print(f"round {rnd}: n={g.n} m={g.m} narrowest core = {size} "
              f"(pivot {pivot}, d2={val}) [{time.time()-t0:.0f}s]", flush=True)

        if size == 1 or (size == 2 and val == Fraction(1, 3)):
            print(f"  *** SPINDLEABLE: core {size} at d2={val} ***", flush=True)
            spun = (spindle_union_auto(g, pivot, core[0])[0] if size == 1
                    else triple_spindle_union(g, pivot))
            ok, _ = is_k_colorable(spun, K, timeout=7200)
            print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
            if ok is False:
                save_certificate(spun, os.path.join(OUT, f"narrowed_k{K}.json"), K,
                                 f"no proper {K}-colouring; from a core of size {size}")
                print(f"  *** {spun.n} vertices, no {K}-colouring ***", flush=True)
                return

        if len(history) >= 4 and min(history[-3:]) >= history[-4]:
            print(f"  core has not narrowed in three rounds ({history[-4:]}) -- "
                  f"a floor, not a slow descent", flush=True)
            return

        tgt = groups_for(g, pivot).get(val)
        try:
            g2, fld = spindle_union_auto(g, pivot, tgt[0])
        except (ValueError, AssertionError) as e:
            print(f"  tightening failed: {e}", flush=True)
            return
        if g2.n > BUDGET:
            print(f"  budget exceeded at n={g2.n}", flush=True)
            return
        print(f"    tightened {g.n} -> {g2.n}", flush=True)
        g = g2
    print(f"finished {ROUNDS} rounds, history {history} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
