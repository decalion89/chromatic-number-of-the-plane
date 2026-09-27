#!/usr/bin/env python3
"""Hill-climb locally, so the machine lasts as long as the mathematics.

The global version doubled the vertex count every round -- 1581, 3008, and by
the fourth it would be scoring for longer than it spent searching.  It hit a
compute wall, not a mathematical one.

Forcing is local.  The pair under test lies within sqrt(40) of the pivot, and
what decides it is the structure nearby, so there is no need to carry the far
side of the graph along.  Each round therefore: take the ball around the
hardest pivot, tighten only that, peel it back with the k-core -- which is
sound for forcing, not merely for colourability, since a vertex of degree
below k is always colourable last -- and score again.
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
from hn.spindle import (SeparationDifficulty, core_preserving_forcing, local_ball,
                        spindle_union_auto)

K = int(os.environ.get("HN_K", "5"))
MAXD2 = Fraction(os.environ.get("HN_MAXD2", "40"))
RADIUS = float(os.environ.get("HN_RADIUS", "3.0"))
SAMPLE = int(os.environ.get("HN_SAMPLE", "20"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "40"))
BUDGET = int(os.environ.get("HN_BUDGET", "40000"))
CONFLICTS = int(os.environ.get("HN_CONFLICTS", "200000"))
SAVE = os.environ.get("HN_SAVE")
OUT = os.environ.get("HN_OUT", "/tmp/hn")


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
        if Fraction(1, 4) <= v <= MAXD2:
            out.setdefault(v, []).append(j)
    return out


def best_rows(g, sample):
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:sample]
    rows = []
    for pivot in order:
        grp = groups_for(g, pivot)
        if not grp:
            continue
        td = SeparationDifficulty(g, K, pivot, sorted({j for js in grp.values() for j in js}))
        try:
            for val, sep, dc, dd, core in td.effort_ranking(grp, conflict_budget=CONFLICTS):
                rows.append((dc, dd, pivot, val, sep, core))
        finally:
            td.close()
    rows.sort(reverse=True, key=lambda r: (r[0], r[1]))
    return rows


def report_forced(g, pivot, val, core):
    print(f"  FORCED at pivot {pivot}, d2={val}, core size {len(core)}", flush=True)
    if len(core) != 1:
        return False
    try:
        spun, fld = spindle_union_auto(g, pivot, core[0])
    except (ValueError, AssertionError) as e:
        print(f"    spindle failed: {e}", flush=True)
        return False
    ok, _ = is_k_colorable(spun, K, timeout=7200)
    print(f"    spindled {spun} over {fld}: {K}-colourable={ok}", flush=True)
    if ok is False:
        save_certificate(spun, os.path.join(OUT, f"found_k{K}.json"), K,
                         f"chi(R^2) >= {K+1}: no proper {K}-colouring")
        print(f"\n*** NON-{K}-COLOURABLE GRAPH ***", flush=True)
        return True
    return False


def main():
    t0 = time.time()
    g = build_G()
    print(f"start: {g}  radius={RADIUS} budget={BUDGET}", flush=True)
    history = []
    for rnd in range(1, ROUNDS + 1):
        rows = best_rows(g, SAMPLE)
        if not rows:
            print("  no scorable pivot; stopping", flush=True)
            return
        forced = [r for r in rows if r[4] is False]
        score, _, pivot, val, _, _ = rows[0]
        history.append(score)
        print(f"round {rnd}: n={g.n} m={g.m} score={score} (pivot {pivot}, d2={val}) "
              f"forced={len(forced)} [{time.time()-t0:.0f}s]", flush=True)
        if forced:
            if report_forced(g, forced[0][2], forced[0][3], forced[0][5]):
                return
        if len(history) >= 4 and max(history[-3:]) <= history[-4]:
            print(f"  score has not improved in three rounds ({history[-4:]}); "
                  f"this family of tightenings plateaus here", flush=True)
            return

        # tighten locally around the hardest pivot, then peel
        ball, bp = local_ball(g, pivot, RADIUS)
        tgts = groups_for(ball, bp).get(val)
        if not tgts:
            print("  hardest distance vanished from the ball; widening", flush=True)
            ball, bp = local_ball(g, pivot, RADIUS * 1.5)
            tgts = groups_for(ball, bp).get(val)
            if not tgts:
                print("  still absent; stopping", flush=True)
                return
        try:
            tight, fld = spindle_union_auto(ball, bp, tgts[0])
        except (ValueError, AssertionError) as e:
            print(f"  tightening failed: {e}", flush=True)
            return
        peeled = core_preserving_forcing(tight, K, 0, [])
        nxt = peeled[0] if peeled else tight
        print(f"    ball n={ball.n} -> tightened n={tight.n} -> peeled n={nxt.n}", flush=True)
        if SAVE:
            save_certificate(nxt, os.path.join(OUT, f"{SAVE}_round{rnd}.json"), K,
                             f"hill-climb round {rnd}, score {score}, not yet {K}-uncolourable")
        if nxt.n > BUDGET:
            print(f"  budget {BUDGET} exceeded at n={nxt.n}; stopping", flush=True)
            return
        g = nxt
    print(f"finished {ROUNDS} rounds [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
