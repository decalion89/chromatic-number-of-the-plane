#!/usr/bin/env python3
"""Tighten the graph and measure whether it worked.

Every search so far was blind: sweep, get "nothing forced", sweep somewhere
else.  With separation effort as a score there is finally a direction to walk
in.  The loop is

    score(G)  = the largest conflict count over a sample of separation queries
    candidates = G u rho(G), spindled on the pairs that scored highest
    keep       = whichever candidate raises the score most

A score that climbs says the tightening is working.  A score that plateaus is
worth knowing too: it says this family of tightenings cannot reach forcing,
which no amount of further sweeping would have revealed.
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
from hn.spindle import SeparationDifficulty, spindle_union_auto

K = int(os.environ.get("HN_K", "5"))
MAXD2 = Fraction(os.environ.get("HN_MAXD2", "40"))
SAMPLE = int(os.environ.get("HN_SAMPLE", "25"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "8"))
FANOUT = int(os.environ.get("HN_FANOUT", "5"))
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


def score(g, sample=SAMPLE):
    """(max conflicts, best (pivot, d2, core) rows). Also reports any forcing."""
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:sample]
    rows = []
    for pivot in order:
        grp = groups_for(g, pivot)
        if not grp:
            continue
        td = SeparationDifficulty(g, K, pivot, sorted({j for js in grp.values() for j in js}))
        try:
            for val, sep, dc, dd, core in td.effort_ranking(grp):
                rows.append((dc, dd, pivot, val, sep, core))
        finally:
            td.close()
    rows.sort(reverse=True, key=lambda r: (r[0], r[1]))
    forced = [r for r in rows if not r[4]]
    return (rows[0][0] if rows else 0), rows, forced


def main():
    t0 = time.time()
    g = build_G()
    print(f"start: {g}", flush=True)
    for rnd in range(1, ROUNDS + 1):
        s, rows, forced = score(g)
        print(f"round {rnd}: {g.n} vertices, {g.m} edges, score={s} "
              f"(hardest d2={rows[0][3] if rows else '-'}), forced={len(forced)} "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if forced:
            dc, dd, pivot, val, sep, core = forced[0]
            print(f"  FORCED at pivot {pivot}, d2={val}, core {len(core)}", flush=True)
            if len(core) == 1:
                try:
                    spun, fld = spindle_union_auto(g, pivot, core[0])
                    ok, _ = is_k_colorable(spun, K, timeout=7200)
                    print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
                    if ok is False:
                        save_certificate(spun, os.path.join(OUT, f"found_k{K}.json"), K,
                                         f"chi(R^2) >= {K+1}: no proper {K}-colouring")
                        print(f"\n*** NON-{K}-COLOURABLE GRAPH ***", flush=True)
                        return
                except (ValueError, AssertionError) as e:
                    print(f"  spindle failed: {e}", flush=True)

        best = None
        for dc, dd, pivot, val, sep, core in rows[:FANOUT]:
            tgt = groups_for(g, pivot).get(val)
            if not tgt:
                continue
            try:
                cand, fld = spindle_union_auto(g, pivot, tgt[0])
            except (ValueError, AssertionError):
                continue
            if cand.n > 60000:
                continue
            cs, _, cforced = score(cand, sample=max(8, SAMPLE // 2))
            print(f"    candidate from pivot {pivot} d2={val}: {cand.n} vertices, "
                  f"score {dc} -> {cs}{'  FORCED' if cforced else ''}", flush=True)
            if best is None or cs > best[0]:
                best = (cs, cand)
        if best is None:
            print("  no candidate built; stopping", flush=True)
            return
        if best[0] <= s:
            print(f"  score did not improve ({s} -> {best[0]}); this family of "
                  f"tightenings plateaus here", flush=True)
        g = best[1]
    print(f"finished {ROUNDS} rounds [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
