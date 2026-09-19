#!/usr/bin/env python3
"""Measure f(4): the smallest unit-distance configuration with a pair forced
monochromatic under four colours.

The stakes, either way.  If a unit-distance graph on m vertices carries such a
pair, spindling it gives a genuine unit-distance graph on at most 2m-1 vertices
with no proper 4-colouring -- a 5-chromatic unit-distance graph.  The smallest
on record is Parts' 509.  So either f(4) >= 255, or this measurement beats that
record.  There is no third case.

Start where forcing is known to live.  The k=4 calibration found it in the ball
of radius 3.0 around the centre of de Grey's graph, at d^2 = 1/3, in two
seconds; the ball of radius 2.5 has none, so scanning it first only burns the
machine.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pysat.formula import CNF
from pysat.solvers import Solver

from hn.certify import save_certificate
from hn.degrey import build_G
from hn.spindle import local_ball

K = int(os.environ.get("HN_K", "4"))
RADIUS = float(os.environ.get("HN_RADIUS", "3.0"))
PIVOTS = int(os.environ.get("HN_PIVOTS", "25"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def build(g, k, pivot, targets):
    n = g.n
    x = lambda v, c: 1 + v * k + c
    a = lambda v: 1 + n * k + v
    sel = 1 + n * k + n
    cnf = CNF()
    for v in range(n):
        cnf.append([-a(v)] + [x(v, c) for c in range(k)])
    for u, v in g.edges():
        for c in range(k):
            cnf.append([-a(u), -a(v), -x(u, c), -x(v, c)])
    for q in targets:
        for c in range(k):
            cnf.append([-sel, -x(pivot, c), -x(q, c)])
    return cnf, a, sel, n * k


def main():
    t0 = time.time()
    g = build_G()
    centre = min(range(g.n), key=lambda v: g.vertices[v].fx ** 2 + g.vertices[v].fy ** 2)
    ball, _ = local_ball(g, centre, RADIUS)
    print(f"ball radius {RADIUS}: {ball}", flush=True)

    # The calibration found forcing here by sweeping the most connected
    # vertices, not the geometric centre.  Searching only the centre reports
    # nothing and contradicts a result we already had -- so sweep the same way.
    pivots = sorted(range(ball.n), key=lambda v: -len(ball.adj[v]))[:PIVOTS]
    for bp in pivots:
        p = ball.vertices[bp]
        groups = {}
        for j in range(ball.n):
            if j == bp:
                continue
            d2 = p.dist2(ball.vertices[j])
            if d2.is_rational() and Fraction(1, 4) <= d2.c[0] <= Fraction(40):
                groups.setdefault(d2.c[0], []).append(j)
        if find_forced(ball, bp, groups, t0):
            return
    print(f"no forced group over {len(pivots)} pivots [{time.time()-t0:.0f}s]", flush=True)


def find_forced(ball, bp, groups, t0):
    for val in sorted(groups):
        targets = groups[val]
        cnf, a, sel, off = build(ball, K, bp, targets)
        s = Solver(name="cd19", bootstrap_with=cnf)
        try:
            keep = set(range(ball.n))
            if s.solve(assumptions=[sel] + [a(v) for v in keep]) is not False:
                continue
            print(f"  pivot {bp} d2={val}: FORCED over {ball.n} vertices — minimising "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            core = {l - 1 - off for l in (s.get_core() or []) if l != sel} | {bp} | set(targets)
            print(f"    first UNSAT core: {len(core)}", flush=True)
            for r in range(40):
                if s.solve(assumptions=[sel] + [a(v) for v in core]) is not False:
                    break
                new = {l - 1 - off for l in (s.get_core() or []) if l != sel} | {bp} | set(targets)
                if not new or len(new) >= len(core):
                    break
                core = new
                print(f"    core round {r+1}: {len(core)}  [{time.time()-t0:.0f}s]", flush=True)
            dropped = 0
            for v in sorted(core, key=lambda v: len(ball.adj[v])):
                if v == bp or v in targets or v not in core:
                    continue
                trial = core - {v}
                if s.solve(assumptions=[sel] + [a(t) for t in trial]) is False:
                    core = trial
                    dropped += 1
                    if dropped % 25 == 0:
                        print(f"    greedy: {len(core)} left  [{time.time()-t0:.0f}s]", flush=True)
            sub = ball.induced(sorted(core))
            print(f"\n  f(4) <= {len(core)}   induced: {sub}   [{time.time()-t0:.0f}s]", flush=True)
            print(f"  spindling it would give a 5-chromatic unit-distance graph on "
                  f"<= {2*len(core)-1} vertices (record: 509)", flush=True)
            if 2 * len(core) - 1 < 509:
                print("  *** THAT BEATS THE RECORD ***", flush=True)
            save_certificate(sub, os.path.join(OUT, "f4_core.json"), K,
                             f"forced monochromatic pair at k={K} on {len(core)} vertices")
            return True
        finally:
            s.delete()
    return False


if __name__ == "__main__":
    main()
