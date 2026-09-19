#!/usr/bin/env python3
"""Push the forced-pair core below 255 vertices.

A unit-distance graph on m vertices with a pair forced monochromatic at k=4
spindles into a genuine 5-chromatic unit-distance graph on at most 2m-1
vertices.  Parts' 509 is the smallest on record, so a core of 254 or fewer
beats it.

A single greedy pass reached 359, which is an upper bound and not a floor:
greedy deletion in one fixed order leaves whatever it happened to pass over.
This runs many randomised passes, each re-solving after every deletion, and
keeps whatever the best pass found.
"""
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pysat.formula import CNF
from pysat.solvers import Solver

from hn.certify import load_certificate, save_certificate
from hn.graph import build_graph

K = int(os.environ.get("HN_K", "4"))
PASSES = int(os.environ.get("HN_PASSES", "60"))
TARGET = int(os.environ.get("HN_TARGET", "254"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
SRC = os.environ.get("HN_SRC", os.path.join(OUT, "f4_core.json"))


def forced_cnf(g, k, pivot, targets):
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
    return cnf, a, sel


def main():
    pts, doc = load_certificate(SRC)
    g = build_graph(pts)
    print(f"loaded {g}", flush=True)

    # recover the forced pair: the pivot is the vertex whose distance-1/3 group
    # carries the forcing, so rediscover it rather than trusting a stored index
    from fractions import Fraction
    best = None
    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:40]:
        p = g.vertices[bp]
        tg = [j for j in range(g.n) if j != bp and p.dist2(g.vertices[j]).is_rational()
              and p.dist2(g.vertices[j]).c[0] == Fraction(1, 3)]
        if not tg:
            continue
        cnf, a, sel = forced_cnf(g, K, bp, tg)
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(assumptions=[sel] + [a(v) for v in range(g.n)]) is False
        s.delete()
        if ok:
            best = (bp, tg)
            break
    if best is None:
        print("no forced pair recovered", flush=True)
        return
    pivot, targets = best
    print(f"forced pair at pivot {pivot}, {len(targets)} targets at d2=1/3", flush=True)

    cnf, a, sel = forced_cnf(g, K, pivot, targets)
    protected = {pivot} | set(targets)
    t0 = time.time()
    overall = set(range(g.n))
    rnd = random.Random(12345)

    for p_ in range(PASSES):
        s = Solver(name="cd19", bootstrap_with=cnf)
        try:
            core = set(overall)
            order = [v for v in core if v not in protected]
            rnd.shuffle(order)
            for v in order:
                if v not in core:
                    continue
                trial = core - {v}
                if s.solve(assumptions=[sel] + [a(t) for t in trial]) is False:
                    core = trial
        finally:
            s.delete()
        if len(core) < len(overall):
            overall = core
            print(f"  pass {p_+1}: {len(overall)} vertices  "
                  f"-> spindles to <= {2*len(overall)-1}  [{time.time()-t0:.0f}s]", flush=True)
            sub = g.induced(sorted(overall))
            save_certificate(sub, os.path.join(OUT, "f4_core_min.json"), K,
                             f"forced monochromatic pair at k={K} on {len(overall)} vertices")
            if len(overall) <= TARGET:
                print(f"\n*** {len(overall)} vertices: spindles to <= {2*len(overall)-1}, "
                      f"BELOW THE RECORD OF 509 ***", flush=True)
                return
        elif p_ % 10 == 9:
            print(f"  pass {p_+1}: no improvement, still {len(overall)}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"finished {PASSES} passes at {len(overall)} vertices "
          f"(spindles to <= {2*len(overall)-1}) [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
