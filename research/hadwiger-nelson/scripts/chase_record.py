#!/usr/bin/env python3
"""Shrink de Grey's graph toward the smallest 5-chromatic unit-distance graph.

The record is Parts' 509 vertices (2020). The route is direct: de Grey's 1581
has no proper 4-colouring -- verified here by drat-trim -- so any subgraph that
still has none is a 5-chromatic unit-distance graph, and the smallest such
subgraph is what the record measures.

Honest odds: this is what Polymath16 did, with Heule spending on the order of
100000 CPU-hours to reach 529 and then 510. Beating 509 on four cores is very
unlikely. An independent minimisation landing anywhere near it still exercises
the whole pipeline on an object that matters, and the number is checkable.

Method: vertex selectors, iterated UNSAT cores, then randomised greedy deletion
passes. Every improvement is certified as it is found, so an interrupted run
still leaves its best result on disk.
"""
import os
import random
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pysat.formula import CNF
from pysat.solvers import Solver

from hn.certify import save_certificate
from hn.degrey import build_G

K = int(os.environ.get("HN_K", "4"))
PASSES = int(os.environ.get("HN_PASSES", "40"))
RECORD = int(os.environ.get("HN_RECORD", "509"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")


def main():
    t0 = time.time()
    g = build_G()
    n = g.n
    print(f"start: {g}  target: beat {RECORD}", flush=True)

    x = lambda v, c: 1 + v * K + c
    a = lambda v: 1 + n * K + v
    cnf = CNF()
    for v in range(n):
        cnf.append([-a(v)] + [x(v, c) for c in range(K)])
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-x(u, c), -x(v, c)])

    s = Solver(name="cd19", bootstrap_with=cnf)
    try:
        keep = set(range(n))
        if s.solve(assumptions=[a(v) for v in keep]) is not False:
            print("the graph IS 4-colourable -- nothing to shrink", flush=True)
            return
        core = {l - 1 - n * K for l in (s.get_core() or [])}
        print(f"  first UNSAT core: {len(core)}  [{time.time()-t0:.0f}s]", flush=True)
        for r in range(30):
            if s.solve(assumptions=[a(v) for v in core]) is not False:
                break
            new = {l - 1 - n * K for l in (s.get_core() or [])}
            if not new or len(new) >= len(core):
                break
            core = new
            print(f"  core round {r+1}: {len(core)}  [{time.time()-t0:.0f}s]", flush=True)

        best = set(core)
        rnd = random.Random(20260919)
        for p in range(PASSES):
            cur = set(best)
            order = list(cur)
            rnd.shuffle(order)
            dropped = 0
            for v in order:
                if v not in cur:
                    continue
                trial = cur - {v}
                if s.solve(assumptions=[a(t) for t in trial]) is False:
                    cur = trial
                    dropped += 1
            if len(cur) < len(best):
                best = cur
                sub = g.induced(sorted(best))
                print(f"  pass {p+1}: {len(best)} vertices, {sub.m} edges  "
                      f"[{time.time()-t0:.0f}s]"
                      f"{'   *** BELOW THE RECORD ***' if len(best) < RECORD else ''}",
                      flush=True)
                save_certificate(
                    sub, os.path.join(OUT, "smallest_5chromatic.json"), K,
                    f"chi(R^2) >= 5: {len(best)}-vertex unit-distance graph with no proper "
                    f"4-colouring", notes={"record_to_beat": RECORD, "pass": p + 1})
            elif p % 5 == 4:
                print(f"  pass {p+1}: no improvement, still {len(best)}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
        print(f"finished at {len(best)} vertices (record {RECORD}) "
              f"[{time.time()-t0:.0f}s]", flush=True)
    finally:
        s.delete()


if __name__ == "__main__":
    main()
