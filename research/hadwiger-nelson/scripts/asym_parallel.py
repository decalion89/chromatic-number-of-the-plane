#!/usr/bin/env python3
"""Narrow the forced core, evaluating candidates in parallel.

Making the search cheaper is not a side-issue in computer-assisted
mathematics, it is usually the work: the four colour theorem needed
discharging and unavoidable sets to make its search finite, the Boolean
Pythagorean triples needed cube-and-conquer, and Parts' 509-vertex graph is a
better minimisation of de Grey's construction rather than a new one.

Here the bottleneck is plain.  Candidate evaluation is independent -- build the
union, locate a forcing pivot, measure the core -- and the loop ran it on one
core out of four.  Each round is a pool map instead.

Workers rebuild their graph from the round's checkpoint rather than receiving
it, since a graph of exact Points carries thousands of Fractions and pickling
them costs more than the solve.
"""
import json
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import load_certificate, save_certificate
from hn.coloring import is_k_colorable
from hn.graph import build_graph
from hn.spindle import (SeparationTest, spindle_union_auto, triple_spindle_union)

K = int(os.environ.get("HN_K", "4"))
ROUNDS = int(os.environ.get("HN_ROUNDS", "14"))
FANOUT = int(os.environ.get("HN_FANOUT", "8"))
SCAN = int(os.environ.get("HN_SCAN", "12"))
BUDGET = int(os.environ.get("HN_BUDGET", "40000"))
WORKERS = int(os.environ.get("HN_WORKERS", "4"))
OUT = os.environ.get("HN_OUT", "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
SRC = os.environ.get("HN_SRC", os.path.join(OUT, "f4_core.json"))
THIRD = Fraction(1, 3)

_cache = {}


def graph_from(path):
    if path not in _cache:
        pts, _ = load_certificate(path)
        _cache[path] = build_graph(pts)
    return _cache[path]


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


def find_pivot(g, scan=SCAN):
    for bp in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:scan]:
        c = forced_core(g, bp, THIRD)
        if c:
            return bp, c
    return None, None


def evaluate(args):
    """One candidate, in its own process."""
    path, q = args
    g = graph_from(path)
    pq = g.vertices[q]
    tg = [j for j in range(g.n) if j != q and pq.dist2(g.vertices[j]).is_rational()
          and Fraction(1, 4) <= pq.dist2(g.vertices[j]).c[0] <= Fraction(8)]
    if not tg:
        return None
    try:
        cand, _ = spindle_union_auto(g, q, tg[0])
    except (ValueError, AssertionError):
        return None
    if cand.n > BUDGET:
        return None
    np_, ncore = find_pivot(cand)
    if np_ is None:
        return None
    return (len(ncore), q, cand.n, np_, sorted(ncore))


def main():
    t0 = time.time()
    path = SRC
    g = graph_from(path)
    pivot, core = find_pivot(g, scan=40)
    if pivot is None:
        print("no forced core in the source", flush=True)
        return
    print(f"start: {g}  pivot {pivot}  core {len(core)}  workers={WORKERS}", flush=True)
    history = [len(core)]

    with Pool(WORKERS) as pool:
        for rnd in range(1, ROUNDS + 1):
            cands = [v for v in sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:40]
                     if v != pivot][:FANOUT]
            results = [r for r in pool.map(evaluate, [(path, q) for q in cands]) if r]
            if not results:
                print(f"round {rnd}: no candidate produced forcing; stopping", flush=True)
                return
            results.sort()
            size, via, n, np_, ncore = results[0]
            print(f"round {rnd}: core {history[-1]} -> {size} via pivot {via}, "
                  f"n={n}  [{time.time()-t0:.0f}s]  "
                  f"(all: {[r[0] for r in results]})", flush=True)

            # rebuild the winner in this process and checkpoint it
            pq = g.vertices[via]
            tg = [j for j in range(g.n) if j != via and pq.dist2(g.vertices[j]).is_rational()
                  and Fraction(1, 4) <= pq.dist2(g.vertices[j]).c[0] <= Fraction(8)]
            g, _ = spindle_union_auto(g, via, tg[0])
            pivot, core = np_, ncore
            history.append(size)
            path = os.path.join(OUT, f"asymp_round{rnd}.json")
            save_certificate(g, path, K, f"asymmetric narrowing round {rnd}, core {size}")
            _cache[path] = g

            if size <= 2:
                print(f"  *** core {size} -- spindleable ***", flush=True)
                spun = (spindle_union_auto(g, pivot, core[0])[0] if size == 1
                        else triple_spindle_union(g, pivot))
                ok, _ = is_k_colorable(spun, K, timeout=7200)
                print(f"  spindled {spun}: {K}-colourable={ok}", flush=True)
                if ok is False:
                    save_certificate(spun, os.path.join(OUT, f"asymp_k{K}.json"), K,
                                     f"no proper {K}-colouring, from a core of size {size}")
                    print("  *** CERTIFICATE SAVED ***", flush=True)
                return
            if len(history) >= 4 and min(history[-3:]) >= history[-4]:
                print(f"  core stalled at {history[-4:]}; stopping", flush=True)
                return
    print(f"finished: {history} [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
