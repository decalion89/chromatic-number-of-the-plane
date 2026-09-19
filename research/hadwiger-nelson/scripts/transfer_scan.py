#!/usr/bin/env python3
"""Look for an interface whose transfer relation runs dry.

The spindle argument glues copies at a point. This glues them along a set: if
W' is congruent to W, a chain of copies with each copy's W' identified with
the next copy's W is a legitimate unit-distance graph, and the colourings that
survive it are exactly the forward-reachable patterns of the transfer
relation. When that set empties after m steps, the chain of m copies has no
k-colouring.

Sound for a reason worth stating: copies placed by isometries may meet
elsewhere and create unit-distance edges the relation never accounted for, but
an edge only removes colourings. The computed reachable set is a superset of
the true one, so emptying is conclusive and non-emptying is not.

Interfaces are scored before they are tested. A W whose realisable patterns
are few is already tightly constrained, and a relation where each pattern
reaches few others is already close to running dry.
"""
import os
import random
import sys
import time
from itertools import combinations

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import load_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G, build_S, build_Sa, build_Y
from hn.graph import build_graph
from hn.transfer import (chain_length, congruent_pairs, realisable_patterns,
                         transfer_relation)

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
SIZE = int(os.environ.get("HN_SIZE", "3"))
TRIES = int(os.environ.get("HN_TRIES", "400"))
MATES = int(os.environ.get("HN_MATES", "6"))
SEED = int(os.environ.get("HN_SEED", "5"))
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y, "G": build_G}


def load(name):
    if name in BUILDERS:
        g = BUILDERS[name]()
        return g if hasattr(g, "vertices") else build_graph(g)
    pts, _doc = load_certificate(name)
    return build_graph(pts)


def main() -> None:
    g = load(BASE)
    ok = is_k_colorable(g, K)[0]
    print(f"{BASE}: {g}  k={K}, {K}-colourable: {ok}", flush=True)
    if not ok:
        print("  vacuous: no k-colouring, every pattern is unrealisable",
              flush=True)
        return
    rnd = random.Random(SEED)
    # tight interfaces first: a set with adjacent points and high degree
    hub = sorted(range(g.n), key=lambda v: -len(g.adj[v]))[:60]
    t0 = time.time()
    best = None
    for _ in range(TRIES):
        seed = rnd.choice(hub)
        nbrs = list(g.adj[seed])
        if len(nbrs) < SIZE - 1:
            continue
        W = [seed] + rnd.sample(nbrs, SIZE - 1)
        pats = realisable_patterns(g, K, W)
        if not pats:
            continue
        mates = congruent_pairs(g, W, limit=MATES)
        for Wp in mates:
            rel = transfer_relation(g, K, W, Wp)
            if not rel:
                print(f"  W={W} W'={list(Wp)}: relation empty -- a single "
                      f"copy already fails  [{time.time() - t0:.0f}s]",
                      flush=True)
                return
            fan = sum(len(v) for v in rel.values()) / len(rel)
            m = chain_length(rel)
            if m is not None:
                print(f"  *** W={W} W'={list(Wp)}: chain empties after {m} "
                      f"copies ***  [{time.time() - t0:.0f}s]", flush=True)
                return
            if best is None or fan < best[0]:
                best = (fan, len(rel), W, list(Wp))
                print(f"  W={W} W'={list(Wp)}: {len(rel)} live patterns, "
                      f"mean fan-out {fan:.2f}  [{time.time() - t0:.0f}s]",
                      flush=True)
    print(f"  tightest seen: fan-out {best[0]:.2f} over {best[1]} patterns"
          if best else "  nothing measurable", flush=True)


if __name__ == "__main__":
    main()
