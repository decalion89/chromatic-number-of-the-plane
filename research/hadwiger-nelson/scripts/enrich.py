#!/usr/bin/env python3
"""Build the most constrained graph that still colours, and keep it.

Three sources of extra points, each added while the graph still k-colours --
past that the forcing question goes vacuous, which is the trap already caught
twice here.

*Unit-triangle centroids.* The two-orbit block needs a target at squared
distance 1/3 from the pivot, and the constructions here are built from unit
steps and rarely land that close. Three points pairwise one apart have a
centroid exactly 1/sqrt(3) from each, so every unit triangle donates a point
with three legs on the classical circle. On Sa: 469 centroids, 288 of them new,
and every pivot tested then has a leg.

*Deep holes.* Points of the plane with many graph vertices exactly one away,
found as unit-circle intersections. A pivot need not be a vertex, and the
useful ones are those whose neighbourhood is most tightly pinned.

*And nothing else*, because the binary search stops where the colouring does.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_G, build_S, build_Sa, build_Y
from hn.graph import build_graph
from hn.mixed import deep_holes, unit_triangle_centroids

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
MIN_DEG = int(os.environ.get("HN_MINDEG", "4"))
LIMIT = int(os.environ.get("HN_LIMIT", "2500"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/enriched.json")
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y, "G": build_G}


def most_that_colours(pts, extra, k):
    """Binary search on how many extras the graph takes while still colouring."""
    lo, hi, best = 0, len(extra), []
    while lo <= hi:
        mid = (lo + hi) // 2
        cand = extra[:mid]
        if is_k_colorable(build_graph(pts + cand), k)[0]:
            best, lo = cand, mid + 1
        else:
            hi = mid - 1
    return best


def main() -> None:
    t0 = time.time()
    base = BUILDERS[BASE]()
    pts = list(base.vertices if hasattr(base, "vertices") else base)
    g = build_graph(pts)
    print(f"{BASE}: {g}", flush=True)

    known = set(g.vertices)
    cents = [c for c in unit_triangle_centroids(g) if c not in known]
    take = most_that_colours(pts, cents, K)
    pts = pts + take
    g = build_graph(pts)
    print(f"  + {len(take)} of {len(cents)} centroids: {g}  "
          f"[{time.time() - t0:.0f}s]", flush=True)

    known = set(g.vertices)
    holes = [x for _d, x in deep_holes(g, min_degree=MIN_DEG, limit=LIMIT)
             if x not in known]
    take = most_that_colours(pts, holes, K)
    pts = pts + take
    g = build_graph(pts)
    print(f"  + {len(take)} of {len(holes)} holes: {g}  "
          f"[{time.time() - t0:.0f}s]", flush=True)

    save_certificate(g, OUT, k=K,
                     claim=f"{BASE} enriched with centroids and holes, "
                           f"still {K}-colourable")
    print(f"  saved to {OUT}", flush=True)


if __name__ == "__main__":
    main()
