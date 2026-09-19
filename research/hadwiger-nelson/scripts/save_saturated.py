#!/usr/bin/env python3
"""Build the hole-saturated graph once and keep it.

Rebuilding costs two and a half minutes of unit-circle intersections every
time, and three separate searches now want the same graph. It is fully
deterministic -- same base, same holes in the same order -- so saving it loses
nothing and the certificate carries the coordinates exactly.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.certify import save_certificate
from hn.coloring import is_k_colorable
from hn.degrey import build_S, build_Sa, build_Y
from hn.graph import build_graph
from hn.mixed import deep_holes

K = int(os.environ.get("HN_K", "4"))
BASE = os.environ.get("HN_BASE", "Sa")
MIN_DEG = int(os.environ.get("HN_MINDEG", "4"))
LIMIT = int(os.environ.get("HN_LIMIT", "2500"))
OUT = os.environ.get(
    "HN_OUT",
    "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415"
    "/scratchpad/saturated_Sa.json")
BUILDERS = {"S": build_S, "Sa": build_Sa, "Y": build_Y}


def main() -> None:
    t0 = time.time()
    pts = list(BUILDERS[BASE]())
    g = build_graph(pts)
    holes = deep_holes(g, min_degree=MIN_DEG, limit=LIMIT)
    known = set(g.vertices)
    fresh = [(d, x) for d, x in holes if x not in known]
    lo, hi, best = 0, len(fresh), []
    while lo <= hi:
        mid = (lo + hi) // 2
        cand = [x for _, x in fresh[:mid]]
        if is_k_colorable(build_graph(pts + cand), K)[0]:
            best, lo = cand, mid + 1
        else:
            hi = mid - 1
    g2 = build_graph(pts + best)
    save_certificate(g2, OUT, k=K,
                     claim=f"{BASE} saturated with {len(best)} deep holes, "
                           f"still {K}-colourable")
    print(f"{BASE} + {len(best)} holes -> {g2}  saved to {OUT}  "
          f"[{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
