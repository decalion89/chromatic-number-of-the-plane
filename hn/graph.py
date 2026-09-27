"""Unit-distance graphs: construction, reduction, and I/O.

Edges are found with a float spatial hash (cheap) and then *confirmed with
exact field arithmetic* (sound).  The float stage can only ever propose a
candidate pair; it never creates an edge on its own.  Coordinates here have
magnitude O(10) and denominators that are products of small primes, so double
precision carries ~1e-14 of error against a 1e-6 acceptance window -- eight
orders of margin.  `hn.certify` re-derives every edge exactly, without the
hash, so a published certificate never rests on this argument.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Set, Tuple

from .geometry import Point

__all__ = ["UnitDistanceGraph", "build_graph"]

_TOL = 1e-6


class UnitDistanceGraph:
    """A finite unit-distance graph: geometry plus the combinatorial shadow."""

    __slots__ = ("vertices", "adj", "_index", "rows", "basis")

    def __init__(self, vertices: Sequence[Point], adj: List[Set[int]]):
        self.vertices = vertices if hasattr(vertices, '_cache') else list(vertices)
        self.adj = adj
        self._index = None

    # -- basic accessors -------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.vertices)

    @property
    def m(self) -> int:
        return sum(len(a) for a in self.adj) // 2

    def edges(self) -> Iterable[Tuple[int, int]]:
        for u, nbrs in enumerate(self.adj):
            for v in nbrs:
                if u < v:
                    yield (u, v)

    def index_of(self, p: Point) -> int:
        if self._index is None:
            self._index = {q: i for i, q in enumerate(self.vertices)}
        return self._index[p]

    def degrees(self) -> List[int]:
        return [len(a) for a in self.adj]

    def __repr__(self) -> str:
        return f"<UnitDistanceGraph n={self.n} m={self.m}>"

    # -- reductions ------------------------------------------------------
    def induced(self, keep: Iterable[int]) -> "UnitDistanceGraph":
        keep = sorted(set(keep))
        remap = {old: new for new, old in enumerate(keep)}
        src = self.vertices
        if hasattr(src, "rows"):  # LazyVertices: stay lazy, slice the int rows
            import numpy as np

            from .fast import LazyVertices

            verts = LazyVertices(src.basis, src.rows[np.asarray(keep, dtype=np.int64)])
        else:
            verts = [src[i] for i in keep]
        adj = [set() for _ in keep]
        for old in keep:
            new = remap[old]
            for w in self.adj[old]:
                if w in remap:
                    adj[new].add(remap[w])
        return UnitDistanceGraph(verts, adj)

    def k_core(self, k: int) -> "UnitDistanceGraph":
        """Drop vertices of degree < k, repeatedly.

        Sound for colouring: a vertex with fewer than k neighbours can always
        be coloured last out of k colours, so G is k-colourable exactly when
        its k-core is.  This is the single most effective reduction we have --
        it routinely removes 90%+ of a generated vertex set.
        """
        deg = [len(a) for a in self.adj]
        alive = [True] * self.n
        stack = [v for v in range(self.n) if deg[v] < k]
        while stack:
            v = stack.pop()
            if not alive[v]:
                continue
            alive[v] = False
            for w in self.adj[v]:
                if alive[w]:
                    deg[w] -= 1
                    if deg[w] < k:
                        stack.append(w)
        return self.induced([v for v in range(self.n) if alive[v]])

    def largest_component(self) -> "UnitDistanceGraph":
        seen = [False] * self.n
        best: List[int] = []
        for s in range(self.n):
            if seen[s]:
                continue
            comp, stack, seen[s] = [], [s], True
            while stack:
                v = stack.pop()
                comp.append(v)
                for w in self.adj[v]:
                    if not seen[w]:
                        seen[w] = True
                        stack.append(w)
            if len(comp) > len(best):
                best = comp
        return self.induced(best)

    # -- heuristics ------------------------------------------------------
    def greedy_coloring(self, order: Sequence[int] | None = None) -> List[int]:
        if order is None:
            order = sorted(range(self.n), key=lambda v: -len(self.adj[v]))
        color = [-1] * self.n
        for v in order:
            used = {color[w] for w in self.adj[v] if color[w] >= 0}
            c = 0
            while c in used:
                c += 1
            color[v] = c
        return color

    def degeneracy(self) -> int:
        """Max over the elimination order of the degree at removal time.
        An upper bound on the chromatic number minus one."""
        deg = [len(a) for a in self.adj]
        alive = [True] * self.n
        buckets = defaultdict(set)
        for v, d in enumerate(deg):
            buckets[d].add(v)
        best = 0
        for _ in range(self.n):
            d = 0
            while d < self.n and not buckets[d]:
                d += 1
            if d >= self.n:
                break
            v = buckets[d].pop()
            alive[v] = False
            best = max(best, d)
            for w in self.adj[v]:
                if alive[w]:
                    buckets[deg[w]].discard(w)
                    deg[w] -= 1
                    buckets[deg[w]].add(w)
        return best

    def triangles_through(self, v: int) -> List[Tuple[int, int, int]]:
        out = []
        nbrs = sorted(self.adj[v])
        for i, a in enumerate(nbrs):
            for b in nbrs[i + 1 :]:
                if b in self.adj[a]:
                    out.append((v, a, b))
        return out

    def find_clique(self, size: int) -> List[int] | None:
        """A small greedy clique, used to break colour symmetry in the SAT encoding."""
        order = sorted(range(self.n), key=lambda v: -len(self.adj[v]))
        for v in order[: min(self.n, 400)]:
            clique = [v]
            cands = set(self.adj[v])
            while cands and len(clique) < size:
                w = max(cands, key=lambda u: len(self.adj[u] & cands))
                clique.append(w)
                cands &= self.adj[w]
            if len(clique) >= size:
                return clique[:size]
        return None

    # -- serialisation ---------------------------------------------------
    def to_dict(self) -> dict:
        f = self.vertices[0].field if self.vertices else None
        return {
            "field": list(f.gens) if f else [],
            "vertices": [
                [[str(c) for c in v.x.c], [str(c) for c in v.y.c]] for v in self.vertices
            ],
            "edges": [list(e) for e in self.edges()],
        }


def build_graph(points: Iterable[Point], verbose: bool = False) -> UnitDistanceGraph:
    """Build the unit-distance graph on `points` (duplicates removed)."""
    verts: List[Point] = []
    seen = set()
    for p in points:
        if p not in seen:
            seen.add(p)
            verts.append(p)
    n = len(verts)
    adj: List[Set[int]] = [set() for _ in range(n)]

    # Bucket by unit cells; a unit-distance pair lies in the same or an
    # adjacent cell, so 9 cells per point suffice.
    cells: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, p in enumerate(verts):
        cells[(math.floor(p.fx), math.floor(p.fy))].append(i)

    exact_checks = 0
    for (cx, cy), bucket in cells.items():
        neigh: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                if (dx, dy) >= (0, 0) or True:
                    neigh.extend(cells.get((cx + dx, cy + dy), ()))
        for i in bucket:
            pi = verts[i]
            xi, yi = pi.fx, pi.fy
            for j in neigh:
                if j <= i:
                    continue
                pj = verts[j]
                d2 = (xi - pj.fx) ** 2 + (yi - pj.fy) ** 2
                if abs(d2 - 1.0) < _TOL:
                    exact_checks += 1
                    if pi.is_unit_apart(pj):
                        adj[i].add(j)
                        adj[j].add(i)
    g = UnitDistanceGraph(verts, adj)
    if verbose:
        print(f"  built {g} ({exact_checks} exact distance checks)")
    return g
