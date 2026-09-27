"""A vectorised integer representation, for searching at scale.

Exact `Fraction` arithmetic is correct but far too slow past ~10^4 points.
Every point we generate is a sum of unit vectors drawn from one finite set, so
all coordinates share a common denominator D.  Multiplying through by D turns
the whole search into int64 arithmetic on a numpy array of shape (N, 2*dim):
the first `dim` columns are the numerators of x over the field basis, the rest
are y.  Generation becomes broadcast addition and `np.unique`.

Soundness note -- the important one.  Edges are found by looking up p + u for
u in a supplied set W of vectors each of which has been *exactly verified* to
have length 1.  If W misses some unit vector, we build a subgraph of the true
unit-distance graph on those points.  That costs power, never correctness: a
subgraph of a unit-distance graph is itself a unit-distance graph, so a
non-k-colourable subgraph proves exactly as much as the full graph would.
Adding a wrong edge would be fatal, and cannot happen -- every u in W is
checked with exact field arithmetic before it is ever used.
"""

from __future__ import annotations

import math
from fractions import Fraction
from collections.abc import Sequence
from typing import Iterable, List, Tuple

import numpy as np

from .field import Field, FieldElement, QSQRT3_11
from .geometry import Point
from .graph import UnitDistanceGraph

__all__ = ["IntBasis", "fast_walk", "fast_edges", "fast_graph",
           "fast_edges_complete", "fast_graph_complete", "fast_graph_lazy",
           "LazyVertices", "dist2_rational"]


class IntBasis:
    """Shared denominator + conversions between exact Points and int64 rows."""

    def __init__(self, field: Field, denom: int):
        self.field = field
        self.dim = field.dim
        self.D = int(denom)
        self.prod = np.array(field._prod, dtype=np.int64)
        self.sqrts = np.array([math.sqrt(p) for p in field._prod], dtype=np.float64)

    @classmethod
    def covering(cls, points) -> "IntBasis":
        """The smallest common denominator for these points."""
        field = points[0].field
        d = 1
        for p in points:
            for coef in list(p.x.c) + list(p.y.c):
                d = d * coef.denominator // math.gcd(d, coef.denominator)
        return cls(field, d)

    def rows(self, points) -> np.ndarray:
        out = np.zeros((len(points), 2 * self.dim), dtype=np.int64)
        for i, p in enumerate(points):
            for m, coef in enumerate(p.x.c):
                v = coef * self.D
                assert v.denominator == 1, "denominator does not cover this point"
                out[i, m] = int(v)
            for m, coef in enumerate(p.y.c):
                v = coef * self.D
                assert v.denominator == 1, "denominator does not cover this point"
                out[i, self.dim + m] = int(v)
        return out

    def points(self, rows: np.ndarray) -> List[Point]:
        f = self.field
        inv = Fraction(1, self.D)
        out = []
        for r in rows:
            x = f.element([Fraction(int(r[m])) * inv for m in range(self.dim)])
            y = f.element([Fraction(int(r[self.dim + m])) * inv for m in range(self.dim)])
            out.append(Point(x, y))
        return out

    def floats(self, rows: np.ndarray) -> np.ndarray:
        """(N,2) float coordinates, for radius filters and plotting only."""
        fx = rows[:, : self.dim] @ self.sqrts
        fy = rows[:, self.dim :] @ self.sqrts
        return np.stack([fx, fy], axis=1) / self.D

    # -- exact field arithmetic on int rows --------------------------------
    def _field_square(self, a: np.ndarray) -> np.ndarray:
        """Square of a field element given by integer numerators (N, dim)."""
        n, dim = a.shape
        out = np.zeros((n, dim), dtype=np.int64)
        for i in range(dim):
            ai = a[:, i]
            for j in range(dim):
                out[:, i ^ j] += ai * a[:, j] * self.prod[i & j]
        return out

    def is_unit_vector(self, row: np.ndarray) -> bool:
        """Exact check that a single row has length exactly 1."""
        r = row.reshape(1, -1)
        s = self._field_square(r[:, : self.dim]) + self._field_square(r[:, self.dim :])
        if s[0, 0] != self.D * self.D:
            return False
        return not np.any(s[0, 1:])

    def overflow_headroom(self, rows: np.ndarray) -> float:
        """Ratio of the worst int64 intermediate to the int64 limit.

        The largest intermediate in `_field_square` is max|coeff|^2 * max(prod),
        summed over dim^2 terms for x and y.  Anything at or above 1.0 means the
        integer path is unsafe and must not be trusted.
        """
        if rows.size == 0:
            return 0.0
        mx = int(np.abs(rows).max())
        worst = 2 * self.dim * self.dim * mx * mx * int(self.prod.max())
        return worst / float(np.iinfo(np.int64).max)


def fast_walk(
    basis: IntBasis,
    unit_rows: np.ndarray,
    steps: int,
    radius: float | None = None,
    cap: int | None = None,
    verbose: bool = False,
) -> np.ndarray:
    """All sums of at most `steps` unit vectors, inside `radius`."""
    cur = np.zeros((1, 2 * basis.dim), dtype=np.int64)
    seen = cur
    for s in range(steps):
        cand = (cur[:, None, :] + unit_rows[None, :, :]).reshape(-1, 2 * basis.dim)
        if radius is not None:
            xy = basis.floats(cand)
            keep = (xy[:, 0] ** 2 + xy[:, 1] ** 2) <= radius * radius + 1e-9
            cand = cand[keep]
        if cand.size == 0:
            break
        allrows = np.unique(np.vstack([seen, cand]), axis=0)
        new_count = len(allrows) - len(seen)
        cur = np.unique(cand, axis=0)
        seen = allrows
        if verbose:
            print(f"    step {s+1}: +{new_count} new, {len(seen)} total", flush=True)
        if cap and len(seen) >= cap:
            if verbose:
                print(f"    cap {cap} reached", flush=True)
            break
    return seen


def fast_edges(
    basis: IntBasis, rows: np.ndarray, unit_rows: np.ndarray
) -> List[Tuple[int, int]]:
    """Edges {p, p+u} with both endpoints present, for every u in `unit_rows`.

    Every u is re-verified here as an exact unit vector; a row that is not is
    dropped rather than trusted.
    """
    good = [u for u in unit_rows if basis.is_unit_vector(u)]
    index = {row.tobytes(): i for i, row in enumerate(np.ascontiguousarray(rows))}
    edges = set()
    for u in good:
        shifted = np.ascontiguousarray(rows + u)
        for i, r in enumerate(shifted):
            j = index.get(r.tobytes())
            if j is not None and j != i:
                edges.add((i, j) if i < j else (j, i))
    return sorted(edges)


def fast_graph(
    basis: IntBasis, rows: np.ndarray, unit_rows: np.ndarray, verbose: bool = False
) -> UnitDistanceGraph:
    """Build a `UnitDistanceGraph` (with exact Point coordinates) from int rows."""
    head = basis.overflow_headroom(rows)
    if head >= 1.0:
        raise OverflowError(
            f"int64 path unsafe: intermediates reach {head:.2f} of the limit; "
            "reduce the denominator or the radius"
        )
    edges = fast_edges(basis, rows, unit_rows)
    n = len(rows)
    adj: List[set] = [set() for _ in range(n)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    pts = basis.points(rows)
    g = UnitDistanceGraph(pts, adj)
    if verbose:
        print(f"  {g}  (int64 headroom {head:.1e})", flush=True)
    return g


def fast_edges_complete(
    basis: IntBasis, rows: np.ndarray, verbose: bool = False
) -> List[Tuple[int, int]]:
    """*Every* pair at distance exactly 1 -- no reliance on a supplied vector set.

    Points are bucketed into cells of side 1, so a unit-distance pair always
    falls in the same cell or in one of the 8 around it.  Four offsets plus the
    diagonal cover each unordered cell pair once.  Float distances only nominate
    candidates; each one is then confirmed by exact integer field arithmetic,
    so the returned edges are exactly the unit-distance pairs.
    """
    n = len(rows)
    if n == 0:
        return []
    xy = basis.floats(rows)
    cell = np.floor(xy).astype(np.int64)
    span = int(cell.max() - cell.min()) + 3
    base = cell.min()
    key = (cell[:, 0] - base) * span + (cell[:, 1] - base)

    order = np.argsort(key, kind="stable")
    key_sorted = key[order]
    uniq, starts = np.unique(key_sorted, return_index=True)
    ends = np.append(starts[1:], len(key_sorted))
    block = {int(k): (int(s), int(e)) for k, s, e in zip(uniq, starts, ends)}

    D2 = basis.D * basis.D
    dim = basis.dim
    cand_i: List[np.ndarray] = []
    cand_j: List[np.ndarray] = []

    offsets = [(0, 0), (1, 0), (0, 1), (1, 1), (-1, 1)]
    for k, (s, e) in block.items():
        idx_a = order[s:e]
        pa = xy[idx_a]
        for dx, dy in offsets:
            k2 = k + dx * span + dy
            if k2 not in block:
                continue
            s2, e2 = block[k2]
            idx_b = order[s2:e2]
            pb = xy[idx_b]
            d2 = (pa[:, None, 0] - pb[None, :, 0]) ** 2 + (pa[:, None, 1] - pb[None, :, 1]) ** 2
            hit = np.abs(d2 - 1.0) < 1e-6
            if dx == 0 and dy == 0:
                hit = np.triu(hit, 1)
            ai, bi = np.nonzero(hit)
            if len(ai):
                cand_i.append(idx_a[ai])
                cand_j.append(idx_b[bi])

    if not cand_i:
        return []
    ci = np.concatenate(cand_i)
    cj = np.concatenate(cand_j)
    if verbose:
        print(f"    {len(ci)} candidate pairs from the grid", flush=True)

    # exact confirmation, in blocks to bound memory
    keep_i: List[np.ndarray] = []
    keep_j: List[np.ndarray] = []
    B = 500_000
    for s in range(0, len(ci), B):
        a, b = ci[s : s + B], cj[s : s + B]
        d = rows[a] - rows[b]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        ok = sq[:, 0] == D2
        for m in range(1, dim):
            ok &= sq[:, m] == 0
        keep_i.append(a[ok])
        keep_j.append(b[ok])
    a = np.concatenate(keep_i)
    b = np.concatenate(keep_j)
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)
    pairs = np.unique(np.stack([lo, hi], axis=1), axis=0)
    if verbose:
        print(f"    {len(pairs)} confirmed exactly", flush=True)
    return [(int(u), int(v)) for u, v in pairs]


def fast_graph_complete(
    basis: IntBasis, rows: np.ndarray, verbose: bool = False
) -> UnitDistanceGraph:
    head = basis.overflow_headroom(rows)
    if head >= 1.0:
        raise OverflowError(f"int64 path unsafe: intermediates reach {head:.2f} of the limit")
    edges = fast_edges_complete(basis, rows, verbose=verbose)
    adj: List[set] = [set() for _ in range(len(rows))]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    g = UnitDistanceGraph(basis.points(rows), adj)
    if verbose:
        print(f"  {g}  (int64 headroom {head:.1e})", flush=True)
    return g


class LazyVertices(Sequence):
    """Exact Points materialised only when someone actually asks for one.

    At 10^5+ vertices, eagerly building Fraction coordinates costs more time
    and memory than the whole SAT solve.  Colouring only ever touches the
    adjacency structure; the geometry is needed for candidate distances and
    for the final certificate, and then only for a handful of vertices.
    """

    def __init__(self, basis: "IntBasis", rows: np.ndarray):
        self.basis = basis
        self.rows = rows
        self._cache: dict = {}

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[j] for j in range(*i.indices(len(self)))]
        if i < 0:
            i += len(self)
        p = self._cache.get(i)
        if p is None:
            p = self.basis.points(self.rows[i : i + 1])[0]
            self._cache[i] = p
        return p

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


def fast_graph_lazy(basis: IntBasis, rows: np.ndarray, verbose: bool = False) -> UnitDistanceGraph:
    """`fast_graph_complete` without materialising exact coordinates up front."""
    head = basis.overflow_headroom(rows)
    if head >= 1.0:
        raise OverflowError(f"int64 path unsafe: intermediates reach {head:.2f} of the limit")
    edges = fast_edges_complete(basis, rows, verbose=verbose)
    adj: List[set] = [set() for _ in range(len(rows))]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    g = UnitDistanceGraph(LazyVertices(basis, rows), adj)
    g.rows = rows
    g.basis = basis
    if verbose:
        print(f"  {g}  (int64 headroom {head:.1e})", flush=True)
    return g


def dist2_rational(basis: IntBasis, rows: np.ndarray, pivot: int) -> np.ndarray:
    """Squared distances from `rows[pivot]` to every row, as Fractions over D^2,
    with NaN where the distance is not rational.

    Returns a float array of the rational values (exact small integers/ratios
    survive float64 here because D^2 and the numerators stay well inside 2^53
    for the sizes we search; callers re-derive the exact value for the few
    candidates they keep).
    """
    dim = basis.dim
    d = rows - rows[pivot]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rational = np.ones(len(rows), dtype=bool)
    for m in range(1, dim):
        rational &= sq[:, m] == 0
    out = np.full(len(rows), np.nan)
    out[rational] = sq[rational, 0] / float(basis.D * basis.D)
    return out
