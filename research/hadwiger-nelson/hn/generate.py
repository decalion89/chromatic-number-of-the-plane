"""Vertex-set generators.

Everything we build lives in the ring  Z[omega, sigma],
    omega = exp(i*pi/3)      = (1 + i*sqrt(3))/2   -- the triangular lattice
    sigma = exp(i*arccos(5/6)) = (5 + i*sqrt(11))/6 -- the Moser spindle hinge

Both have modulus 1, so every product omega^a sigma^m is a *unit vector*, and
a sum of n such vectors is a point reachable by an n-step walk of unit steps.
The 1581-vertex graph of de Grey (2018) and the 509-vertex graph of Parts
(2020) both live inside this ring; generating it systematically is how we get
a superset to search.
"""

from __future__ import annotations

import math
from typing import Iterable, List, Sequence, Set

from .field import Field, QSQRT3_11
from .geometry import Point, Rotation, ROT60, SPINDLE, eisenstein, origin

__all__ = [
    "unit_vectors",
    "unit_vectors_multi",
    "lattice_rotations",
    "walk_ball",
    "walk_ball_from",
    "hex_ball",
    "minkowski",
    "spindled",
    "trim_to_radius",
]


def unit_vectors(m_max: int = 1, field: Field = QSQRT3_11) -> List[Point]:
    """{ omega^a sigma^m : 0 <= a < 6, |m| <= m_max }, as points on the unit circle."""
    out: List[Point] = []
    seen: Set[Point] = set()
    e1 = Point(field.one(), field.zero())
    for m in range(-m_max, m_max + 1):
        base = (SPINDLE**m)(e1)
        for a in range(6):
            p = (ROT60**a)(base)
            if p not in seen:
                seen.add(p)
                out.append(p)
    return out


def walk_ball(
    steps: int,
    m_max: int = 1,
    radius: float | None = None,
    field: Field = QSQRT3_11,
    verbose: bool = False,
) -> List[Point]:
    """All points reachable from the origin in at most `steps` unit steps drawn
    from `unit_vectors(m_max)`, optionally kept inside a disk of `radius`."""
    U = unit_vectors(m_max, field)
    frontier = {origin(field)}
    seen = set(frontier)
    for s in range(steps):
        nxt = set()
        for p in frontier:
            for u in U:
                q = p + u
                if q in seen:
                    continue
                if radius is not None and (q.fx * q.fx + q.fy * q.fy) > radius * radius + 1e-9:
                    continue
                nxt.add(q)
        seen |= nxt
        frontier = nxt
        if verbose:
            print(f"    step {s+1}: +{len(nxt)} new, {len(seen)} total")
        if not frontier:
            break
    return list(seen)


def hex_ball(r: int, field: Field = QSQRT3_11) -> List[Point]:
    """The triangular-lattice ball of graph-radius r (1, 7, 19, 37, ... points)."""
    out = []
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            if abs(a + b) <= r:
                out.append(eisenstein(a, b, field))
    return out


def minkowski(A: Sequence[Point], B: Sequence[Point], radius: float | None = None) -> List[Point]:
    """{ a + b }, deduplicated."""
    seen: Set[Point] = set()
    for a in A:
        for b in B:
            p = a + b
            if radius is not None and (p.fx * p.fx + p.fy * p.fy) > radius * radius + 1e-9:
                continue
            seen.add(p)
    return list(seen)


def spindled(points: Sequence[Point], pivot: Point, rot: Rotation = SPINDLE) -> List[Point]:
    """S union rot(S) about `pivot` -- the move that builds the Moser spindle."""
    turn = rot.about(pivot)
    seen: Set[Point] = set(points)
    for p in points:
        seen.add(turn(p))
    return list(seen)


def trim_to_radius(points: Iterable[Point], radius: float, center: Point | None = None) -> List[Point]:
    cx, cy = (0.0, 0.0) if center is None else (center.fx, center.fy)
    r2 = radius * radius + 1e-9
    return [p for p in points if (p.fx - cx) ** 2 + (p.fy - cy) ** 2 <= r2]


# --- richer rotation families -------------------------------------------

def lattice_rotations(max_d2: int = 100, field: Field = QSQRT3_11) -> List[tuple]:
    """Every rotation arccos(1 - 1/(2 d2)) with d2 an Eisenstein norm whose sine
    stays inside `field`.

    Rotating a point at distance sqrt(d2) from the pivot by this angle sends it
    to distance exactly 1 from where it started -- the spindling move, once per
    realisable lattice distance.  For the centred hexagonal numbers
    d2 = 3k^2+3k+1 (1, 7, 19, 37, 61, 91, ...) the sine is (2k+1)sqrt(3)/(2 d2),
    so that whole family is available without extending the field.

    Returns [(d2, Rotation), ...].
    """
    from .geometry import rotation_joining

    norms = sorted({a * a + a * b + b * b for a in range(-15, 16) for b in range(-15, 16)} - {0})
    out = []
    for d2 in norms:
        if d2 > max_d2:
            break
        try:
            out.append((d2, rotation_joining(d2, field)))
        except ValueError:
            continue  # the sine would need a bigger field
    return out


def unit_vectors_multi(
    rotations: Sequence[Rotation],
    depth: int = 1,
    field: Field = QSQRT3_11,
    cap: int | None = None,
) -> List[Point]:
    """Unit vectors reachable by composing at most `depth` of `rotations`
    (and their inverses) starting from (1, 0), closed under 60-degree turns."""
    e1 = Point(field.one(), field.zero())
    seen = {e1}
    frontier = [e1]
    gens = []
    for r in rotations:
        gens.append(r)
        gens.append(r.inverse())
    for _ in range(depth):
        nxt = []
        for p in frontier:
            for r in gens:
                q = r(p)
                if q not in seen:
                    seen.add(q)
                    nxt.append(q)
        frontier = nxt
        if cap and len(seen) > cap:
            break
    out = set()
    for p in seen:
        q = p
        for _ in range(6):
            out.add(q)
            q = ROT60(q)
    return list(out)


def walk_ball_from(
    unit_set: Sequence[Point],
    steps: int,
    radius: float | None = None,
    field: Field = QSQRT3_11,
    cap: int | None = None,
    verbose: bool = False,
) -> List[Point]:
    """`walk_ball`, but over an explicitly supplied set of unit steps."""
    frontier = {origin(field)}
    seen = set(frontier)
    for s in range(steps):
        nxt = set()
        for p in frontier:
            for u in unit_set:
                q = p + u
                if q in seen:
                    continue
                if radius is not None and (q.fx * q.fx + q.fy * q.fy) > radius * radius + 1e-9:
                    continue
                nxt.add(q)
        seen |= nxt
        frontier = nxt
        if verbose:
            print(f"    step {s+1}: +{len(nxt)} new, {len(seen)} total", flush=True)
        if cap and len(seen) > cap:
            break
        if not frontier:
            break
    return list(seen)
