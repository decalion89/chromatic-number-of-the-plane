"""Rotations built from the conflicts they are meant to create.

Every bound in `hn.transversal` is a *same-distance* bound.  They share one
cause: rotating about the pivot keeps each target at its own distance, so
same-distance targets put all their images on one circle, where a point has at
most two neighbours at distance 1.  Degree 2 leaves paths and cycles, Niven
leaves only d^2 = 1/3 among the rational radii, ramification removes the
prime-power circles, and the covering condition caps capacity at 5.

Targets at different distances escape all of it, and it is worth a great deal:
on the same graph at the same k, grouping by distance gives a minimal core of
34 while mixing gives 5.

What does not come free is the *rotations*.  For a single circle they are
forced -- the images must close into a cycle, so the angle is 2 pi t / n.  For
mixed distances nothing forces them, and picking them by hand or by field
order leaves the conflict graph too thin to block: orders 3, 4, 6, 12 reached
a cross-degree of 13 and blocked nothing.

So build them from what they have to do instead.  A rotation rho creates a
conflict between rho(q_a) and q_b exactly when

    |rho(q_a) - q_b|^2 = 1,

and since |rho(q_a)| = |q_a| that is

    <rho(q_a), q_b> = (|q_a|^2 + |q_b|^2 - 1) / 2 =: R.

Writing rho = (c, s) and expanding, <rho(q_a), q_b> = c P + s Q with
P = <q_a, q_b> and Q = the cross product, and P^2 + Q^2 = |q_a|^2 |q_b|^2.  A
line meeting the unit circle: two exact solutions,

    c = (P R +- Q D) / (P^2 + Q^2),   s = (Q R -+ P D) / (P^2 + Q^2),

with D = sqrt(P^2 + Q^2 - R^2).  Every quantity is a field element and the one
square root either lives in the field or names the extension that holds it.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from .geometry import Point, Rotation


def _conflict_terms(p: Point, qa: Point, qb: Point):
    """P, Q, R and the discriminant, for a conflict between qa's image and qb."""
    ax, ay = qa.x - p.x, qa.y - p.y
    bx, by = qb.x - p.x, qb.y - p.y
    da2 = ax * ax + ay * ay
    db2 = bx * bx + by * by
    P = ax * bx + ay * by
    Q = ax * by - ay * bx
    two = P.field.rational(2)
    R = (da2 + db2 - P.field.rational(1)) / two
    disc = P * P + Q * Q - R * R
    return P, Q, R, disc


def conflict_rotations(graph, pivot: int, a: int, b: int) -> List[Rotation]:
    """The rotations about the pivot putting a's image one from b.

    Returns both solutions when the discriminant is a square in the field, and
    nothing when it is not -- the rotation exists in the plane either way, but
    a wider field is needed to name it exactly, and this package never leaves
    exact arithmetic to guess.
    """
    p = graph.vertices[pivot]
    qa, qb = graph.vertices[a], graph.vertices[b]
    P, Q, R, disc = _conflict_terms(p, qa, qb)
    denom = P * P + Q * Q
    if denom == 0:
        return []
    field = P.field
    try:
        D = field.sqrt_element(disc)
    except (ValueError, AttributeError, ZeroDivisionError):
        D = _sqrt_in_field(disc)
        if D is None:
            return []
    out = []
    for sign in (1, -1):
        sg = field.rational(sign)
        c = (P * R + sg * Q * D) / denom
        s = (Q * R - sg * P * D) / denom
        try:
            out.append(Rotation(c, s))
        except ValueError:
            continue
    return out


def _sqrt_in_field(v) -> Optional[object]:
    """A square root of v inside its own field, or None.

    Only the rational case is decided here: a general multiquadratic square
    root is a separate problem, and returning None loses a rotation rather
    than inventing one.
    """
    field = v.field
    if not v.is_rational():
        return None
    q = v.c[0]
    if q < 0:
        return None
    num, den = q.numerator, q.denominator
    try:
        root = field.sqrt(num * den)
    except ValueError:
        return None
    from fractions import Fraction
    return root * field.rational(Fraction(1, den))


def conflict_rotation_set(graph, pivot: int, targets: Sequence[int],
                          limit: int = 400) -> List[Rotation]:
    """Every rotation that creates at least one conflict among the targets.

    These are the copies worth taking.  A rotation creating no conflict adds a
    free choice to every colouring and can only make blocking harder.
    """
    seen = {}
    for a in targets:
        for b in targets:
            for rot in conflict_rotations(graph, pivot, a, b):
                key = (rot.cos, rot.sin)
                if key not in seen:
                    seen[key] = rot
                    if len(seen) >= limit:
                        return list(seen.values())
    return list(seen.values())


def count_cross_transversals(graph, pivot: int, targets: Sequence[int],
                             rotations: Sequence[Rotation], cap: int = 10000) -> int:
    """How many escapes the copies still leave, counted up to `cap`.

    `cross_blocks` answers yes or no, and a no says nothing about how near a
    configuration came.  This counts the independent systems of
    representatives instead -- the colourings that still get away -- so a
    search has something to descend.  Zero is the lemma firing.

    Counting stops at `cap`, since a configuration leaving thousands of
    escapes is not a near miss and the exact figure is worthless.
    """
    from .multispindle import cross_conflict_graph

    targets = list(targets)
    adj = cross_conflict_graph(graph, pivot, targets, rotations)
    m = len(rotations)
    total = 0
    chosen: List[Tuple[int, int]] = []

    def rec(i: int) -> bool:
        """True when the cap is reached and the walk should stop."""
        nonlocal total
        if i == m:
            total += 1
            return total >= cap
        for q in targets:
            key = (i, q)
            if any(c in adj[key] for c in chosen):
                continue
            chosen.append(key)
            if rec(i + 1):
                chosen.pop()
                return True
            chosen.pop()
        return False

    rec(0)
    return total


def in_field_rotations(graph, pivot: int, targets: Sequence[int],
                       extra_orders: Sequence[int] = (3, 4, 6, 8, 12, 24),
                       depth: int = 2, limit: int = 60) -> List[Rotation]:
    """Every rotation this field can name, from the configuration itself.

    The conflict rotations of `conflict_rotation_set` are the ones that do the
    work, but two thirds of them need a square root the multiquadratic field
    does not have -- 52 of 81 pairs at one pivot -- and a tower holding all of
    them would have degree 2^52.  So take what the field does hold: the
    spindle rotation of every distance present, the rotations of finite order,
    and short products of those.
    """
    from .geometry import Rotation as _R
    from .geometry import rotation_joining
    from .multispindle import _rotation_of_order

    p = graph.vertices[pivot]
    field = p.x.field
    seeds = list(conflict_rotation_set(graph, pivot, targets, limit=limit))
    for j in targets:
        d2 = p.dist2(graph.vertices[j])
        if not d2.is_rational():
            continue
        try:
            seeds.append(rotation_joining(d2.c[0], field))
        except (ValueError, ZeroDivisionError):
            continue
    for n in extra_orders:
        r = _rotation_of_order(n, field)
        if r is not None:
            seeds.append(r)
    seen = {(r.cos, r.sin): r for r in seeds}
    cur = list(seen.values())
    for _ in range(depth - 1):
        for a in list(cur):
            for b in list(seen.values()):
                c = a.cos * b.cos - a.sin * b.sin
                s = a.cos * b.sin + a.sin * b.cos
                if (c, s) not in seen:
                    seen[(c, s)] = _R(c, s)
                    if len(seen) >= limit:
                        return list(seen.values())
        cur = list(seen.values())
    return list(seen.values())
