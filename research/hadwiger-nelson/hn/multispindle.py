"""The spindle argument in its general form.

Both spindles used so far are special cases of one lemma, and stating the
lemma properly shows that the usual "three copies at most" ceiling is an
artefact of asking for more than the argument needs.

    Multi-copy spindle.  Let G be a unit-distance graph, p a vertex of G and
    Q a set of vertices.  Suppose every proper k-colouring of G gives p the
    same colour as at least one member of Q.  Let rho_1..rho_m be rotations
    about p, and let W be the union of the rotated copies rho_i(G).

    For each q in Q form the conflict graph H_q on {1..m}:

        i ~ j   iff   |rho_i(q) - rho_j(q)| = 1.

    If no function f : {1..m} -> Q has every f-class independent in the
    corresponding H -- that is, no f avoids all i != j with f(i) = f(j) and
    i ~_{H_f(i)} j -- then W has no proper k-colouring at all.

    Proof.  Each rho_i(G) is congruent to G, so any proper colouring c of W
    restricted to it is a proper colouring of G, giving some f(i) in Q with
    c(rho_i(f(i))) = c(rho_i(p)) = c(p), the last step because a rotation
    about p fixes p.  By hypothesis f is not an independent transversal, so
    some i != j have f(i) = f(j) = q and rho_i(q) ~ rho_j(q).  Both carry
    c(p) and they are adjacent in W.  Contradiction.

Where the ceiling actually comes from.  The familiar arguments take every
H_q to be complete, which forces m > |Q| rotated images to be pairwise at
distance 1.  Two such points fit on any circle of radius >= 1/2; three fit
only on the circumcircle of a unit equilateral triangle, radius 1/sqrt(3) --
equivalently, the spindle angle there is exactly 120 degrees, and three of
them close the circle.  A fourth is impossible, so complete conflict graphs
cap the argument at |Q| <= 2.

The lemma does not ask for complete conflict graphs.  It asks only that the
H_q *jointly* block every assignment.  That is a much weaker demand, and it
opens two doors the complete-graph version cannot reach:

  * A base rotation of odd order n gives H_q = an n-cycle rather than a
    triangle.  Odd cycles are not 2-colourable, so |Q| = 2 is blocked with
    n copies for any odd n, not just n = 3.

  * Targets at *different* radii give *different* conflict graphs, since the
    step length depends on the radius.  Blocking then means that {1..m}
    cannot be covered by one independent set from each H_q -- a genuine
    combinatorial condition that can hold for |Q| >= 3, which no single
    complete graph on a circle can ever achieve.

The condition is a small finite check (|Q|^m assignments, with m and |Q| in
single digits), so it costs nothing to test once a forced target set is in
hand.  Forced sets of size 3 or more, which the earlier searches discarded
as unusable, become usable here.
"""

from __future__ import annotations

from itertools import combinations_with_replacement, product
from typing import Dict, List, Optional, Sequence, Tuple

from .field import Field
from .geometry import Point, Rotation
from .graph import UnitDistanceGraph, build_graph

__all__ = [
    "MULTIQUADRATIC_ORDERS",
    "available_rotation_orders",
    "spindle_catalogue",
    "squared_distance_for_step",
    "conflict_graphs",
    "independent_transversal",
    "blocks_all_assignments",
    "multi_spindle_union",
    "rotation_powers",
]


def rotation_powers(rot: Rotation, m: int) -> List[Rotation]:
    """{rot^0, ..., rot^(m-1)} -- the natural family of copies about a pivot."""
    out, cur = [], Rotation(rot.field.one(), rot.field.zero(), check=False)
    for _ in range(m):
        out.append(cur)
        cur = cur * rot
    return out


def conflict_graphs(
    graph: UnitDistanceGraph,
    pivot: int,
    targets: Sequence[int],
    rotations: Sequence[Rotation],
) -> Dict[int, List[set]]:
    """H_q for each target q, computed with exact arithmetic."""
    p = graph.vertices[pivot]
    m = len(rotations)
    out: Dict[int, List[set]] = {}
    for q_idx in targets:
        q = graph.vertices[q_idx]
        images = [r.about(p)(q) for r in rotations]
        adj = [set() for _ in range(m)]
        for i in range(m):
            for j in range(i + 1, m):
                if images[i].is_unit_apart(images[j]):
                    adj[i].add(j)
                    adj[j].add(i)
        out[q_idx] = adj
    return out


def independent_transversal(
    m: int, targets: Sequence[int], H: Dict[int, List[set]]
) -> Optional[Tuple[int, ...]]:
    """An assignment of copies to targets with every class independent, or None.

    Returning None is the good case: it means the copies cannot be shared out
    among the targets without some pair colliding, which is exactly the
    hypothesis the lemma needs.
    """
    targets = list(targets)
    # depth-first with pruning; |targets|^m is small but pruning makes it trivial
    assign: List[int] = []

    def ok(i: int, q: int) -> bool:
        adj = H[q][i]
        return not any(assign[j] == q and j in adj for j in range(i))

    def rec(i: int) -> Optional[Tuple[int, ...]]:
        if i == m:
            return tuple(assign)
        for q in targets:
            if ok(i, q):
                assign.append(q)
                r = rec(i + 1)
                if r is not None:
                    return r
                assign.pop()
        return None

    return rec(0)


def blocks_all_assignments(
    graph: UnitDistanceGraph,
    pivot: int,
    targets: Sequence[int],
    rotations: Sequence[Rotation],
) -> bool:
    """True when this family of rotations makes the lemma fire."""
    H = conflict_graphs(graph, pivot, targets, rotations)
    return independent_transversal(len(rotations), targets, H) is None


def multi_spindle_union(
    graph: UnitDistanceGraph, pivot: int, rotations: Sequence[Rotation]
) -> UnitDistanceGraph:
    """W = union of rho_i(G) over the given rotations about `pivot`."""
    p = graph.vertices[pivot]
    pts, seen = [], set()
    for r in rotations:
        turn = r.about(p)
        for v in graph.vertices:
            q = turn(v)
            if q not in seen:
                seen.add(q)
                pts.append(q)
    return build_graph(pts)


# --- which spindles a given field can actually support ---------------------
#
# A rotation of order n needs both cos(2pi/n) and sin(2pi/n) in the field, so
# zeta_n = cos + i sin lies in F(i).  When F is multiquadratic so is F(i), and
# a subfield of a multiquadratic field is multiquadratic, so Gal(Q(zeta_n)/Q)
# = (Z/n)* must have exponent at most 2.  That holds exactly when n divides 24.
#
#     Realisable orders:  1, 2, 3, 4, 6, 8, 12, 24.
#     The only odd one above 1 is 3.
#
# That order bound is classical, not a discovery here: (Z/n)* has exponent 2
# exactly when n | 24 is a standard fact, and the rest is bookkeeping.  What it
# buys is a precise statement of which rotations are on the table at all.
#
# Each conflict graph has a single connection element, so it is a union of
# cycles of that element's order (or of paths, for infinite order) and is
# bipartite unless that order is odd.  Since 3 is the only odd order available,
# the only *non-bipartite* conflict graph over a multiquadratic field is the
# triangle -- which is exactly the d^2 = 1/3 case, recovered from arithmetic
# rather than from the picture of a triangle inscribed in a circle.
#
# Two cautions, both of which cost earlier versions of this file a wrong claim:
#
#   * Non-bipartiteness is NOT required to block two targets.  When the two
#     conflict graphs differ, a pair of bipartite graphs can be jointly
#     unblockable: order 4 with steps (1, 2) does it.  What χ(H) > r governs is
#     the case where every H_q is the *same* graph.  With distinct radii the
#     question is a joint covering problem, and the answer is not read off any
#     single chromatic number.
#
#   * Three targets have NOT been shown impossible.  `spindle_catalogue`
#     searches families of the form {rho^0 .. rho^(n-1)} for a single rho of
#     order n, and over those it finds no three-target configuration.  The
#     lemma allows arbitrary finite families, including a torsion rotation
#     mixed with an infinite-order one, and that space has not been settled
#     here.  Treat "two targets" as the reach of the catalogue, not as a
#     theorem.
#
# Leaving multiquadratic fields is still the natural direction for more odd
# orders -- Q(sin 2pi/7) carries order-7 rotations -- but that is a lead, not
# a result.

import math
from fractions import Fraction

MULTIQUADRATIC_ORDERS = (1, 2, 3, 4, 6, 8, 12, 24)


def _rotation_of_order(n: int, field: Field) -> Optional[Rotation]:
    """The rotation by 2*pi/n, if this field contains it.

    Built one case at a time: asking a field for a radical it does not have
    raises, and that must disqualify only the order that needed it.
    """
    half = Fraction(1, 2)
    quarter = Fraction(1, 4)
    try:
        if n == 1:
            c, sn = field.one(), field.zero()
        elif n == 2:
            c, sn = field.rational(-1), field.zero()
        elif n == 3:
            c, sn = field.rational(-half), field.sqrt(3) * field.rational(half)
        elif n == 4:
            c, sn = field.zero(), field.one()
        elif n == 6:
            c, sn = field.rational(half), field.sqrt(3) * field.rational(half)
        elif n == 8:
            r2 = field.sqrt(2) * field.rational(half)
            c, sn = r2, r2
        elif n == 12:
            c, sn = field.sqrt(3) * field.rational(half), field.rational(half)
        elif n == 24:
            r6, r2 = field.sqrt(6), field.sqrt(2)
            c = (r6 + r2) * field.rational(quarter)
            sn = (r6 - r2) * field.rational(quarter)
        else:
            return None
    except ValueError:
        return None          # this field lacks the radical that order needs
    try:
        return Rotation(c, sn)
    except ValueError:
        return None


def available_rotation_orders(field: Field) -> List[int]:
    """The finite rotation orders this field realises exactly."""
    out = []
    for n in MULTIQUADRATIC_ORDERS:
        try:
            r = _rotation_of_order(n, field)
        except ValueError:
            continue
        if r is None:
            continue
        p = Rotation(field.one(), field.zero(), check=False)
        for _ in range(n):
            p = p * r
        if p.cos == 1 and p.sin.is_zero():
            out.append(n)
    return out


def _circulant(n: int, s: int) -> List[set]:
    return [{(i + s) % n, (i - s) % n} - {i} for i in range(n)]


def _coverable(n: int, steps: Sequence[int]) -> bool:
    H = [_circulant(n, s) for s in steps]
    r = len(steps)
    a: List[int] = []

    def rec(i: int) -> bool:
        if i == n:
            return True
        for t in range(r):
            if not any(a[j] == t and j in H[t][i] for j in range(i)):
                a.append(t)
                if rec(i + 1):
                    return True
                a.pop()
        return False

    return rec(0)


def squared_distance_for_step(n: int, s: int, field: Field):
    """The d^2 whose spindle angle is s turns of 2*pi/n, as an exact field element.

    From cos(angle) = 1 - 1/(2 d^2).  d^2 = 1 means the target sits at distance
    one from the pivot, hence adjacent to it and never monochromatic with it, so
    such a step is useless and is reported as None.
    """
    rot = _rotation_of_order(n, field)
    if rot is None:
        return None
    p = Rotation(field.one(), field.zero(), check=False)
    for _ in range(s):
        p = p * rot
    one_minus = field.one() - p.cos
    if one_minus.is_zero():
        return None
    d2 = (one_minus * field.rational(2)).inverse()
    if d2 == 1:
        return None
    return d2


def spindle_catalogue(field: Field, max_targets: int = 3) -> List[dict]:
    """Every multi-copy spindle this field supports.

    Each entry is {n, steps, d2, rotations}: use those n rotated copies about
    the pivot, with a forced target set holding one member at each listed
    squared distance.

    The rotations are carried explicitly because a catalogued d^2 need not be
    rational -- 2 + sqrt(3) occurs at order 12 -- and `rotation_joining` only
    accepts rational squared distances.  Here the rotation is known by
    construction as a power of the order-n rotation, so nothing has to take a
    square root inside the field to recover it.
    """
    out = []
    for n in available_rotation_orders(field):
        if n < 3:
            continue
        for r in range(1, max_targets + 1):
            for steps in combinations_with_replacement(range(1, n // 2 + 1), r):
                d2s = [squared_distance_for_step(n, s, field) for s in steps]
                if any(d is None for d in d2s):
                    continue
                if not _coverable(n, steps):
                    base = _rotation_of_order(n, field)
                    out.append({
                        "n": n,
                        "steps": steps,
                        "d2": d2s,
                        "rotations": rotation_powers(base, n),
                    })
    return out
