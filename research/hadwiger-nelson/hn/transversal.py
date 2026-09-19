"""Why the spindle method stalls, and the one number that explains it: 24.

The spindle method blocks a forced disjunction by taking rotated copies of the
configuration.  Each copy contributes "one of my targets carries the pivot's
colour", so a colouring that escapes picks one target image per copy with no
two picked images adjacent -- an independent system of representatives.  The
union is uncolourable exactly when no such choice exists.

Every search in this package grouped targets by their distance from the pivot
before testing them, and that grouping turns out to decide the answer.

**Same distance means one circle.**  Rotating about the pivot preserves the
distance to it, so all images of all same-distance targets lie on a single
circle C.  A point x of C is at distance 1 from at most *two* points of C,
since the unit circle about x and C are distinct circles and two circles meet
in at most two points.  So the unit-distance graph H on the distinct images
has maximum degree at most 2: a disjoint union of paths and cycles, nothing
more.

Measured, because the first version of this file got it wrong.  Counting
degrees in the conflict graph gave 4, 10, 22, 46, 94 as copies were added and
looked like an escape route.  It was multiplicity: at squared distance 1/3,
384 images collapse onto 12 distinct points, and a point hit by 32 different
(copy, target) pairs has 32 times the degree while being one place to stand.
Among distinct points the degree is 2 at every copy count tried.

**Only odd cycles trap.**  A rotation of finite order n spaces each target's
images by 2*pi/n, and adjacency joins images t apart, so H is the circulant
C_n(t) -- cycles of length n / gcd(n, t).  An even cycle or a path is
bipartite, and one side is an independent set meeting every copy.  An odd
cycle C_q has independence number (q-1)/2, so it can trap target sets of size
at most q - (q-1)/2.

**And a multiquadratic field offers exactly one odd cycle.**  A rotation of
order n needs zeta_n inside K(i); for K multiquadratic that Galois group is an
elementary abelian 2-group, which forces (Z/n)* to be one too -- and that
happens precisely when **n divides 24**.  The odd divisors of 24 = 2^3 * 3 are
1 and 3.  The only odd cycle available on any circle, over any multiquadratic
field, with any number of copies, is the triangle.

A triangle has independence number 1, so it traps at most 2 targets.  That is
the ceiling this package spent a day climbing towards from below: a disjunction
over 34 targets narrowed to 11 and stopped, needing 2.  It was not a search
that failed.  Over ``Q(sqrt3, sqrt5, sqrt7, sqrt11)`` -- de Grey's field, and
the field of every configuration examined here -- no search could have arrived.

**The escape is the same sentence read forwards.**  Leave n | 24.  A rotation
of order 5 lives in ``Q(zeta_5)``, whose Galois group is cyclic of order 4 and
therefore not multiquadratic; it makes H a union of pentagons, independence
number 2, trapping up to 3 targets.  Order 7 gives heptagons and 4.  The
ceiling is not a fact about the plane, it is a fact about which fields the
search was willing to build points in.
"""
from __future__ import annotations

from math import gcd
from typing import Iterable, Sequence

# A rotation of order n is available over a multiquadratic field exactly when
# (Z/n)* is an elementary abelian 2-group, i.e. when n divides 24.
MULTIQUADRATIC_ORDER_BOUND = 24


def multiquadratic_orders() -> list:
    """Rotation orders a multiquadratic field admits: the divisors of 24."""
    n = MULTIQUADRATIC_ORDER_BOUND
    return [d for d in range(1, n + 1) if n % d == 0]


def cycle_lengths(order: int) -> list:
    """Lengths of the cycles H can have, for a rotation of the given order.

    Adjacency joins images t steps apart, so H is the circulant C_n(t) and its
    cycles have length n / gcd(n, t).
    """
    if order < 1:
        raise ValueError("order must be positive")
    return sorted({order // gcd(order, t) for t in range(1, order + 1)})


def trapping_bound(cycle_length: int) -> int:
    """Largest same-distance target set an odd cycle of this length can trap.

    A bipartite component never traps: one side is an independent set meeting
    every copy, so the escape always exists and the bound is 0.  An odd cycle
    C_q has independence number (q-1)/2, leaving q - (q-1)/2 targets trappable.
    """
    if cycle_length < 3 or cycle_length % 2 == 0:
        return 0
    return cycle_length - (cycle_length - 1) // 2


def same_distance_ceiling(orders: Iterable[int]) -> int:
    """The most targets any of these rotation orders could ever trap.

    Over the multiquadratic orders this returns 2, for any number of copies.
    """
    best = 0
    for n in orders:
        for q in cycle_lengths(n):
            best = max(best, trapping_bound(q))
    return best


def count_circles(graph, pivot: int, targets: Iterable[int]) -> int:
    """How many distinct distances from the pivot the targets occupy.

    One means the ceiling above applies.  More than one means the images are
    spread over several circles and this analysis says nothing.
    """
    p = graph.vertices[pivot]
    return len({p.dist2(graph.vertices[j]) for j in targets})


def diagnose(graph, pivot: int, targets: Sequence[int], orders=None) -> dict:
    """Say whether a forced disjunction is worth trying to block by rotation."""
    targets = list(targets)
    orders = list(orders) if orders is not None else multiquadratic_orders()
    circles = count_circles(graph, pivot, targets)
    ceiling = same_distance_ceiling(orders)
    single = circles == 1
    return {
        "targets": len(targets),
        "circles": circles,
        "single_circle": single,
        "ceiling": ceiling,
        # Only a single-circle disjunction is ruled out; several circles put
        # the images outside the reach of this argument, in either direction.
        "hopeless": single and len(targets) > ceiling,
    }


# -- the magic circles --------------------------------------------------------
#
# The analysis above assumed H was a cycle, and that assumption pins the radius.
# Adjacency on the circle joins images whose angular separation is
# beta = 2*arcsin(1/(2d)), and H is the circulant C_n(t) only when beta is
# exactly t steps of the rotation, beta = 2*pi*t/n.  Then
#
#     cos beta = 1 - 1/(2 d^2) = cos(2 pi t / n)   =>   d = 1/(2 sin(pi t / n)),
#
# so each rotation order comes with its own radius, and targets anywhere else
# see a path, which is bipartite, and cannot be trapped at all.
#
# n = 3 gives d = 1/sqrt(3): the circle where two points are adjacent exactly
# when they are 120 degrees apart, the one every spindle argument in the
# literature uses, and -- since n | 24 -- the only one a multiquadratic field
# has.  The rest of the hierarchy has never been available to look at.


def largest_odd_divisor(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


def magic_radius(n: int, t: int = 1) -> float:
    """Where targets must sit for an order-n rotation to close them into cycles.

    Exactly, d^2 = 1 / (2 - zeta^t - zeta^-t), an element of Q(zeta_n).  The
    float is for reading; nothing in a decision path should use it.
    """
    from math import pi, sin

    if n < 3 or not 0 < t < n:
        raise ValueError("need n >= 3 and 0 < t < n")
    return 1.0 / (2.0 * sin(pi * t / n))


def trapping_capacity(n: int) -> int:
    """Targets an order-n rotation can trap, once they sit on its circle.

    The cycles have length n / gcd(n, t); the odd ones are what trap, and the
    longest available is the largest odd divisor q of n, holding (q + 1) / 2.
    """
    q = largest_odd_divisor(n)
    return trapping_bound(q) if q >= 3 else 0


def order_for_capacity(r: int, limit: int = 200) -> int:
    """Smallest rotation order whose circle can trap r targets, or 0.

    The narrowing runs here stalled at 11 same-distance targets while the
    multiquadratic capacity is 2.  Eleven needs an odd cycle of length 21.
    """
    for n in range(3, limit + 1):
        if trapping_capacity(n) >= r:
            return n
    return 0
