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


# Capacities computed exhaustively by `_capacity` below, which is exponential
# in the cycle length.  Beyond 15 the enumeration stops being cheap; nothing
# here assumes a value it has not computed.
_CAPACITY = {3: 2, 5: 3, 7: 4, 9: 4, 11: 5, 13: 5, 15: 4}

# Past this the enumeration is not affordable: C_33 alone has more independent
# sets than the whole table cost.  Asking beyond it raises rather than guessing,
# because a silent 0 would read as "traps nothing" when it means "not computed".
CAPACITY_LIMIT = 15


def _capacity(n: int) -> int:
    """Largest r for which some r-subset of C_n resists every escape.

    The copies are the n rotations of the cycle, so their image sets are the n
    shifts of the target set T, and a colouring escapes by choosing one image
    per copy with no two adjacent -- an independent set S of C_n meeting every
    shift.  S meets T + k exactly when k lies in S - T, so escape means

        S - T = Z_n   for some independent S,

    and blocking means no independent S covers.  That is the real condition,
    and it is stricter than the independence number suggests.
    """
    from itertools import combinations

    if n < 3 or n % 2 == 0:
        return 0
    indep = [()]
    for size in range(1, n // 2 + 1):
        for S in combinations(range(n), size):
            if all((a - b) % n not in (1, n - 1) for a in S for b in S if a != b):
                indep.append(S)
    full = set(range(n))
    best = 0
    for r in range(1, n + 1):
        for T in combinations(range(n), r):
            if not any({(a - b) % n for a in S for b in T} == full for S in indep):
                best = r
                break
    return best


def trapping_bound(cycle_length: int) -> int:
    """Largest same-distance target set an odd cycle of this length can trap.

    A bipartite component never traps: one side is an independent set meeting
    every copy, so the escape always exists and the bound is 0.

    An earlier version of this returned ``q - (q-1)/2`` here, reading the
    independence number as the answer.  That is an over-estimate, and from
    length 9 upwards a wrong one: the capacities are 2, 3, 4, 4, 5, 5, 4 for
    lengths 3 to 15, against 2, 3, 4, 5, 6, 7, 8.  They do not grow with the
    cycle -- they peak and come back down, because a longer cycle also has
    larger independent sets and covering gets easier faster than trapping does.
    """
    if cycle_length < 3 or cycle_length % 2 == 0:
        return 0
    if cycle_length not in _CAPACITY:
        if cycle_length > CAPACITY_LIMIT:
            raise ValueError(
                f"capacity of C_{cycle_length} is not computed "
                f"(limit {CAPACITY_LIMIT}); the enumeration is exponential")
        _CAPACITY[cycle_length] = _capacity(cycle_length)
    return _CAPACITY[cycle_length]


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
    """Best a field with zeta_n offers, across all the radii it can reach.

    The cycles have length n / gcd(n, t), so every odd divisor of n is a
    reachable cycle length and the best of them is what the field is worth.
    Not the largest: capacity is not monotone in the cycle length, so a field
    with zeta_15 is worth its C_5 (3), not its C_15 (4)... and in fact C_15
    gives 4 while C_5 gives 3, so both have to be checked.
    """
    best = 0
    q = largest_odd_divisor(n)
    for d in range(3, min(q, CAPACITY_LIMIT) + 1, 2):
        if q % d == 0:
            best = max(best, trapping_bound(d))
    return best


def order_for_capacity(r: int, limit: int = 15) -> int:
    """Smallest rotation order whose circle can trap r targets, or 0.

    Bounded by the computed capacities, which stop at cycle length 15 because
    the enumeration is exponential.  Capacity peaks at 5 over that range, so
    anything larger returns 0 -- not "unknown, try harder", but "no cycle up to
    15 does it", and the trend is downwards after 11.
    """
    for n in range(3, limit + 1):
        if trapping_capacity(n) >= r:
            return n
    return 0


# -- Niven: the rational radii, and why there are only four -------------------
#
# The Galois argument above is about which fields admit which rotations.  There
# is a second bound, elementary and field-independent, and it turns out to be
# the binding one for everything this package searched.
#
# A rotation traps only if it closes the images into a cycle, which needs its
# angle commensurable with 2*pi.  And cos t = 1 - 1/(2 d^2) is rational exactly
# when d^2 is.  Niven's theorem says the only rational cosines of rational
# multiples of pi are 0, +-1/2, +-1, so a rational squared distance admits a
# finite-order rotation for just four values:
#
#     d^2 = 1     -> 60 degrees,  order 6
#     d^2 = 1/2   -> 90 degrees,  order 4
#     d^2 = 1/3   -> 120 degrees, order 3     <- the only odd one
#     d^2 = 1/4   -> 180 degrees, order 2
#
# At every other rational d^2 the images form a path, which is bipartite, and
# nothing is trapped at all -- over any field, with any number of copies.
#
# Twelve call sites in this package skip a target whose squared distance is
# irrational.  Every search run here was therefore inside a space capped at 2
# before it started, and no field would have rescued it: the escape needs
# targets at an *irrational* squared distance, namely 1 / (4 sin^2(pi t/n)) for
# an n with an odd divisor of at least 5.

NIVEN_RADII = {
    "1": 6,
    "1/2": 4,
    "1/3": 3,
    "1/4": 2,
}


def niven_order(d2) -> int:
    """Rotation order available at a rational squared distance, or 0.

    `d2` is a Fraction or anything Fraction accepts.  Returns 0 when the angle
    is incommensurable with 2*pi, which is the generic case and means the
    images form a path and cannot be trapped.
    """
    from fractions import Fraction

    return NIVEN_RADII.get(str(Fraction(d2)), 0)


def niven_capacity(d2) -> int:
    """Targets trappable at a rational squared distance: 2 at 1/3, else 0.

    The cycle length here is the rotation order itself, not its largest odd
    divisor.  A radius fixes the adjacency angle, and the images close into a
    cycle whose length is the denominator of that angle over 2*pi in lowest
    terms -- 6 at d^2 = 1, which is even and traps nothing, however many
    triangles an order-6 rotation could reach at some *other* radius.
    `trapping_capacity` answers that other question: what a field offers across
    all its radii.
    """
    n = niven_order(d2)
    return trapping_bound(n) if n else 0
