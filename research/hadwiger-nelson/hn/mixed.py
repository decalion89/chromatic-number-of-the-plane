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

from fractions import Fraction
from typing import List, Optional, Sequence, Tuple

from .field import Field
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


def _rational_sqrt(field, q):
    """sqrt of a non-negative rational, inside the field, or None."""
    from fractions import Fraction

    if q < 0:
        return None
    if q == 0:
        return field.zero()
    num, den = q.numerator, q.denominator
    try:
        root = field.sqrt(num * den)
    except ValueError:
        return None
    return root * field.rational(Fraction(1, den))


def _sqrt_in_field(v) -> Optional[object]:
    """A square root of v inside its own field, or None.

    The rational case is easy.  The two-term case is the classical denesting

        sqrt(a + b sqrt(d)) = sqrt(x) + sqrt(y),
        x + y = a,  4 x y = b^2 d,   so   x, y = (a +- sqrt(a^2 - b^2 d)) / 2,

    which lands inside a multiquadratic field exactly when a^2 - b^2 d is a
    rational square and the squarefree parts of x and y are already generators.
    A negative a^2 - b^2 d means the element is not totally positive, so its
    square root is not in any real multiquadratic field at all -- that is
    arithmetic refusing, not this function giving up.

    It matters more than it looks.  Restricted to the rational case this
    returned nothing at all for a core of three on Sa, whose conflict
    discriminants are all of the form a + b sqrt(33); three of the six denest,
    and they are exactly the cross-target pairs that mixed-distance blocking
    needs.  The whole rotation set was being thrown away by a weak square root.
    """
    field = v.field
    terms = [(i, c) for i, c in enumerate(v.c) if c]
    if not terms:
        return field.zero()
    if len(terms) == 1 and terms[0][0] == 0:
        return _rational_sqrt(field, terms[0][1])
    if len(terms) != 2 or terms[0][0] != 0:
        return None                       # more than two terms: not handled
    a = terms[0][1]
    mask, b = terms[1]
    d = 1
    for i, g in enumerate(field.gens):
        if mask >> i & 1:
            d *= g
    inner = _rational_sqrt(field, a * a - b * b * d)
    if inner is None or not inner.is_rational():
        return None
    sgn = inner.c[0]
    from fractions import Fraction
    x = (a + sgn) / 2
    y = (a - sgn) / 2
    rx, ry = _rational_sqrt(field, x), _rational_sqrt(field, y)
    if rx is None or ry is None:
        return None
    root = rx + ry if b > 0 else rx - ry
    return root if root * root == v else None


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


def search_block(graph, pivot: int, targets: Sequence[int],
                 candidates: Sequence[Rotation], cap: int = 40000,
                 nodes: int = 20000, max_depth: int = 14, width: int = 4):
    """Search for a set of copies leaving no escape at all.

    Adding copies wholesale makes things worse, not better: seven conflict
    rotations left 432 escapes on a core of three, and closing them under
    products to nineteen pushed it past 200000. Every copy multiplies the
    choices by |targets| while pruning only in proportion to its own conflict
    degree, so the count is not monotone and the set has to be chosen.

    An escape is a choice of one target per copy whose images are pairwise
    non-adjacent. The surviving escapes are kept explicitly, and a new copy
    extends each of them by its own targets that clash with nothing already
    chosen.

    The count rises before it falls -- with no copies there is exactly one
    escape, the empty assignment, and the first copy always takes it to
    |targets| -- so refusing any copy that increases it stalls at the root,
    which is what a first attempt did. Every candidate under the cap is tried
    instead, cheapest first, keeping the best `width` at each level and
    backtracking. Reaching zero is the lemma firing.
    """
    from .multispindle import cross_conflict_graph

    targets = list(targets)
    cands = list(candidates)
    adj = cross_conflict_graph(graph, pivot, targets, cands)
    budget = [nodes]
    best = [None, None]

    def extend(esc, i):
        out = []
        for e in esc:
            for q in targets:
                if any(c in adj[(i, q)] for c in e):
                    continue
                out.append(e + ((i, q),))
                if len(out) > cap:
                    return None            # runaway branch, not worth keeping
        return out

    def rec(used, esc, depth):
        if not esc:
            best[0], best[1] = list(used), 0
            return True
        if budget[0] <= 0 or depth == 0:
            return False
        if used and (best[1] is None or len(esc) < best[1]):
            best[0], best[1] = list(used), len(esc)
        options = []
        for i in range(len(cands)):
            if i in used:
                continue
            budget[0] -= 1
            if budget[0] <= 0:
                break
            nxt = extend(esc, i)
            if nxt is not None:
                options.append((len(nxt), i, nxt))
        options.sort(key=lambda t: t[0])
        for _, i, nxt in options[:width]:
            if rec(used + [i], nxt, depth - 1):
                return True
        return False

    rec([], [()], max_depth)
    return ([cands[i] for i in (best[0] or [])], best[1])


def unit_circle_intersections(u: Point, v: Point) -> List[Point]:
    """The points at distance exactly 1 from both u and v.

    On the perpendicular bisector, h from the midpoint, with
    h^2 = 1 - D/4 for D = |u - v|^2.  So

        x = (u + v)/2 +- sqrt((4 - D) / (4 D)) * (-(dy), dx),   d = v - u,

    and the square root goes through the same denesting as everything else --
    present when the field holds it, and honestly absent otherwise.
    """
    d = v - u
    D = d.x * d.x + d.y * d.y
    field = D.field
    four = field.rational(4)
    if D == 0:
        return []
    t = (four - D) / (four * D)
    root = _sqrt_in_field(t)
    if root is None:
        return []
    half = field.rational(Fraction(1, 2))
    mx, my = (u.x + v.x) * half, (u.y + v.y) * half
    out = []
    for sign in (1, -1):
        sg = field.rational(sign)
        x = Point(mx - sg * root * d.y, my + sg * root * d.x)
        if x.dist2(u) == 1 and x.dist2(v) == 1:
            out.append(x)
    return out


def deep_holes(graph, min_degree: int = 6, limit: int = 4000):
    """Points of the plane with many graph vertices exactly one away.

    Every pivot tried in this package was already a vertex, which is a
    restriction the argument never asked for: the pivot is a point whose
    colour is being constrained, and any point of the plane will do. The ones
    worth adding are those with the most neighbours, since forcing comes from
    how tightly a pivot's own neighbourhood is pinned.

    Candidates are the intersections of unit circles about pairs of vertices,
    which is where a point can have two neighbours at all; the count is then
    exact against the whole vertex set.
    """
    from collections import defaultdict

    pts = list(graph.vertices)
    zs = [(complex(float(p.x), float(p.y)) if hasattr(p.x, "__float__")
           else None) for p in pts]
    seen = {}
    for i, u in enumerate(pts):
        for j in range(i + 1, len(pts)):
            v = pts[j]
            D = u.dist2(v)
            if not D.is_rational():
                continue
            if not (0 < D.c[0] <= 4):
                continue
            for x in unit_circle_intersections(u, v):
                if x not in seen:
                    seen[x] = None
                    if len(seen) >= limit:
                        break
            if len(seen) >= limit:
                break
        if len(seen) >= limit:
            break
    out = []
    for x in seen:
        deg = sum(1 for q in pts if x.is_unit_apart(q))
        if deg >= min_degree:
            out.append((deg, x))
    out.sort(key=lambda t: -t[0])
    return out


class Reflection:
    """Reflection in a line through the origin, at half-angle (c, s).

    The spindle argument needs its copies to fix the pivot and preserve
    distances, and rotations are not the only isometries that do. A reflection
    in any line through the pivot fixes it, so its images carry the pivot's
    colour exactly as a rotation's do, and the lemma applies unchanged.

    They are also free. A rotation by an angle that creates a conflict needs a
    square root, and two thirds of those turned out to be outside the field.
    Reflecting in the line at angle t needs cos 2t and sin 2t, which is the
    same data as a rotation by 2t -- so every rotation already in hand gives a
    reflection for nothing, and the composition of two reflections is a
    rotation, so the group they generate is larger than the rotations alone.
    """

    __slots__ = ("cos", "sin")

    def __init__(self, cos, sin):
        self.cos = cos
        self.sin = sin

    @property
    def field(self):
        return self.cos.field

    def __call__(self, p: Point) -> Point:
        # (x, y) -> (x cos + y sin, x sin - y cos), the reflection in the line
        # at half the angle of (cos, sin)
        return Point(p.x * self.cos + p.y * self.sin,
                     p.x * self.sin - p.y * self.cos)

    def about(self, pivot: Point):
        def reflect(p: Point, _r=self, _c=pivot) -> Point:
            return _r(p - _c) + _c

        return reflect


def reflections_from(rotations: Sequence) -> List[Reflection]:
    """One reflection per rotation, reusing its cos and sin exactly."""
    return [Reflection(r.cos, r.sin) for r in rotations]


def conflict_reflections(graph, pivot: int, a: int, b: int) -> List[Reflection]:
    """The reflections about the pivot putting a's image one from b.

    The same derivation as for rotations, with the reflection's expansion in
    place of the rotation's:

        <sigma(a), b> = c (a_x b_x - a_y b_y) + s (a_y b_x + a_x b_y),

    so P and Q change but P^2 + Q^2 is still |a|^2 |b|^2, and R is unchanged.
    **The discriminant is therefore identical.** Whenever a conflict rotation
    can be named in the field, a conflict reflection can be named too, from the
    same square root -- the copies double for nothing, and no new radical is
    needed anywhere.
    """
    p = graph.vertices[pivot]
    qa, qb = graph.vertices[a], graph.vertices[b]
    ax, ay = qa.x - p.x, qa.y - p.y
    bx, by = qb.x - p.x, qb.y - p.y
    da2 = ax * ax + ay * ay
    db2 = bx * bx + by * by
    field = ax.field
    P = ax * bx - ay * by
    Q = ay * bx + ax * by
    R = (da2 + db2 - field.rational(1)) / field.rational(2)
    denom = P * P + Q * Q
    if denom == 0:
        return []
    D = _sqrt_in_field(denom - R * R)
    if D is None:
        return []
    out = []
    for sign in (1, -1):
        sg = field.rational(sign)
        c = (P * R + sg * Q * D) / denom
        s = (Q * R - sg * P * D) / denom
        if c * c + s * s != 1:
            continue
        r = Reflection(c, s)
        if r.about(p)(qa).is_unit_apart(qb):
            out.append(r)
    return out


def conflict_isometries(graph, pivot: int, targets: Sequence[int],
                        limit: int = 400) -> List:
    """Every rotation *and* reflection creating a conflict among the targets.

    Reflections were left out of this package entirely, and they cost nothing:
    they fix the pivot, so the lemma applies to them unchanged, and they share
    their square root with the rotations.
    """
    seen = {}
    for a in targets:
        for b in targets:
            for r in conflict_rotations(graph, pivot, a, b):
                seen[("rot", r.cos, r.sin)] = r
            for r in conflict_reflections(graph, pivot, a, b):
                seen[("ref", r.cos, r.sin)] = r
            if len(seen) >= limit:
                return list(seen.values())
    return list(seen.values())


def forbidden_patterns(graph, k: int, W: Sequence[int],
                       solver: str = "cd19") -> List[Tuple]:
    """Which colour patterns on W no k-colouring of the graph realises.

    Forcing, as this package has used it, is the smallest case of a much more
    general question. "Some target carries the pivot's colour" says one
    particular partition of {pivot} u T -- the one where the pivot is alone --
    is unrealisable. Nothing restricts the question to that shape.

    So: for a small set W, ask of every set partition of W whether some
    k-colouring of the graph induces it. The ones that do not are forced
    constraints, and a set with several of them is carrying far more
    information than a single disjunction does.

    On a graph with no k-colouring at all every partition comes back
    forbidden, which is vacuous rather than informative -- the same trap as
    measuring forcing on a 5-chromatic graph at k=4. Check colourability first.

    Returned as canonical partitions: a tuple of blocks, each a sorted tuple,
    ordered by their least element.
    """
    from itertools import product

    from pysat.formula import CNF
    from pysat.solvers import Solver

    W = list(W)
    n, w = graph.n, len(W)

    def x(v, c):
        return 1 + v * k + c

    base = CNF()
    for v in range(n):
        base.append([x(v, c) for c in range(k)])
    for u, v in graph.edges():
        for c in range(k):
            base.append([-x(u, c), -x(v, c)])

    out = []
    seen = set()
    for assign in product(range(k), repeat=w):
        blocks = {}
        for i, c in enumerate(assign):
            blocks.setdefault(c, []).append(i)
        part = tuple(sorted((tuple(b) for b in blocks.values()),
                            key=lambda t: t[0]))
        if part in seen:
            continue
        seen.add(part)
        cnf = CNF(from_clauses=base.clauses)
        for i, c in enumerate(assign):
            cnf.append([x(W[i], c)])
        with Solver(name=solver, bootstrap_with=cnf) as s:
            if not s.solve():
                out.append(part)
    return out


def pattern_pressure(graph, k: int, W: Sequence[int], **kw) -> float:
    """Fraction of W's partitions the graph forbids.

    Zero means W is unconstrained and carries nothing; one would mean no
    colouring survives at all. A single forced disjunction shows up here as
    one forbidden partition among many, which is how little of the available
    information the spindle argument uses.
    """
    forb = forbidden_patterns(graph, k, W, **kw)
    total = _bell(len(W))
    return len(forb) / total if total else 0.0


def _bell(n: int) -> int:
    row = [1]
    for _ in range(n):
        new = [row[-1]]
        for x in row:
            new.append(new[-1] + x)
        row = new
    return row[0]


def blocks_two_targets(graph, pivot: int, targets: Sequence[int],
                       isometries: Sequence) -> bool:
    """Exact blocking test for a core of two, in linear time, by 2-SAT.

    With exactly two targets each copy makes a binary choice, and a conflict
    between two chosen images is a forbidden pair -- which is a 2-SAT clause.
    An escape is a satisfying assignment, so the copies block exactly when the
    instance is unsatisfiable.

    That matters for reach, not elegance. `cross_blocks` searches the choices
    directly, which is 2^m in the number of copies and confines the test to a
    handful of them; 2-SAT decides it in linear time, so a core of two can be
    thrown against hundreds of copies at once. And since blocking is monotone
    upward -- an escape for a larger set restricts to one for any subset --
    more copies can only help here, which is not true for larger cores where
    each copy also multiplies the choices.

    Implemented by implication-graph strongly connected components: the
    instance is unsatisfiable exactly when some variable shares a component
    with its negation.
    """
    from .multispindle import cross_conflict_graph

    targets = list(targets)
    if len(targets) != 2:
        raise ValueError("this test is for exactly two targets")
    m = len(isometries)
    adj = cross_conflict_graph(graph, pivot, targets, isometries)

    # variable i true  <=> copy i chooses targets[0]
    def lit(i, first):
        return 2 * i + (0 if first else 1)

    def neg(l):
        return l ^ 1

    imp = [[] for _ in range(2 * m)]
    for i in range(m):
        for a, qa in enumerate(targets):
            for j in range(m):
                if j == i:
                    continue
                for b, qb in enumerate(targets):
                    if (j, qb) in adj[(i, qa)]:
                        # not (i chooses qa and j chooses qb)
                        imp[lit(i, a == 0)].append(neg(lit(j, b == 0)))
                        imp[lit(j, b == 0)].append(neg(lit(i, a == 0)))

    n = 2 * m
    index, low, on, stack, comp = [0] * n, [0] * n, [False] * n, [], [-1] * n
    counter = [1, 0]

    def strong(v0):
        work = [(v0, 0)]
        while work:
            v, pi = work[-1]
            if pi == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on[v] = True
            recurse = False
            for i in range(pi, len(imp[v])):
                w = imp[v][i]
                if index[w] == 0:
                    work[-1] = (v, i + 1)
                    work.append((w, 0))
                    recurse = True
                    break
                if on[w]:
                    low[v] = min(low[v], index[w])
            if recurse:
                continue
            if low[v] == index[v]:
                while True:
                    w = stack.pop()
                    on[w] = False
                    comp[w] = counter[1]
                    if w == v:
                        break
                counter[1] += 1
            work.pop()
            if work:
                u = work[-1][0]
                low[u] = min(low[u], low[v])

    for v in range(n):
        if index[v] == 0:
            strong(v)
    return any(comp[2 * i] == comp[2 * i + 1] for i in range(m))


def two_orbit_block(graph, pivot: int, a: int, b: int):
    """Six copies that close any forced pair with one leg on d^2 = 1/3.

    A theorem rather than a search. Let p be the pivot, |p - a|^2 = 1/3, and
    let d be the distance from p to b. Take rho, the 120 degree rotation about
    p, and sigma, the rotation by the unit-chord angle of b's circle, so that
    two points of that circle sigma apart are exactly one apart. The six
    copies are the two orbits

        rho^0, rho^1, rho^2   and   sigma rho^0, sigma rho^1, sigma rho^2.

    They block, for three reasons that fit together:

    * On the circle of radius 1/sqrt(3) two points are adjacent exactly when
      they are 120 degrees apart, so an orbit's three images of `a` form a
      unit triangle. Every chosen image carries the pivot's colour, so **at
      most one copy per orbit may choose a**, and at least two must choose b.

    * Two subsets of size at least two of a three-element set intersect. So
      there is a j for which rho^j and sigma rho^j both choose b.

    * Those two b-images differ by sigma, which is the unit-chord angle of
      their own circle, so they are one apart -- and both carry the pivot's
      colour. Contradiction.

    Returns the six copies, or None when sigma cannot be named in the field.
    The caller still has to supply the forcing; this only closes it.
    """
    from fractions import Fraction

    from .geometry import Rotation, rotation_joining

    p = graph.vertices[pivot]
    field = p.x.field
    da2 = p.dist2(graph.vertices[a])
    if not da2.is_rational() or da2.c[0] != Fraction(1, 3):
        raise ValueError("a must sit on the classical circle, d^2 = 1/3")
    db2 = p.dist2(graph.vertices[b])
    if not db2.is_rational():
        return None
    try:
        rho = rotation_joining(Fraction(1, 3), field)     # 120 degrees
        sigma = rotation_joining(db2.c[0], field)
    except (ValueError, ZeroDivisionError):
        return None

    def compose(x, y):
        return Rotation(x.cos * y.cos - x.sin * y.sin,
                        x.cos * y.sin + x.sin * y.cos, check=False)

    ident = Rotation(field.rational(1), field.zero())
    orbit = [ident, rho, compose(rho, rho)]
    return orbit + [compose(sigma, r) for r in orbit]


def unit_triangle_centroids(graph, limit: int = 20000) -> List[Point]:
    """The centroid of every unit triangle, which is 1/sqrt(3) from its corners.

    The two-orbit block needs a target at squared distance 1/3 from the pivot,
    and points that close together are scarce in the constructions here -- most
    of them are built from unit steps and land no nearer than that. But the
    distance is manufacturable. Three points pairwise one apart have a centroid
    exactly 1/sqrt(3) from each of them, so every unit triangle in a graph
    donates a point with three legs on the classical circle at once.

    The centroid is (u + v + w)/3, exact in the field, and the distance is
    checked rather than assumed.
    """
    from fractions import Fraction

    third = Fraction(1, 3)
    out, seen = [], set()
    n = graph.n
    for u in range(n):
        nbrs = sorted(graph.adj[u])
        for i, v in enumerate(nbrs):
            if v < u:
                continue
            for w in nbrs[i + 1:]:
                if w < u or w not in graph.adj[v]:
                    continue
                a, b, c = graph.vertices[u], graph.vertices[v], graph.vertices[w]
                g = Point((a.x + b.x + c.x) * a.x.field.rational(third),
                          (a.y + b.y + c.y) * a.y.field.rational(third))
                if g in seen:
                    continue
                d = g.dist2(a)
                if not (d.is_rational() and d.c[0] == third):
                    continue
                seen.add(g)
                out.append(g)
                if len(out) >= limit:
                    return out
    return out


def odd_orbit_block(graph, pivot: int, a: int, b: int, order: int = 3,
                    step: int = 1):
    """Two orbits of an odd-order rotation close a forced pair. 2n copies.

    The six-copy version is the case n = 3. Nothing in its argument needed the
    triangle specifically -- only that the cycle be odd.

    Let rho be the rotation by 2*pi/order about the pivot, and suppose a lies
    on that order's magic circle, radius 1/(2 sin(pi*step/order)), so two of
    its images `step` apart are exactly one apart. An orbit's images of a then
    form the circulant C_order(step); for step coprime to order that is an
    odd cycle, whose independent sets have at most (order-1)/2 elements. Every
    chosen image carries the pivot's colour, so **at most (order-1)/2 copies
    per orbit may choose a**, and at least (order+1)/2 must choose b.

    Take a second orbit offset by sigma, the rotation by the unit-chord angle
    of b's own circle. Two subsets of {0..order-1} each of size at least
    (order+1)/2 have total size at least order+1, so they intersect. At the
    shared index j both rho^j and sigma rho^j choose b, their images differ by
    sigma, and so they are one apart while both carry the pivot's colour.

    Returns the 2*order copies, or None when a rotation cannot be named in the
    field. `order` must be odd and coprime to `step`.
    """
    from math import gcd

    from .geometry import Rotation, rotation_joining

    if order < 3 or order % 2 == 0 or gcd(order, step) != 1:
        raise ValueError("order must be odd, at least 3, and coprime to step")
    p = graph.vertices[pivot]
    field = p.x.field
    db2 = p.dist2(graph.vertices[b])
    if not db2.is_rational():
        return None

    # rho: the rotation by 2*pi/order. For order 3 that is the joining
    # rotation of d^2 = 1/3; in general it is built from the same identity,
    # cos(2 pi/order) = 1 - 1/(2 d^2) at the magic radius d.
    if order == 3:
        try:
            rho = rotation_joining(Fraction(1, 3), field)
        except (ValueError, ZeroDivisionError):
            return None
    else:
        rho = _rotation_of_order_in(field, order)
        if rho is None:
            return None
    try:
        sigma = rotation_joining(db2.c[0], field)
    except (ValueError, ZeroDivisionError):
        return None

    def compose(x, y):
        return Rotation(x.cos * y.cos - x.sin * y.sin,
                        x.cos * y.sin + x.sin * y.cos, check=False)

    ident = Rotation(field.rational(1), field.zero())
    orbit, cur = [ident], ident
    for _ in range(order - 1):
        cur = compose(cur, rho)
        orbit.append(cur)
    return orbit + [compose(sigma, r) for r in orbit]


def _rotation_of_order_in(field, order: int):
    """A rotation of exactly this order over the field, or None."""
    from .multispindle import _rotation_of_order

    try:
        return _rotation_of_order(order, field)
    except (ValueError, ZeroDivisionError):
        return None


def circle_intersections(c1: Point, r1sq, c2: Point, r2sq) -> List[Point]:
    """The points at squared distance r1sq from c1 and r2sq from c2.

    `unit_circle_intersections` is the case r1sq = r2sq = 1. The general form
    is what a *constructed* configuration needs: placing a pivot at a chosen
    distance from two chosen points, rather than taking whatever the graph
    already offers.

    With d = c2 - c1 and D = |d|^2, the foot of the perpendicular sits at
    t = (D + r1sq - r2sq) / (2 D) along d, and the offset is h with
    h^2 = r1sq - t^2 D, carried on d's perpendicular scaled by 1/sqrt(D).
    Both square roots go through the same denesting as everything else, and an
    empty list means the field cannot name the point, not that none exists.
    """
    d = c2 - c1
    D = d.x * d.x + d.y * d.y
    if D == 0:
        return []
    field = D.field
    two = field.rational(2)
    t = (D + r1sq - r2sq) / (two * D)
    h2 = r1sq - t * t * D
    scale = _sqrt_in_field(h2 / D)
    if scale is None:
        return []
    fx, fy = c1.x + t * d.x, c1.y + t * d.y
    out = []
    for sign in (1, -1):
        s = field.rational(sign) * scale
        q = Point(fx - s * d.y, fy + s * d.x)
        if q.dist2(c1) == r1sq and q.dist2(c2) == r2sq:
            out.append(q)
    return out


def joining_rotation_exists(field, d2) -> bool:
    """Whether the field can name the rotation closing a circle of this radius.

    Two points of a circle of squared radius d2 are one apart at the angle
    with cos = 1 - 1/(2 d2); the rotation needs its sine too. The radius must
    also be at least 1/4, or the circle has no unit chord at all.
    """
    one = field.rational(1)
    two = field.rational(2)
    if isinstance(d2, (int, Fraction)):
        d2 = field.rational(d2)
    if d2 == 0:
        return False
    c = one - one / (two * d2)
    return _sqrt_in_field(one - c * c) is not None


# -- a core of two that is forced only jointly ----------------------------
#
# Every forced pair found by searching this package's graphs turned out to be
# forced one leg at a time: the 409-vertex certificate blocks a pair whose 1/3
# leg is already forced on its own, a core of one, so two copies of the
# ordinary spindle would close it and the block was unnecessary. The case the
# classical argument cannot reach is a pair where neither leg is forced by
# itself -- there is then no single target to spindle, and the colours are
# exhausted only by the pair together. This is such a configuration, built
# rather than found.

def joint_core_configuration():
    """Four points and two rotations realising a genuinely joint core of 2.

    Returns (field, [p, x, y, z], rho, sigma).

        x = (0,0)   y = (1,0)   z = (1/2, sqrt(3)/2)      unit triangle
        p = (5/6, -sqrt(11)/6)

        |p-x|^2 = 1     |p-y|^2 = 1/3     |p-z|^2 = (7 + sqrt(33))/6

    In any 3-colouring x, y, z take three colours and p differs from x, so p
    repeats y's or z's -- and both completions exist, so neither leg alone.
    The core is mixed, which is what carries it past the capacity ceiling:
    same-distance cores are capped at 2 and same-distance blocks at 3 colours.

    rho is the 120-degree rotation, the one odd order Niven's theorem leaves on
    a rational radius, closing the 1/3 circle. sigma closes the other:

        cos t = (-5 + 3 sqrt 33)/16,   sin^2 t = (-66 + 30 sqrt 33)/256

    That sin^2 has a negative conjugate, so its square root lies in no totally
    real field, and every multiquadratic field is totally real -- hence the one
    real quadratic extension. The rotation is perfectly ordinary in the plane,
    about 40.1 degrees; it is the arithmetic that needs room.
    """
    from .realext import RealQuadExt

    base = Field((3, 11))
    v = base.rational(Fraction(-66, 256)) \
        + base.sqrt(33) * base.rational(Fraction(30, 256))
    field = RealQuadExt(base, v)

    def pt(ax, bx, rx, ay, by, ry):
        cx = base.rational(ax) + (base.sqrt(rx) if rx else base.rational(1)) \
            * base.rational(bx)
        cy = base.rational(ay) + (base.sqrt(ry) if ry else base.rational(1)) \
            * base.rational(by)
        return Point(field.embed(cx), field.embed(cy))

    pts = [pt(Fraction(5, 6), 0, 0, 0, Fraction(-1, 6), 11),   # p
           pt(0, 0, 0, 0, 0, 0),                               # x
           pt(1, 0, 0, 0, 0, 0),                               # y
           pt(Fraction(1, 2), 0, 0, 0, Fraction(1, 2), 3)]     # z
    rho = Rotation(field.rational(Fraction(-1, 2)),
                   field.embed(base.sqrt(3) * base.rational(Fraction(1, 2))))
    sigma = Rotation(
        field.embed(base.rational(Fraction(-5, 16))
                    + base.sqrt(33) * base.rational(Fraction(3, 16))),
        field.radical())
    return field, pts, rho, sigma


def compose_rotations(u, w) -> Rotation:
    """u after w, without re-checking cos^2 + sin^2 on every product."""
    return Rotation(u.cos * w.cos - u.sin * w.sin,
                    u.cos * w.sin + u.sin * w.cos, check=False)


def joint_core_copies(field, rho, sigma):
    """The six isometries of the block: the 120-degree orbit, twice."""
    ident = Rotation(field.rational(1), field.zero())
    orbit = [ident, rho, compose_rotations(rho, rho)]
    return orbit + [compose_rotations(sigma, r) for r in orbit]


def joint_core_union():
    """The union of the six copies: 19 points with no 3-colouring."""
    field, pts, rho, sigma = joint_core_configuration()
    p = pts[0]
    copies = joint_core_copies(field, rho, sigma)
    return list(dict.fromkeys(r.about(p)(q) for r in copies for q in pts))
