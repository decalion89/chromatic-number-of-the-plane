"""The field of a 5-chromatic unit-distance graph, as a parameter.

Composing a forced pair with a rotated copy of its own carrier moves the pair
to any distance in [0, 2D]:

    tau(p) = q + R_phi (p - v),   |v - tau(q)|^2 = 2 D^2 (1 + cos phi)

so a target d^2 needs cos phi = d^2/(2 D^2) - 1, which is RATIONAL whenever the
two squared distances are.  Only sin phi carries a radical.  With D^2 = 64/9
and d^2 = p/q,

    cos phi = (9p - 128q) / (128q),     sin phi = 3 sqrt(p (256q - 9p)) / (128q)

and sqrt(p (256q - 9p)) is the entire arithmetic cost of the move.  Choosing p
and q therefore CHOOSES the field, which is the opposite of how every other
construction here works -- in those the radii pick the radicals and one takes
what one is given.

Two things fall out.

  * d^2 = 1 makes the composite pair ADJACENT, so the chained union is already
    5-chromatic and no spindle is needed at all: half the vertices.  Its
    radical is sqrt(247), the same one the ordinary spindle at 64/9 needs --
    a sanity check, since tuning to 1 is that spindle wearing a different hat.

  * every other target is genuinely new.  d^2 = 2 needs sqrt(119) = sqrt7
    sqrt17, d^2 = 3 needs sqrt(687) = sqrt3 sqrt229, d^2 = 1/3 needs
    sqrt(759) = sqrt3 sqrt11 sqrt23 -- while d^2 = 4 and d^2 = 2/3 land back
    inside de Grey's own field by a route he did not take.

Each case is built, its composite distance checked exactly, its forcing
confirmed by the solver rather than assumed, and the result tested for four
colours.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
from math import isqrt
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()

def squarefree(m):
    out, r = m, 1
    d = 2
    while d * d <= out:
        while out % (d * d) == 0:
            out //= d * d; r *= d
        d += 1
    return out, r

# target d^2 -> the generators the radical needs, on top of Sa's own
TARGETS = [Fr(1), Fr(1, 3), Fr(2), Fr(4)]
D2 = Fr(64, 9)
for d2 in TARGETS:
    p, q = d2.numerator, d2.denominator
    rad = p * (256 * q - 9 * p)
    sf, mult = squarefree(rad)
    # Sa already lives in Q(sqrt3, sqrt5, sqrt7, sqrt11), so only the part of
    # the radical those four do not already supply costs a new generator.  sf
    # is squarefree, so each of 3, 5, 7, 11 divides it at most once and can be
    # divided straight out; whatever survives becomes one generator of its own
    # (it need not be prime -- Field((3, 11, 247)) is already in the project).
    extra = sf
    for base in (3, 5, 7, 11):
        if extra % base == 0:
            extra //= base
    gens = tuple(sorted({3, 5, 7, 11} | ({extra} if extra > 1 else set())))
    print(f"\n  target d^2 = {d2}: radical sqrt({rad}) = {mult} sqrt({sf}), "
          f"field {gens}   [{time.time()-t0:.0f}s]", flush=True)
    if len(gens) > 6:
        print("    (too many generators, skipped)", flush=True); continue
    F = Field(gens)
    rot60 = _rot60(F)
    Sa = build_Sa(F)
    def orb(pt):
        out, z = [], pt
        for _ in range(6):
            if z not in out: out.append(z)
            z = rot60(z)
        for r in [Point(w.x, -w.y) for w in list(out)]:
            if r not in out: out.append(r)
        return out
    seen, U = set(Sa), list(Sa)
    for w in orb(Sa[25]):
        rot = rotation_joining(Fr(1), F).about(w)
        for pt in Sa:
            z = rot(pt)
            if z not in seen: seen.add(z); U.append(z)
    g = build_graph(U); n = g.n
    K = 4
    X = lambda u, c: 1 + u * K + c
    cnf = [[X(u, c) for c in range(K)] for u in range(n)]
    for u in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cnf.append([-X(u, a), -X(u, b)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="m22", bootstrap_with=cnf); assert s.solve()
    rng = random.Random(3)
    def read():
        pz = set(l for l in s.get_model() if l > 0)
        return [next(c for c in range(K) if X(u, c) in pz) for u in range(n)]
    cols = [read()]
    while len(cols) < 10:
        s.add_clause([-X(u, cols[-1][u]) for u in rng.sample(range(n), 30)])
        s.set_phases([(1 if rng.random() < 0.25 else -1) * X(u, c)
                      for u in range(n) for c in range(K)])
        if not s.solve(): break
        cols.append(read())
    buck = defaultdict(list)
    for u in range(n):
        buck[tuple(c[u] for c in cols)].append(u)
    DD = F.rational(D2)
    pair = None
    for vs in buck.values():
        for i, x in enumerate(vs):
            for y in vs[i+1:]:
                if (g.vertices[x] - g.vertices[y]).norm2() == DD and \
                   not s.solve(assumptions=[X(x, 0), X(y, 1)]):
                    pair = (x, y); break
            if pair: break
        if pair: break
    s.delete()
    assert pair
    v, qq = pair
    cosphi = F.rational(d2 / (2 * D2) - 1)
    # sin^2 = 1 - cos^2 = [(128q)^2 - (9p - 128q)^2] / (128q)^2
    #                   = 9 p (256q - 9p) / (128q)^2
    # so the 3 is part of the formula, not of the radical.
    sinphi = F.sqrt(rad) * F.rational(Fr(3, 128 * q))
    assert cosphi * cosphi + sinphi * sinphi == F.rational(1), "not a rotation"
    R = Rotation(cosphi, sinphi)
    V, Q = g.vertices[v], g.vertices[qq]
    def tau(pt):
        dd = Point(pt.x - V.x, pt.y - V.y)
        r = R(dd)
        return Point(Q.x + r.x, Q.y + r.y)
    assert tau(V) == Q
    tq = tau(Q)
    assert (V - tq).norm2() == F.rational(d2), "composite distance missed"
    seen2, W = set(g.vertices), list(g.vertices)
    for pt in g.vertices:
        z = tau(pt)
        if z not in seen2: seen2.add(z); W.append(z)
    g2 = build_graph(W); pos = {pt: i for i, pt in enumerate(g2.vertices)}
    def cnf_for(gr, KK):
        N = gr.n
        Y = lambda u, c: 1 + u * KK + c
        cl = [[Y(u, c) for c in range(KK)] for u in range(N)]
        for u in range(N):
            for a in range(KK):
                for b in range(a + 1, KK):
                    cl.append([-Y(u, a), -Y(u, b)])
        for x, y in gr.edges():
            for c in range(KK):
                cl.append([-Y(x, c), -Y(y, c)])
        return cl, Y
    cl2, Y2 = cnf_for(g2, 4)
    s2 = Solver(name="cd19", bootstrap_with=cl2)
    if d2 == 1:
        ok = s2.solve(); s2.delete()
        m2 = sum(len(a) for a in g2.adj) // 2
        print(f"    chain alone n={g2.n} m={m2}: 4-colourable {ok}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
        final = g2
    else:
        forced = not s2.solve(assumptions=[Y2(pos[V], 0), Y2(pos[tq], 1)])
        s2.delete()
        print(f"    composite forced at four: {forced}", flush=True)
        assert forced
        rho = rotation_joining(d2, F).about(V)
        assert (tq - rho(tq)).norm2() == F.rational(1)
        seen3, Z = set(g2.vertices), list(g2.vertices)
        for pt in g2.vertices:
            z = rho(pt)
            if z not in seen3: seen3.add(z); Z.append(z)
        g3 = build_graph(Z)
        cl3, _ = cnf_for(g3, 4)
        t1 = time.time()
        s3 = Solver(name="cd19", bootstrap_with=cl3); ok = s3.solve(); s3.delete()
        m3 = sum(len(a) for a in g3.adj) // 2
        print(f"    spindled n={g3.n} m={m3}: 4-colourable {ok}"
              f"   [{time.time()-t1:.0f}s]", flush=True)
        final = g3
    if not ok:
        mm = sum(len(a) for a in final.adj) // 2
        name = f"five_tuned_{p}_{q}.json"
        json.dump({"field_generators": list(F.gens), "n": final.n, "m": mm,
                   "target_d2": [p, q], "radical": rad,
                   "mechanism": "composed forcing tuned, then spindled"
                                if d2 != 1 else "composed forcing tuned to 1",
                   "points": [[[[c.numerator, c.denominator] for c in pt.x.c],
                               [[c.numerator, c.denominator] for c in pt.y.c]]
                              for pt in final.vertices]},
                  open(f"{ROOT}/data/{name}", "w"))
        print(f"    *** refuses four -- written data/{name} ***", flush=True)
