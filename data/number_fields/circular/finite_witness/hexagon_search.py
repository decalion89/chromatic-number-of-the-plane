"""Induced unit-distance realisations of the subdivided hexagon M over Q(sqrt d).

usage: python3 hexagon_search.py d D [count]     (for instance 7 80 400; about a minute)

M is the hexagon v_0 ... v_5 with its three long diagonals v_j v_{j+3} subdivided by vertices m_0, m_1, m_2 (9
vertices, 12 edges); it is the triangle-free graph with 9 vertices and chi_c = 3 that has the fewest edges
(nine_vertices.py).  The search takes the hexagon's edge vectors u_1, ..., u_6 among the unit vectors with denominator
D (u_6 = -(u_1 + ... + u_5) must be one too), with v_0 = 0 and v_i = v_{i-1} + u_i, requires each long diagonal
v_{j+3} - v_j = u_{j+1} + u_{j+2} + u_{j+3} to be a sum of two unit vectors c + c' (then m_j = v_j + c is at
distance 1 from v_j and v_{j+3}), and keeps the sets of nine distinct points found in this way.  It prints the first
`count` of them whose only pairs at distance 1 are the twelve edges of M, as lists of [a, b, c, e] (the point
((a + b sqrt d)/D, (c + e sqrt d)/D)), and at the end the number of point sets for each set of further pairs at
distance 1, with an example when these only join midpoints: M plus one or two such edges are the other two
nine-vertex graphs with chi_c = 3 (three would close a triangle).  Exact integer arithmetic throughout.
"""
import itertools, json, math, sys


def unit_vectors(d, D):
    """All (a, b, c, e) with a^2 + d b^2 + c^2 + d e^2 = D^2 and a b + c e = 0 (unit vectors with denominator D)."""
    out = set()
    B = math.isqrt(D * D // d)
    for b in range(-B, B + 1):
        for e in range(-B, B + 1):
            r = D * D - d * (b * b + e * e)
            if r < 0:
                continue
            for a in range(-D, D + 1):
                c2 = r - a * a
                if c2 < 0:
                    continue
                c = math.isqrt(c2)
                if c * c == c2:
                    for cc in {c, -c}:
                        if a * b + cc * e == 0:
                            out.add((a, b, cc, e))
    return sorted(out)


def main():
    d, D = int(sys.argv[1]), int(sys.argv[2])
    count = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    U = unit_vectors(d, D)
    Us = set(U)
    add = lambda a, b: tuple(x + y for x, y in zip(a, b))
    neg = lambda a: tuple(-x for x in a)

    def unit(p, q):
        a, b, c, e = (q[k] - p[k] for k in range(4))
        return a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0

    pairs = {}
    for u in U:
        for w in U:
            pairs.setdefault(add(u, w), []).append(u)
    print(f'd = {d}, D = {D}: {len(U)} unit vectors', flush=True)
    edges = {(i, (i + 1) % 6) for i in range(6)} | {(j, 6 + j) for j in range(3)} | {(j + 3, 6 + j) for j in range(3)}
    edges = {tuple(sorted(e)) for e in edges}
    names = ['v0', 'v1', 'v2', 'v3', 'v4', 'v5', 'm0', 'm1', 'm2']
    seen, found, patterns, example = set(), 0, {}, {}
    for u1, u2, u3 in itertools.product(U, repeat=3):
        s0 = add(add(u1, u2), u3)
        if s0 not in pairs or s0 == (0, 0, 0, 0):
            continue
        for u4 in U:
            s1 = add(add(u2, u3), u4)
            if s1 not in pairs:
                continue
            for u5 in U:
                s2 = add(add(u3, u4), u5)
                if s2 not in pairs or neg(add(add(s0, u4), u5)) not in Us:
                    continue
                v = [(0, 0, 0, 0)]
                for u in (u1, u2, u3, u4, u5):
                    v.append(add(v[-1], u))
                if len(set(v)) < 6:
                    continue
                for c0 in pairs[s0]:
                    for c1 in pairs[s1]:
                        for c2 in pairs[s2]:
                            pts = v + [add(v[0], c0), add(v[1], c1), add(v[2], c2)]
                            if len(set(pts)) < 9 or tuple(sorted(pts)) in seen:
                                continue
                            seen.add(tuple(sorted(pts)))
                            extra = tuple(names[i] + names[j] for i in range(9) for j in range(i + 1, 9)
                                          if (i, j) not in edges and unit(pts[i], pts[j]))
                            patterns[extra] = patterns.get(extra, 0) + 1
                            example.setdefault(extra, pts)
                            if not extra and found < count:
                                found += 1
                                print('M:', json.dumps([list(p) for p in pts]), flush=True)
    print(f'{len(seen)} sets of nine points realise the twelve edges of M; by further pairs at distance 1:')
    for extra, n in sorted(patterns.items(), key=lambda x: (len(x[0]), x[0])):
        print(f'  {n:6d}  {" ".join(extra) if extra else "none (M itself)"}')
        if extra and all(s[0] == 'm' and s[2] == 'm' for s in extra):
            print('          for example', json.dumps([list(p) for p in example[extra]]))


if __name__ == '__main__':
    main()
