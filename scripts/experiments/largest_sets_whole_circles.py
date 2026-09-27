"""Do the largest independent sets of G_11 and G_13 contain a point together with a whole circle?

G_q = Cay(F_q^2, {x^2 - n y^2 = 1}), n the least non-residue mod q, vertex v = q x + y. A circle
C_c = {z : N(z) = c} is independent when no two of its points are adjacent; a set S contains a point together
with a whole circle when p + C_c lies in S for some p in S and some independent circle C_c. Written for the
correction of 27 September in the research log (section "`α(G₁₇)`: the best sets are rosettes").

  largest_sets_whole_circles.py 11 28   lists the orbits of independent sets with at least 28 points of G_11 under
      its automorphisms z -> l z + s, l conj(z) + s (N(l) = 1), by a SAT solver that blocks every image of each
      set it finds until none is left (a complete enumeration, without a proof certificate; about 8 minutes),
      checks that each set found is maximal (so no independent set has more than t points: one would contain
      an image of a set found), and prints the whole circles in each;
  largest_sets_whole_circles.py 13 --known   checks the 15 independent sets of 36 points of G_13 listed below,
      found by an earlier search: each is independent and maximal, no two are in the same orbit, and it prints
      the whole circles in each. That these are all the orbits of 36-point sets is not claimed.
"""
import sys
import time
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

KNOWN_G13 = [
    [0, 1, 2, 3, 4, 5, 7, 11, 21, 23, 27, 30, 33, 35, 37, 47, 49, 59, 61, 63, 66, 69, 105, 108, 124, 126, 128, 138, 140, 144, 147, 150, 152, 154, 164, 166],
    [0, 1, 2, 3, 5, 9, 19, 25, 28, 29, 45, 49, 57, 59, 61, 63, 67, 68, 77, 78, 83, 84, 87, 88, 92, 93, 95, 103, 107, 122, 126, 149, 150, 152, 155, 164],
    [0, 1, 2, 3, 5, 11, 17, 21, 25, 27, 28, 45, 49, 57, 59, 60, 61, 63, 66, 67, 83, 84, 88, 89, 92, 93, 95, 103, 122, 128, 138, 147, 150, 152, 155, 164],
    [0, 1, 2, 3, 7, 9, 17, 21, 25, 27, 28, 57, 59, 61, 63, 83, 84, 86, 88, 89, 92, 93, 95, 103, 104, 107, 122, 128, 138, 147, 150, 152, 155, 161, 164, 167],
    [0, 1, 2, 4, 5, 11, 21, 23, 27, 30, 33, 35, 37, 44, 47, 49, 59, 61, 63, 66, 69, 105, 108, 124, 126, 128, 135, 138, 140, 144, 147, 150, 152, 154, 164, 166],
    [0, 1, 2, 4, 6, 8, 18, 24, 27, 28, 29, 34, 36, 38, 44, 48, 56, 60, 62, 63, 67, 83, 86, 90, 92, 106, 121, 125, 128, 131, 148, 149, 151, 155, 163, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 62, 63, 65, 67, 79, 121, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 62, 63, 65, 67, 92, 121, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 62, 63, 65, 92, 104, 121, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 62, 63, 79, 104, 106, 121, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 62, 63, 92, 104, 106, 121, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 63, 65, 67, 92, 121, 122, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 57, 59, 60, 63, 65, 92, 104, 121, 122, 124, 125, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 59, 60, 63, 65, 67, 92, 121, 122, 124, 125, 127, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167],
    [0, 1, 2, 17, 24, 29, 31, 33, 34, 36, 38, 40, 45, 48, 56, 59, 60, 63, 79, 104, 106, 121, 122, 124, 125, 127, 128, 131, 136, 139, 146, 148, 153, 155, 160, 167]
]


def plane(q):
    n = next(a for a in range(2, q) if pow(a, (q - 1) // 2, q) == q - 1)
    P = [(x, y) for x in range(q) for y in range(q)]
    I = {z: q * z[0] + z[1] for z in P}
    N = lambda z: (z[0] * z[0] - n * z[1] * z[1]) % q
    add = lambda a, b: ((a[0] + b[0]) % q, (a[1] + b[1]) % q)
    mul = lambda a, b: ((a[0] * b[0] + n * a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)
    U = [z for z in P if N(z) == 1]
    adj = {z: {add(z, u) for u in U} for z in P}
    C = {c: [z for z in P if N(z) == c] for c in range(1, q)}
    indep = [c for c in C if all(b not in adj[a] for a in C[c] for b in C[c])]
    maps = [[I[add(mul(l, (z[0], (-z[1]) % q) if cj else z), s)] for z in P] for l in U for cj in (0, 1) for s in P]
    return P, I, adj, C, indep, maps, add


def circles(S, P, C, indep, add):
    Z = {P[v] for v in S}
    found = {tuple(c for c in indep if all(add(p, w) in Z for w in C[c])) for p in Z} - {()}
    return sorted(found) or "none"


def enumerate_orbits(q, t):
    P, I, adj, C, indep, maps, add = plane(q)
    edges = {tuple(sorted((I[z], I[w]))) for z in P for w in adj[z]}
    print(f"G_{q}: {len(P)} vertices, {len(edges)} edges, {len(maps)} automorphisms; independent circles {indep}")
    sol = Solver(name="cd19")
    for a, b in edges:
        sol.add_clause([-(a + 1), -(b + 1)])
    for cl in CardEnc.atleast([v + 1 for v in range(len(P))], bound=t, top_id=len(P), encoding=EncType.seqcounter).clauses:
        sol.add_clause(cl)
    sol.add_clause([1])                      # by translation, the set contains vertex 0
    t0, orbits = time.time(), []
    while sol.solve():
        m = sol.get_model()
        S = [v for v in range(len(P)) if m[v] > 0]
        images = {frozenset(p[v] for v in S) for p in maps}
        for J in images:
            if 0 in J:
                sol.add_clause([-(v + 1) for v in J])
        orbits.append((S, len(maps) // len(images)))
    print(f"complete: {len(orbits)} orbits of independent sets with at least {t} points ({time.time() - t0:.0f} s)")
    for S, stab in orbits:
        Z = {P[v] for v in S}
        assert all(z in Z or adj[z] & Z for z in P), "a set found is not maximal"
        print(f"  {len(S)} points, maximal, stabiliser of order {stab}: whole circles around a point: "
              f"{circles(S, P, C, indep, add)}")


def check_known():
    P, I, adj, C, indep, maps, add = plane(13)
    canon = set()
    for S in KNOWN_G13:
        Z = {P[v] for v in S}
        assert len(Z) == 36 and all(not (adj[a] & Z) for a in Z), "not an independent set of 36 points"
        assert all(z in Z or adj[z] & Z for z in P), "not maximal"
        canon.add(min(tuple(sorted(m[v] for v in S)) for m in maps))
        print(f"  36 points, maximal: whole circles around a point: {circles(S, P, C, indep, add)}")
    print(f"G_13: {len(KNOWN_G13)} sets in {len(canon)} distinct orbits; independent circles {indep}")


if __name__ == "__main__":
    if sys.argv[1:] == ["13", "--known"]:
        check_known()
    else:
        enumerate_orbits(int(sys.argv[1]), int(sys.argv[2]))
