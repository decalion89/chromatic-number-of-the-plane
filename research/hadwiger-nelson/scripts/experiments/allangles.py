"""Every relative angle that gives the confined set an edge at all, enumerated.

Two auxiliaries q = u + v and q' = u' + v' drawn from hexagons a and b differ
by D1 + D2 with D1 = r_a(w^i - w^i') and D2 = r_b(w^j - w^j'), w = exp(i pi/3).
Writing z = D1 conj(D2) / (r_a conj(r_b)) -- an Eisenstein integer -- and
t = r_a conj(r_b) for the relative rotation of the two hexagons,

    |q - q'|^2 = |D1|^2 + |D2|^2 + 2 Re(t z) = 1,

so Re(t z) = (1 - |D1|^2 - |D2|^2)/2, a LINE in t meeting the unit circle in at
most two points.  There are 31 choices of each D, so at most 1922 candidate
angles: the design space is finite and enumerable, not a continuum.

One case is free and useless.  When D2 = 0 -- same v, adjacent u -- the
condition is |D1| = 1 and holds at every angle; but confinement of u_i + v_j
depends on i + j, and adjacent i differ by one, so those two auxiliaries never
have the same parity.  Every always-available edge joins a confined auxiliary
to a non-confined one.  The confined set's own edges must come from the
discrete solutions, which is why Pythagorean rotations give it none.

This enumerates them, measures the confined degeneracy at each, and reports
any that beats 2 -- the only thing that would reopen the route.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, cmath, itertools, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

W = cmath.exp(1j * math.pi / 3)
HEX = [W ** k for k in range(6)]
EPS = 1e-9

# every difference within a hexagon, as an Eisenstein integer
DIFFS = sorted({(round((HEX[i] - HEX[j]).real, 9),
                 round((HEX[i] - HEX[j]).imag, 9))
                for i in range(6) for j in range(6)})
DIFFS = [complex(*d) for d in DIFFS]
print(f"{len(DIFFS)} hexagon differences (including zero)")

cands = set()
for D1 in DIFFS:
    for D2 in DIFFS:
        if abs(D2) < EPS:
            continue                      # the free case, never confined-confined
        z = D1 * D2.conjugate()
        if abs(z) < EPS:
            continue
        R = (1 - abs(D1) ** 2 - abs(D2) ** 2) / 2
        if abs(R) > abs(z) + EPS:
            continue                      # the line misses the circle
        base = cmath.phase(z)
        c = max(-1.0, min(1.0, R / abs(z)))
        for sign in (1, -1):
            cands.add(round((sign * math.acos(c) - base) % (2 * math.pi), 9))

# hexagons are invariant under 60 degrees, so angles reduce mod pi/3
red = sorted({round(a % (math.pi / 3), 7) for a in cands})
print(f"{len(cands)} raw solutions -> {len(red)} distinct angles mod 60 degrees")
theta = math.acos(5 / 6)
print(f"  Moser theta/2 = {math.degrees(theta/2):.4f} deg is "
      f"{'PRESENT' if any(abs(a - theta/2 % (math.pi/3)) < 1e-5 for a in red) else 'absent'}")


def degeneracy(adj, verts):
    rem = {v: len(adj[v] & verts) for v in verts}
    d = 0
    while rem:
        v = min(rem, key=rem.get)
        d = max(d, rem[v])
        for u in adj[v]:
            if u in rem:
                rem[u] -= 1
        del rem[v]
    return d


def measure(angles):
    """Confined degeneracy of hexagons at these offsets, worst orientation."""
    # Hexagons are invariant under 60 degrees, so two offsets that agree
    # modulo 60 are the SAME six points.  A tolerance of 1e-9 on the points
    # missed near-coincidences and reported a spurious degeneracy 4 at an
    # angle that was 60 degrees to within 1e-7 -- the coincident-hexagon trap
    # again, this time in floating point.  Compare the OFFSETS, generously.
    for x, y in itertools.combinations(angles, 2):
        if min((x - y) % (math.pi / 3), (y - x) % (math.pi / 3)) < 1e-5:
            return None
    hexes = [[cmath.exp(1j * a) * h for h in HEX] for a in angles]
    aux = {}
    for a in range(len(hexes)):
        for b in range(a, len(hexes)):
            for i, u in enumerate(hexes[a]):
                for j, v in enumerate(hexes[b]):
                    if (a, i) >= (b, j):
                        continue
                    s = u + v
                    if abs(abs(s) - 1) < EPS or abs(s) < EPS:
                        continue
                    key = (round(s.real, 7), round(s.imag, 7))
                    aux.setdefault(key, set()).add((a, i))
                    aux[key].add((b, j))
    pts = list(aux)
    adj = {i: set() for i in range(len(pts))}
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(abs(complex(*pts[i]) - complex(*pts[j])) - 1) < 1e-7:
                adj[i].add(j)
                adj[j].add(i)
    best = 0
    for bits in itertools.product((0, 1), repeat=len(hexes)):
        col = {(a, i): (i + bits[a]) % 2 for a in range(len(hexes))
               for i in range(6)}
        conf = {i for i in range(len(pts))
                if {col[n] for n in aux[pts[i]]} == {0, 1}}
        best = max(best, degeneracy(adj, conf))
    return best, len(pts)


print("\nTWO hexagons, one at each enumerable angle.  The enumeration is")
print("complete only for t = 2: with three hexagons an edge can join")
print("auxiliaries from DIFFERENT pairs, giving two free angles against one")
print("equation -- a curve, not a discrete set.  de Grey's theta/2 is exactly")
print("such a case, which is why it is absent above and why chi(confined)")
print("jumps from 1 to 3 when the third hexagon is added.\n")
tally, winners = {}, []
for a in red:
    r = measure([0.0, a])
    if r is None:
        continue
    d, n = r
    tally[d] = tally.get(d, 0) + 1
    if d >= 3:
        winners.append((a, d, n))
for d in sorted(tally):
    print(f"  degeneracy {d}: {tally[d]} angles")
print(f"  beating 2: {len(winners)}")
for a, d, n in winners[:10]:
    print(f"    {math.degrees(a):.5f} deg -> degeneracy {d} on {n} auxiliaries")
