"""Is every confined degree even?  If so the target is harder than degeneracy 3.

On de Grey's configuration every confined auxiliary has degree 0, 2 or 4 --
never odd -- with 42 of 72 isolated.  If that is a property of the
construction rather than of one angle, then a subgraph of minimum degree 3 is
impossible outright and reaching degeneracy 3 needs minimum degree FOUR: a
4-regular subgraph, which is a far stronger demand than the list-colouring
bound by itself suggests.

The likely cause is a symmetry.  Reflecting across the line through the pivot
and an auxiliary maps the configuration to itself when the hexagon offsets are
symmetric about it, pairing each neighbour with another, so degrees come in
twos.  This checks the parity over the whole scanned space rather than
assuming the symmetry holds everywhere.
"""
import sys, cmath, math, itertools, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")

W = cmath.exp(1j * math.pi / 3)
HEX = [W ** k for k in range(6)]
SIXTY = math.pi / 3


def build(angles):
    for x, y in itertools.combinations(angles, 2):
        if min((x - y) % SIXTY, (y - x) % SIXTY) < 1e-4:
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
                    if abs(abs(s) - 1) < 1e-9 or abs(s) < 1e-9:
                        continue
                    k = (round(s.real, 9), round(s.imag, 9))
                    aux.setdefault(k, set()).add((a, i))
                    aux[k].add((b, j))
    keys = list(aux)
    pts = [complex(*k) for k in keys]
    adj = {i: set() for i in range(len(pts))}
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(abs(pts[i] - pts[j]) - 1) < 1e-7:
                adj[i].add(j)
                adj[j].add(i)
    return keys, aux, adj


def candidates(a1):
    out = collections.Counter()
    e1 = cmath.exp(1j * a1)
    for i in range(6):
        for ip in range(6):
            if i == ip:
                continue
            d = HEX[i] - HEX[ip]
            for j in range(6):
                A = d + e1 * HEX[j]
                m = abs(A)
                if m > 2 or m < 1e-12:
                    continue
                for l in range(6):
                    C = A * HEX[l].conjugate()
                    base, off = cmath.phase(C), math.acos(min(1.0, m / 2))
                    for s in (1, -1):
                        out[round((base + s * off) % SIXTY, 6)] += 1
    return out


odd_seen, checked, maxdeg, best_min = 0, 0, 0, 0
N = 300
for n in range(N):
    a1 = SIXTY * (n + 0.5) / N
    c = candidates(a1)
    usable = [(a, m) for a, m in c.most_common()
              if min(a % SIXTY, (-a) % SIXTY) > 1e-4
              and min((a - a1) % SIXTY, (a1 - a) % SIXTY) > 1e-4]
    for a2, _ in usable[:4]:
        b = build([0.0, a1, a2])
        if b is None:
            continue
        keys, aux, adj = b
        checked += 1
        for bits in itertools.product((0, 1), repeat=3):
            col = {(a, i): (i + bits[a]) % 2 for a in range(3) for i in range(6)}
            conf = {i for i in range(len(keys))
                    if {col[x] for x in aux[keys[i]]} == {0, 1}}
            degs = [len(adj[v] & conf) for v in conf]
            if any(d % 2 for d in degs):
                odd_seen += 1
            maxdeg = max([maxdeg] + degs)
            # the largest minimum degree over any subgraph, which is the
            # degeneracy -- reported to confirm nothing reaches 3
            rem = {v: len(adj[v] & conf) for v in conf}
            d = 0
            while rem:
                v = min(rem, key=rem.get)
                d = max(d, rem[v])
                for u in adj[v]:
                    if u in rem:
                        rem[u] -= 1
                del rem[v]
            best_min = max(best_min, d)

print(f"{checked} configurations x 8 orientations")
print(f"  orientations with any ODD confined degree: {odd_seen}")
print(f"  largest confined degree seen: {maxdeg}")
print(f"  largest degeneracy seen: {best_min}")
