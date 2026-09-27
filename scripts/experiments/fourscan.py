"""Four hexagons, three free angles: the last lever this route has.

The peak method extends one hexagon at a time.  Given offsets already fixed,
an edge joining an auxiliary that uses the NEW hexagon to one that does not
reads |r_new w^m - X| = 1, where X is the old auxiliary minus the old circle
point the new one is paired with.  That is the same line-meets-circle problem
as before, Re(r_new w^m conj(X)) = |X|^2 / 2, solvable whenever |X| <= 2, so
the rich fourth offsets are again the histogram peaks.

Three hexagons gave 9600 orientations with degeneracy never above 2.  A fourth
adds a free angle and roughly doubles the auxiliaries, which is the last
straightforward way to find the minimum degree 3 the list argument needs.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, cmath, math, itertools, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

W = cmath.exp(1j * math.pi / 3)
HEX = [W ** k for k in range(6)]
SIXTY = math.pi / 3


def auxiliaries(angles):
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
    return hexes, aux


def evaluate(angles):
    r = auxiliaries(angles)
    if r is None:
        return None
    hexes, aux = r
    keys = list(aux)
    pts = [complex(*k) for k in keys]
    adj = {i: set() for i in range(len(pts))}
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(abs(pts[i] - pts[j]) - 1) < 1e-7:
                adj[i].add(j)
                adj[j].add(i)
    best = 0
    for bits in itertools.product((0, 1), repeat=len(hexes)):
        col = {(a, i): (i + bits[a]) % 2 for a in range(len(hexes))
               for i in range(6)}
        conf = {i for i in range(len(pts))
                if {col[x] for x in aux[keys[i]]} == {0, 1}}
        rem = {v: len(adj[v] & conf) for v in conf}
        d = 0
        while rem:
            v = min(rem, key=rem.get)
            d = max(d, rem[v])
            for u in adj[v]:
                if u in rem:
                    rem[u] -= 1
            del rem[v]
        best = max(best, d)
    return best, len(pts)


def next_offsets(angles):
    """Histogram peaks for one more hexagon, given the ones already placed."""
    r = auxiliaries(angles)
    if r is None:
        return []
    hexes, aux = r
    olds = [complex(*k) for k in aux]
    out = collections.Counter()
    for m in range(6):
        wm = HEX[m]
        for a, h in enumerate(hexes):
            for i, u in enumerate(h):
                for q in olds:
                    X = q - u
                    n = abs(X)
                    if n > 2 or n < 1e-12:
                        continue
                    C = wm * X.conjugate()
                    base, off = cmath.phase(C), math.acos(min(1.0, n / 2))
                    for s in (1, -1):
                        out[round((-base + s * off) % SIXTY, 6)] += 1
    bad = [a % SIXTY for a in angles]
    return [(a, m) for a, m in out.most_common()
            if all(min((a - b) % SIXTY, (b - a) % SIXTY) > 1e-4 for b in bad)]


theta = math.acos(5 / 6)
print("validation -- adding a fourth hexagon to de Grey's three:")
top = next_offsets([0.0, theta / 2, theta])[:3]
for a, m in top:
    print(f"  peak {math.degrees(a):8.4f} deg, multiplicity {m} -> "
          f"{evaluate([0.0, theta / 2, theta, a])}", flush=True)

print("\nscanning three free angles by successive peaks:")
best = (0, None)
count = 0
N = 60
for n in range(N):
    a1 = SIXTY * (n + 0.5) / N
    for a2, _ in next_offsets([0.0, a1])[:2]:
        for a3, _ in next_offsets([0.0, a1, a2])[:3]:
            r = evaluate([0.0, a1, a2, a3])
            if r is None:
                continue
            count += 1
            if r[0] > best[0]:
                best = (r[0], (a1, a2, a3, r[1]))
                print(f"  {math.degrees(a1):7.3f} {math.degrees(a2):7.3f} "
                      f"{math.degrees(a3):7.3f} -> degeneracy {r[0]} "
                      f"on {r[1]} auxiliaries"
                      + ("   *** 3 OR MORE ***" if r[0] >= 3 else ""),
                      flush=True)
print(f"\n{count} four-hexagon configurations; best degeneracy {best[0]}")
