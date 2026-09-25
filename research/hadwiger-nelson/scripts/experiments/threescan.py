"""Three hexagons, two free angles: scan the space de Grey's family is a curve in.

With three hexagons an edge of the confined set can join auxiliaries from
different pairs -- q = u_i + v_j from (0,1) and q' = u_i' + w_l from (0,2) --
so for a fixed alpha1 the condition

    |(w^i - w^i') + e^(i alpha1) w^j - e^(i alpha2) w^l| = 1

solves for alpha2: writing A for the first two terms and C = A conj(w^l),
Re(C e^(-i alpha2)) = |A|^2 / 2, so alpha2 = arg C -+ arccos(|A|/2), real
whenever |A| <= 2.  Each (i, i', j, l) contributes up to two values, 2160 in
all, and a configuration is rich exactly where MANY of them coincide.

So the scan does not evaluate a grid -- it histograms the candidate alpha2 at
each alpha1 and evaluates only the peaks, which is where simultaneous edges
are.  de Grey's own configuration (alpha1 = theta/2, alpha2 = theta) has to
show up as one, and that is the check that the method sees anything at all.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, cmath, math, itertools, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

W = cmath.exp(1j * math.pi / 3)
HEX = [W ** k for k in range(6)]
SIXTY = math.pi / 3
TOL = 1e-7   # safe once keys are rounded to 1e-9, not 1e-6


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
    for x, y in itertools.combinations(angles, 2):
        if min((x - y) % SIXTY, (y - x) % SIXTY) < 1e-5:
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
                    key = (round(s.real, 9), round(s.imag, 9))
                    aux.setdefault(key, set()).add((a, i))
                    aux[key].add((b, j))
    pts = [complex(*k) for k in aux]
    keys = list(aux)
    adj = {i: set() for i in range(len(pts))}
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(abs(pts[i] - pts[j]) - 1) < TOL:
                adj[i].add(j)
                adj[j].add(i)
    best, edges = 0, sum(len(v) for v in adj.values()) // 2
    for bits in itertools.product((0, 1), repeat=len(hexes)):
        col = {(a, i): (i + bits[a]) % 2 for a in range(len(hexes))
               for i in range(6)}
        conf = {i for i in range(len(pts))
                if {col[n] for n in aux[keys[i]]} == {0, 1}}
        best = max(best, degeneracy(adj, conf))
    return best, len(pts), edges


theta = math.acos(5 / 6)
print("validation -- de Grey's own configuration:")
print(f"  alpha1 = theta/2 = {math.degrees(theta/2):.4f} deg, "
      f"alpha2 = theta = {math.degrees(theta):.4f} deg")
c = candidates(theta / 2)
peak = c[round(theta % SIXTY, 6)]
top = c.most_common(1)[0]
print(f"  its multiplicity {peak}, the largest at this alpha1 is {top[1]} "
      f"at {math.degrees(top[0]):.4f} deg")
print(f"  measured: {measure([0.0, theta / 2, theta])}\n")

print("scanning alpha1, evaluating the three richest alpha2 at each:")
best_overall = (0, None)
N = 240
for n in range(N):
    a1 = SIXTY * (n + 0.5) / N
    c = candidates(a1)
    # The largest peak is always alpha2 = alpha1 -- the same six points --
    # so coincident offsets are dropped before ranking, not after.
    usable = [(a, m) for a, m in c.most_common()
              if min(a % SIXTY, (-a) % SIXTY) > 1e-4
              and min((a - a1) % SIXTY, (a1 - a) % SIXTY) > 1e-4]
    for a2, mult in usable[:3]:
        r = measure([0.0, a1, a2])
        if r is None:
            continue
        d, npts, ne = r
        if d > best_overall[0]:
            best_overall = (d, (a1, a2, mult, npts, ne))
            print(f"  alpha1 {math.degrees(a1):8.4f}  alpha2 "
                  f"{math.degrees(a2):8.4f}  mult {mult:3}  "
                  f"{npts} aux, {ne} edges -> degeneracy {d}"
                  + ("   *** 3 OR MORE ***" if d >= 3 else ""), flush=True)
print(f"\nbest degeneracy over the scan: {best_overall[0]}")
