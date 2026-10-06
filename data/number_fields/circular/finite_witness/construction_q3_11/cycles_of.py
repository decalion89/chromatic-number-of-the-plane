"""cycles_of.py H.json OUT.json: replay the homotopies of per_build.py and list the cycles whose tightness Lemma P
needs to exclude: every square of a swap (both orientations) and one simple sub-cycle of each relation walk (both
orientations; only lengths divisible by 4 can be tight).  Writes {D, units, points, cycles} (cycles as point indices)
in the format of certify4.py."""
import json, sys
W = json.load(open(sys.argv[1])); U = W["units"]; dim = len(U[0])
pts = [tuple(p) for p in W["points"]]; idx = {p: i for i, p in enumerate(pts)}
R = [tuple(r) for r in W["relations"]]; T = W["triples"]


def steps_of(v, sign=1):
    st = []
    for j, c in enumerate(v):
        st += [(j, 1 if c > 0 else -1)] * abs(c)
    if sign < 0:
        st = [(j, -s) for (j, s) in reversed(st)]
    return st


def points_of(steps):
    p = [0] * dim; out = [tuple(p)]
    for j, s in steps:
        for t in range(dim):
            p[t] += s * U[j][t]
        out.append(tuple(p))
    return out


cyc = set(); nsq = 0
def put(c):
    """c: list of point indices of a simple cycle; store both orientations, canonically rotated"""
    for cc in (c, c[::-1]):
        k = cc.index(min(cc)); cyc.add(tuple(cc[k:] + cc[:k]))


for v in R:
    if not any(v):
        continue
    ps = points_of(steps_of(v))[:-1]
    # first simple sub-cycle: walk until a vertex repeats
    seen = {}; sub = None
    for i, p in enumerate(ps + [ps[0]]):
        if p in seen:
            sub = ps[seen[p]:i]; break
        seen[p] = i
    if sub is not None and len(sub) % 4 == 0:
        put([idx[p] for p in sub])
for (i, s, j, k) in T:
    st = steps_of(R[i]) + steps_of(R[j], s); ps = points_of(st)
    changed = True
    while changed:
        changed = False
        a = 0
        while a < len(st) - 1:
            x, y = st[a], st[a + 1]
            if x[0] == y[0] and x[1] == -y[1]:
                del st[a:a + 2]; del ps[a + 1:a + 3]; changed = True; a = max(a - 1, 0); continue
            if x[0] > y[0]:
                b = ps[a]; c1 = ps[a + 1]; c2 = ps[a + 2]
                new = tuple(b[t] + y[1] * U[y[0]][t] for t in range(dim))
                sq = [idx[b], idx[c1], idx[c2], idx[new]]
                if len(set(sq)) == 4:
                    put(sq); nsq += 1
                st[a], st[a + 1] = y, x; ps[a + 1] = new; changed = True
            a += 1
    assert st == steps_of(R[k])
print(f"{nsq} swaps, {len(cyc)} directed cycles listed")
json.dump({"D": W["D"], "units": U, "points": [list(p) for p in pts], "cycles": [list(c) for c in sorted(cyc)]},
          open(sys.argv[2], "w"))
