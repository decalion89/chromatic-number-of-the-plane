#!/usr/bin/env python3
"""Referee checker (written from scratch, shares no code with the repository).

Usage:  python3 referee_check.py W.json.gz OUTDIR [--nofix]

Checks, in exact integer arithmetic:
  (1) d is squarefree > 1 (so {1, sqrt d} is a Q-basis and the representation of a point is unique),
      the points are distinct, every listed edge has length exactly 1, and the edge list is
      exactly the set of all unordered pairs at distance 1 (all pairs compared);
  (2) the colouring is a (7,2)-colouring: c(y)-c(x) mod 7 in {2,3,4,5} on every edge;
  (3) the listed cycles are closed walks of H (and reports simplicity, lengths, lengths mod 7),
      and writes a CNF in a NEW encoding (colour-major numbering, one-sided non-tight indicators):
        x(v,k) = k*n + v + 1          "vertex v has colour k"          (k in 0..6)
        N(a)   = 7n + idx(a) + 1      "arc a=(u,w) of a listed cycle is NOT tight" (one direction only:
                                       N(a) -> not(x(u,k) and x(w,k+2)) for every k)
      clauses: at-least-one and at-most-one colour per vertex; for every edge {u,w} and colour k,
      not(x(u,k) and x(w,k+delta)) for delta in {0,1,6}; the N-implications; for every listed cycle
      the clause OR_a N(a); and (unless --nofix) the unit clause x(fixed_vertex, 0).
  A satisfying assignment gives a (7,2)-colouring in which every listed cycle has a non-tight arc
  (choose any true colour per vertex; at-most-one makes it unique anyway); conversely such a
  colouring c, rotated so that c(fixed_vertex)=0, satisfies the formula with N(a) := [a non-tight].
"""
import gzip, json, sys, os, hashlib
from collections import Counter, defaultdict, deque


def is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def squarefree(d):
    if d < 2:
        return False
    k = 2
    while k * k <= d:
        if d % (k * k) == 0:
            return False
        k += 1
    return True


def fail(msg):
    print("FAIL:", msg)
    sys.exit(1)


def main():
    path, outdir = sys.argv[1], sys.argv[2]
    nofix = '--nofix' in sys.argv[3:]
    name = os.path.basename(path).replace('.json.gz', '')
    with gzip.open(path, 'rt') as f:
        W = json.load(f)
    d, D = W['d'], W['denominator']
    P, E, col, cyc = W['points'], W['edges'], W['colouring'], W['cycles']
    fv = W.get('fixed_vertex')
    report = []

    def say(s):
        print(s)
        report.append(s)

    say(f"== {name}: d={d}, D={D}")
    # ---------------- (1) points and edges ----------------
    if not (is_int(d) and squarefree(d)):
        fail("d not a squarefree integer > 1")
    if not (is_int(D) and D > 0):
        fail("bad denominator")
    n = len(P)
    for p in P:
        if not (isinstance(p, list) and len(p) == 4 and all(is_int(t) for t in p)):
            fail(f"bad point {p}")
    if len(set(map(tuple, P))) != n:
        fail("repeated point")
    # gcd information (not needed for correctness, reported only)
    say(f"(1) {n} points, all distinct integer 4-tuples; d squarefree, so distinct tuples = distinct points")
    # all pairs at distance 1, exact
    unit = set()
    vecs = Counter()
    for i in range(n):
        a1, b1, c1, e1 = P[i]
        for j in range(i + 1, n):
            a2, b2, c2, e2 = P[j]
            A, B, C, Ee = a2 - a1, b2 - b1, c2 - c1, e2 - e1
            # ((A + B r)^2 + (C + Ee r)^2) / D^2 = 1 with r = sqrt d  <=>  rational and irrational parts
            if A * A + d * B * B + C * C + d * Ee * Ee == D * D and A * B + C * Ee == 0:
                unit.add((i, j))
                vecs[(A, B, C, Ee)] += 1
                vecs[(-A, -B, -C, -Ee)] += 1
    elist = []
    for e in E:
        if not (isinstance(e, list) and len(e) == 2 and all(is_int(t) for t in e)):
            fail(f"bad edge {e}")
        i, j = e
        if not (0 <= i < n and 0 <= j < n) or i == j:
            fail(f"edge out of range or loop {e}")
        elist.append((min(i, j), max(i, j)))
    eset = set(elist)
    if len(eset) != len(elist):
        fail("duplicate edge")
    bad = [e for e in eset if e not in unit]
    missing = [e for e in unit if e not in eset]
    say(f"    {len(elist)} listed edges; unit-distance pairs among all {n*(n-1)//2} pairs: {len(unit)}; "
        f"listed edges not at distance 1: {len(bad)}; unit pairs not listed: {len(missing)}")
    if bad or missing:
        fail("edge set != unit-distance pairs")
    say(f"    => every edge has length exactly 1 and H is the induced unit-distance graph on its points; "
        f"{len(vecs)} distinct unit vectors occur as edge differences")
    adj = defaultdict(set)
    for i, j in elist:
        adj[i].add(j)
        adj[j].add(i)
    deg = [len(adj[v]) for v in range(n)]
    # connectivity
    seen = {0}
    dq = deque([0])
    while dq:
        x = dq.popleft()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                dq.append(y)
    say(f"    min degree {min(deg)}, max degree {max(deg)}, connected: {len(seen) == n}")
    # ---------------- (2) colouring ----------------
    if len(col) != n or not all(is_int(c) and 0 <= c < 7 for c in col):
        fail("colouring has wrong length or values outside 0..6")
    badc = [(i, j) for i, j in elist if (col[j] - col[i]) % 7 not in (2, 3, 4, 5)]
    say(f"(2) colouring values in 0..6; edges violating the (7,2) condition: {len(badc)}")
    if badc:
        fail("not a (7,2)-colouring")
    # ---------------- (3) cycles ----------------
    lens = Counter()
    nonsimple = 0
    arcs = set()
    for C in cyc:
        if not (isinstance(C, list) and len(C) >= 3 and all(is_int(t) and 0 <= t < n for t in C)):
            fail(f"bad cycle {C}")
        m = len(C)
        for k in range(m):
            u, w = C[k], C[(k + 1) % m]
            if (min(u, w), max(u, w)) not in eset:
                fail(f"cycle step {u}->{w} is not an edge")
            arcs.add((u, w))
        if len(set(C)) != m:
            nonsimple += 1
        lens[m] += 1
    say(f"(3) {len(cyc)} cycles; every step is an edge of H; non-simple: {nonsimple}; "
        f"lengths {dict(sorted(lens.items()))}; lengths not divisible by 7: {sum(v for k, v in lens.items() if k % 7)}")
    distinct_cycles = len(set(tuple(C) for C in cyc))
    # canonical up to rotation
    def canon(C):
        m = len(C)
        return min(tuple(C[k:] + C[:k]) for k in range(m))
    distinct_rot = len(set(canon(C) for C in cyc))
    say(f"    distinct as sequences: {distinct_cycles}; distinct up to rotation: {distinct_rot}; "
        f"directed arcs used: {len(arcs)}")
    tight_given = sum(1 for C in cyc if all((col[C[(k + 1) % len(C)]] - col[C[k]]) % 7 == 2 for k in range(len(C))))
    say(f"    listed cycles tight under the stored colouring: {tight_given} (must be >= 1 if the claim holds)")
    if fv is not None:
        mx = max(deg)
        say(f"    fixed_vertex {fv}: degree {deg[fv]}; max degree {mx}; first vertex of max degree: {deg.index(mx)}")
    # ---------------- CNF ----------------
    def x(v, k):
        return k * n + v + 1
    arcl = sorted(arcs)
    aidx = {a: 7 * n + i + 1 for i, a in enumerate(arcl)}
    nv = 7 * n + len(arcl)
    cls = []
    for v in range(n):
        cls.append([x(v, k) for k in range(7)])
        for k in range(7):
            for l in range(k + 1, 7):
                cls.append([-x(v, k), -x(v, l)])
    for (u, w) in elist:
        for k in range(7):
            for dl in (0, 1, 6):
                cls.append([-x(u, k), -x(w, (k + dl) % 7)])
    for (u, w), N in aidx.items():
        for k in range(7):
            cls.append([-N, -x(u, k), -x(w, (k + 2) % 7)])
    for C in cyc:
        m = len(C)
        cls.append([aidx[(C[k], C[(k + 1) % m])] for k in range(m)])
    if fv is not None and not nofix:
        cls.append([x(fv, 0)])
    os.makedirs(outdir, exist_ok=True)
    suffix = '_nofix' if (nofix or fv is None) else ''
    out = os.path.join(outdir, f"{name}{suffix}.ref.cnf")
    txt = f"p cnf {nv} {len(cls)}\n" + ''.join(' '.join(map(str, c)) + ' 0\n' for c in cls)
    with open(out, 'w') as f:
        f.write(txt)
    h = hashlib.sha256(txt.encode()).hexdigest()
    say(f"    CNF {out}: {nv} variables, {len(cls)} clauses; sha256 {h}"
        + ("" if (fv is None or nofix) else f"; unit clause fixes vertex {fv} to colour 0"))
    with open(os.path.join(outdir, f"{name}{suffix}.ref.log"), 'w') as f:
        f.write('\n'.join(report) + '\n')


if __name__ == '__main__':
    main()
