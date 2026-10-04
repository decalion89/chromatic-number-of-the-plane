#!/usr/bin/env python3
"""Clause-by-clause semantic validation of a CNF in the referee's E1 layout (also used by referee_cegar.py).

Usage:  python3 validate_e1.py W.json.gz F.cnf

Interpretation: x(v,k) = k*n + v + 1 ("v has colour k"); a variable N > 7n is an arc indicator whose arc (a,b) is read
off its clauses (-N | -x(a,k) | -x(b,k+2)), meaning N -> "a->b is not tight".
Hypothesis H0: c is a (7,2)-colouring of H whose tight digraph is ACYCLIC, rotated so that the unit clause's vertex
has colour 0.  Every clause must be satisfied by x(v,k) = [c(v)=k], N = [arc(N) not tight] under H0:
  ALO  {x(v,0..6)};  AMO {-x(v,k), -x(v,l)};  EDGE {-x(i,k), -x(j,k')} with ij an edge, k'-k in {0,1,6} mod 7;
  NDEF {-N, -x(a,k), -x(b,k+2)} with ab an edge (consistent arc per N);
  CYC  {N_1, ..., N_m} whose arcs are edges of H and CONTAIN A DIRECTED CYCLE (so under H0 one of them is not tight);
  UNIT {x(v0,0)} (at most one).
So if F is UNSAT, no (7,2)-colouring of H has an acyclic tight digraph, i.e. every (7,2)-colouring has a tight cycle.
"""
import gzip, json, sys
from collections import Counter, defaultdict

if not __debug__:
    sys.exit('run without -O: this script relies on assert statements')


def has_directed_cycle(arcs):
    out = defaultdict(list)
    nodes = set()
    for a, b in arcs:
        out[a].append(b); nodes.add(a); nodes.add(b)
    color = {}
    for s in nodes:
        if s in color:
            continue
        stack = [(s, iter(out[s]))]
        color[s] = 1
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                color[v] = 2; stack.pop()
            elif color.get(nxt) == 1:
                return True
            elif nxt not in color:
                color[nxt] = 1; stack.append((nxt, iter(out[nxt])))
    return False


def is_simple_cycle(arcs):
    outd, ind = Counter(a for a, b in arcs), Counter(b for a, b in arcs)
    nodes = set(outd) | set(ind)
    if any(outd[v] != 1 or ind[v] != 1 for v in nodes):
        return False
    succ = dict(arcs)
    s = next(iter(nodes)); v = succ[s]; k = 1
    while v != s:
        v = succ[v]; k += 1
    return k == len(nodes)


def main():
    with gzip.open(sys.argv[1], 'rt') as f:
        W = json.load(f)
    n = len(W['points'])
    E = {(min(i, j), max(i, j)) for i, j in W['edges']}
    lines = [l for l in open(sys.argv[2]).read().split('\n') if l and not l.startswith('c')]
    nv, nc = map(int, lines[0].split()[2:4])
    cls = [list(map(int, l.split()))[:-1] for l in lines[1:]]
    assert len(cls) == nc
    X = 7 * n
    xv = lambda lit: ((abs(lit) - 1) % n, (abs(lit) - 1) // n)   # (vertex, colour)
    arc = {}
    kinds = Counter(); bad = []; units = []; cyc = []; simple = 0
    for c in cls:
        if len(c) == 7 and all(0 < l <= X for l in c) and len({xv(l)[0] for l in c}) == 1 \
                and {xv(l)[1] for l in c} == set(range(7)):
            kinds['ALO'] += 1
        elif len(c) == 2 and all(l < 0 and -l <= X for l in c):
            (i, k), (j, k2) = xv(c[0]), xv(c[1])
            if i == j and k != k2:
                kinds['AMO'] += 1
            elif (min(i, j), max(i, j)) in E and (k2 - k) % 7 in (0, 1, 6):
                kinds['EDGE'] += 1
            else:
                bad.append(c)
        elif len(c) == 3 and all(l < 0 for l in c) and sum(1 for l in c if -l > X) == 1:
            N = -[l for l in c if -l > X][0]
            (a, k), (b, k2) = [xv(l) for l in c if -l <= X]
            ok = False
            for (p, kp), (q, kq) in (((a, k), (b, k2)), ((b, k2), (a, k))):
                if (kq - kp) % 7 == 2 and (min(p, q), max(p, q)) in E and arc.setdefault(N, (p, q)) == (p, q):
                    ok = True; break
            if ok:
                kinds['NDEF'] += 1
            else:
                bad.append(c)
        elif all(l > X for l in c):
            cyc.append(c)
        elif len(c) == 1 and 0 < c[0] <= X:
            units.append(c)
        else:
            bad.append(c)
    for c in cyc:
        arcs = [arc.get(l) for l in c]
        if None not in arcs and has_directed_cycle(arcs):
            kinds['CYC'] += 1
            simple += is_simple_cycle(arcs)
        else:
            bad.append(c)
    if len(units) == 1 and xv(units[0][0])[1] == 0:
        kinds['UNIT'] += 1
    elif units:
        bad += units
    print(f"{sys.argv[2].split('/')[-1]}: {nv} vars, {nc} clauses: {dict(kinds)}; CYC clauses that are exactly one simple "
          f"directed cycle: {simple}; unexplained clauses: {len(bad)}" + (f"; e.g. {bad[:3]}" if bad else ""))


if __name__ == '__main__':
    main()
