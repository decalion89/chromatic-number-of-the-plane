#!/usr/bin/env python3
"""Semantic validation of the repository's STORED formula W.cnf.gz (referee code, independent of the repository).

Usage:  python3 validate_stored.py W.json.gz W.cnf.gz W.drat.xz OUTDIR

Writes OUTDIR/W.stored.cnf and OUTDIR/W.stored.drat (decompressed) and reports their sha256 and sizes.
Then checks that EVERY clause of the stored CNF is satisfied by the 'intended assignment' of any (7,2)-colouring c of H
in which every listed cycle has a non-tight arc (rotated so that the fixed vertex has colour 0), under the
interpretation x(v,k) = 7v + k + 1 and t = [arc(t) is tight], where arc(t) is read off the clauses that mention t.
Accepted clause shapes (anything else is reported):
  ALO   {x(v,0..6)} positive
  EDGE  {-x(i,k), -x(j,k')}, {i,j} an edge of H (as verified by referee_check.py: the induced unit-distance graph),
        k'-k mod 7 in {0,1,6}
  TDEF  {-x(a,k), -x(b,k+2), +t}  with ab an edge; all TDEF clauses of t name the same arc
  CYC   {-t_1, ..., -t_m}: the arcs of the t's form exactly the arc set of a listed cycle (hence contain a directed cycle)
  UNIT  {x(v0,0)}: at most one such clause, v0 = fixed_vertex
If all clauses pass, UNSAT of the stored CNF implies: every (7,2)-colouring of H has a tight listed cycle.
"""
import gzip, json, sys, os, hashlib, lzma
from collections import Counter, defaultdict


if not __debug__:
    sys.exit('run without -O: this script relies on assert statements')


def main():
    wpath, cpath, dpath, outdir = sys.argv[1:5]
    name = os.path.basename(wpath).replace('.json.gz', '')
    with gzip.open(wpath, 'rt') as f:
        W = json.load(f)
    n = len(W['points'])
    E = {(min(i, j), max(i, j)) for i, j in W['edges']}
    cycarcs = set()
    for C in W['cycles']:
        m = len(C)
        cycarcs.add(frozenset((C[k], C[(k + 1) % m]) for k in range(m)))
    fv = W.get('fixed_vertex')
    with gzip.open(cpath, 'rt') as f:
        txt = f.read()
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, name + '.stored.cnf'), 'w') as f:
        f.write(txt)
    print(f"{name}: stored CNF {len(txt)} bytes uncompressed, sha256 {hashlib.sha256(txt.encode()).hexdigest()}")
    with lzma.open(dpath, 'rb') as f:
        dr = f.read()
    with open(os.path.join(outdir, name + '.stored.drat'), 'wb') as f:
        f.write(dr)
    fmt = 'binary' if (dr[:1] in (b'a', b'd') and bytes([0]) in dr[:200]) else 'text?'
    print(f"{name}: stored DRAT {len(dr)} bytes uncompressed, sha256 {hashlib.sha256(dr).hexdigest()}, {fmt} format")
    lines = [l for l in txt.split('\n') if l and not l.startswith('c')]
    assert lines[0].startswith('p cnf')
    nv, nc = map(int, lines[0].split()[2:4])
    cls = []
    for l in lines[1:]:
        t = list(map(int, l.split()))
        assert t[-1] == 0 and 0 not in t[:-1]
        cls.append(t[:-1])
    assert len(cls) == nc, (len(cls), nc)
    X = 7 * n

    def xv(lit):
        v = abs(lit) - 1
        return v // 7, v % 7

    arc = {}
    kinds = Counter()
    bad = []
    units = []
    cycle_clauses = []
    for c in cls:
        if all(l > 0 and l <= X for l in c) and len(c) == 7 and len({xv(l)[0] for l in c}) == 1 \
                and {xv(l)[1] for l in c} == set(range(7)):
            kinds['ALO'] += 1
        elif len(c) == 2 and all(l < 0 and -l <= X for l in c):
            (i, k), (j, k2) = xv(c[0]), xv(c[1])
            if (min(i, j), max(i, j)) in E and (k2 - k) % 7 in (0, 1, 6):
                kinds['EDGE'] += 1
            else:
                bad.append(c)
        elif len(c) == 3 and sum(1 for l in c if l > X) == 1:
            t = [l for l in c if l > X][0]
            xs = [l for l in c if abs(l) <= X]
            if len(xs) == 2 and all(l < 0 for l in xs):
                (a, k), (b, k2) = xv(xs[0]), xv(xs[1])
                ok = False
                for (p, kp), (q, kq) in (((a, k), (b, k2)), ((b, k2), (a, k))):
                    if (kq - kp) % 7 == 2 and (min(p, q), max(p, q)) in E:
                        if arc.setdefault(t, (p, q)) == (p, q):
                            ok = True
                            break
                if ok:
                    kinds['TDEF'] += 1
                else:
                    bad.append(c)
            else:
                bad.append(c)
        elif all(l < 0 and -l > X for l in c):
            cycle_clauses.append(c)
        elif len(c) == 1 and 0 < c[0] <= X:
            units.append(c)
        else:
            bad.append(c)
    # TDEF arcs: each t has 7 TDEF clauses (one per k) -- informational
    for c in cycle_clauses:
        if all(-l in arc for l in c) and frozenset(arc[-l] for l in c) in cycarcs:
            kinds['CYC'] += 1
        else:
            bad.append(c)
    if len(units) == 1 and xv(units[0][0]) == (fv, 0):
        kinds['UNIT'] += 1
    elif units:
        bad += units
    covered = set()
    for c in cycle_clauses:
        covered.add(frozenset(arc[-l] for l in c))
    print(f"{name}: {nv} variables, {nc} clauses: {dict(kinds)}; unexplained clauses: {len(bad)}; "
          f"listed cycles covered by CYC clauses: {len(covered & cycarcs)} of {len(cycarcs)}")
    if bad:
        print("   first unexplained:", bad[:5])


if __name__ == '__main__':
    main()
