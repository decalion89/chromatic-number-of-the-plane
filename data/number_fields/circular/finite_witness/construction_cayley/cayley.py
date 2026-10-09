"""cayley.py: the witnesses H4 (over Q(sqrt3, sqrt11)) and H4' (over Q(sqrt2, sqrt3)) restricted to the Cayley graph of
the unit vectors U of their construction, with chi_c = 4 still (every proper 4-colouring has a tight cycle).

  python3 cayley.py init WHICH OUT.json   WHICH = q3_11 (the 27 vectors of ../../at_four/q3_11.py) or q2_3 (the 60
            vectors zeta_24^j w^l, j mod 24, |l| <= 2, w = (1 + 2 sqrt(-2))/3, one of each pair +-u): the points of
            ../witness_WHICH.json.gz, the pairs that differ by an element of U u -U as edges, and the listed cycles; then,
            while kissat ($KISSAT) finds a proper 4-colouring with no tight listed cycle, every tight directed cycle found
            by a depth-first search in its tight arcs is added (with its four rotations of colours as clauses); stops
            when the formula is unsatisfiable (or exits with code 2 if a colouring has an acyclic tight digraph)
  python3 cayley.py cnf IN.json OUT.cnf    the formula of ../check_witness4.py
  python3 cayley.py reduce IN.json CORE OUT.json   keeps the cycles with a clause in the clausal core CORE (drat-trim -c)
  python3 cayley.py final IN.json NAME     NAME.json.gz and NAME.cnf.gz (reproducible bytes)
Exact integers and fractions only."""
import sys, os, json, gzip, subprocess
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'at_four'))
from check_witness4 import cnf_text, load

x = lambda v, k: 4 * v + k + 1


def units_q2_3(D):
    """the 120 unit vectors zeta_24^j w^l (j mod 24, |l| <= 2) of Q(sqrt2, sqrt3)^2, as integer 8-tuples over D"""
    def mul(p, q):      # Q(sqrt2, sqrt3) over (1, s2, s3, s6): s2 s3 = s6, s2 s6 = 2 s3, s3 s6 = 3 s2
        return (p[0] * q[0] + 2 * p[1] * q[1] + 3 * p[2] * q[2] + 6 * p[3] * q[3],
                p[0] * q[1] + p[1] * q[0] + 3 * (p[2] * q[3] + p[3] * q[2]),
                p[0] * q[2] + p[2] * q[0] + 2 * (p[1] * q[3] + p[3] * q[1]),
                p[0] * q[3] + p[3] * q[0] + p[1] * q[2] + p[2] * q[1])
    add = lambda p, q: tuple(s + t for s, t in zip(p, q))
    neg = lambda p: tuple(-s for s in p)
    cmul = lambda z, w: (add(mul(z[0], w[0]), neg(mul(z[1], w[1]))), add(mul(z[0], w[1]), mul(z[1], w[0])))
    f = lambda *t: tuple(Fr(s) for s in t)
    zeta = (f(0, Fr(1, 4), 0, Fr(1, 4)), f(0, Fr(-1, 4), 0, Fr(1, 4)))      # cos 15 = (s2 + s6)/4, sin 15 = (s6 - s2)/4
    w = (f(Fr(1, 3), 0, 0, 0), f(0, Fr(2, 3), 0, 0))                         # (1 + 2 sqrt(-2))/3
    wbar = (w[0], neg(w[1]))
    one = (f(1, 0, 0, 0), f(0, 0, 0, 0))
    powers = {0: one, 1: w, 2: cmul(w, w), -1: wbar, -2: cmul(wbar, wbar)}
    U, z = set(), one
    for j in range(24):
        for l, wl in powers.items():
            v = cmul(z, wl)
            t = [s * D for s in v[0] + v[1]]
            assert all(s.denominator == 1 for s in t)
            U.add(tuple(int(s) for s in t))
        z = cmul(z, zeta)
    assert len(U) == 120, len(U)
    return list(U)


def init(which, outp):
    W = load(os.path.join(HERE, '..', f'witness_{which}.json.gz'))
    if which == 'q3_11':
        from q3_11 import units_311
        assert W['denominator'] == 84
        U = [tuple(u) for u in units_311()]
    else:
        U = units_q2_3(W['denominator'])
    P = W['points']
    dirs = set(U) | {tuple(-t for t in u) for u in U}
    idx = {tuple(p): k for k, p in enumerate(P)}
    cay = set()
    for k, p in enumerate(P):
        for d in dirs:
            q = tuple(p[t] + d[t] for t in range(8))
            if q in idx:
                cay.add((min(k, idx[q]), max(k, idx[q])))
    E = sorted(cay)
    n = len(P); nv = 4 * n
    print(f'{which}: {n} points, {len(E)} of the {len(W["edges"])} unit pairs differ by an element of U u -U '
          f'({len(dirs)} vectors)', flush=True)
    clauses = [l for l in cnf_text(n, E, W['cycles'], W['fixed_vertex']).split('\n')[1:] if l.strip()]
    adj = [set() for _ in range(n)]
    for i, j in E:
        adj[i].add(j); adj[j].add(i)
    kissat = os.environ.get('KISSAT', 'kissat')
    added = []
    for it in range(200):
        cnf = outp + '.lazy.cnf'
        with open(cnf, 'w') as f:
            f.write(f'p cnf {nv} {len(clauses) + len(added)}\n' + '\n'.join(clauses) + '\n'
                    + ''.join(' '.join(map(str, c)) + ' 0\n' for c in added))
        r = subprocess.run([kissat, cnf], capture_output=True, text=True)
        if r.returncode == 20:
            print(f'round {it}: unsatisfiable with {len(added) // 4} added cycles', flush=True)
            break
        assert r.returncode == 10, r.returncode
        val = set()
        for l in r.stdout.split('\n'):
            if l.startswith('v '):
                val.update(int(t) for t in l.split()[1:] if int(t) > 0)
        c = [next(k for k in range(4) if x(v, k) in val) for v in range(n)]
        out = [[w for w in adj[u] if (c[w] - c[u]) % 4 == 1] for u in range(n)]
        state = [0] * n; cycles = []
        for s in range(n):
            if state[s]:
                continue
            stack = [(s, 0)]; path = [s]; onpath = {s: 0}; state[s] = 1
            while stack:
                u, i = stack[-1]
                if i < len(out[u]):
                    stack[-1] = (u, i + 1)
                    w = out[u][i]
                    if w in onpath:
                        cycles.append(path[onpath[w]:])
                    elif state[w] == 0:
                        state[w] = 1; onpath[w] = len(path); path.append(w); stack.append((w, 0))
                else:
                    stack.pop(); path.pop(); del onpath[u]; state[u] = 2
        if not cycles:
            print(f'round {it}: a proper 4-colouring without tight cycles; chi_c < 4', flush=True)
            sys.exit(2)
        for cyc in cycles:
            for k in range(4):
                added.append([-x(cyc[t], (k + t) % 4) for t in range(len(cyc))])
        print(f'round {it}: {len(cycles)} tight cycles added, {len(added) // 4} in all', flush=True)
    os.remove(cnf)
    newcyc = [[(-added[q][t] - 1) // 4 for t in range(len(added[q]))] for q in range(0, len(added), 4)]
    out = {k: W[k] for k in ('field', 'basis', 'denominator', 'points', 'colouring', 'fixed_vertex')}
    out['generators'] = [list(u) for u in sorted({min(tuple(u), tuple(-t for t in u)) for u in U})]
    out['edges'] = [list(e) for e in E]
    out['cycles'] = W['cycles'] + newcyc
    json.dump(out, open(outp, 'w'))
    print(f'{outp}: {len(out["generators"])} generators, {len(E)} edges, {len(out["cycles"])} cycles')


def reduce(inp, core, outp):
    W = json.load(open(inp))
    keep_cl = set()
    for l in open(core):
        if l.startswith('p') or l.startswith('c') or not l.strip():
            continue
        keep_cl.add(frozenset(int(t) for t in l.split()[:-1]))
    keep = [cyc for cyc in W['cycles']
            if any(frozenset(-x(cyc[t], (k + t) % 4) for t in range(len(cyc))) in keep_cl for k in range(4))]
    print(f'cycles {len(W["cycles"])} -> {len(keep)}')
    W['cycles'] = keep
    json.dump(W, open(outp, 'w'))


def final(inp, name):
    W = json.load(open(inp))
    out = {k: W[k] for k in ('field', 'basis', 'denominator', 'generators', 'points', 'edges', 'colouring', 'cycles',
                             'fixed_vertex')}
    with open(name + '.json.gz', 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as f:
        f.write(json.dumps(out).encode())
    with open(name + '.cnf.gz', 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as f:
        f.write(cnf_text(len(W['points']), W['edges'], W['cycles'], W['fixed_vertex']).encode())
    print(f'{name}.json.gz, {name}.cnf.gz')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'init':
        init(sys.argv[2], sys.argv[3])
    elif cmd == 'cnf':
        W = json.load(open(sys.argv[2]))
        open(sys.argv[3], 'w').write(cnf_text(len(W['points']), W['edges'], W['cycles'], W['fixed_vertex']))
    elif cmd == 'reduce':
        reduce(sys.argv[2], sys.argv[3], sys.argv[4])
    elif cmd == 'final':
        final(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
