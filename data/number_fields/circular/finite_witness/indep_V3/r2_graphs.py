"""Referee checks 2-5 on the two nine-point witnesses: structure, colourings, tight cycles, the hand proof,
K_{8/3}, and the criticality certificates.  Own code; networkx only for isomorphism tests."""
import gzip, json, sys, itertools
from fractions import Fraction as Fr
import networkx as nx
from pysat.solvers import Solver

BASE = sys.argv[1]


def load(name):
    with gzip.open(f'{BASE}/{name}', 'rt') as f:
        return json.load(f)


def ek(a, b):
    return (a, b) if a < b else (b, a)


def has_directed_cycle(n, arcs):
    """Kahn's algorithm (different from the DFS of check_small.py)."""
    indeg = [0] * n; out = [[] for _ in range(n)]
    for a, b in arcs:
        out[a].append(b); indeg[b] += 1
    queue = [v for v in range(n) if indeg[v] == 0]; seen = 0; order = []
    while queue:
        v = queue.pop(); seen += 1; order.append(v)
        for w in out[v]:
            indeg[w] -= 1
            if indeg[w] == 0:
                queue.append(w)
    return seen < n, order


def tight_arcs(c, E, p, q, skip=()):
    arcs = []
    for a, b in E:
        if a in skip or b in skip:
            continue
        if (c[b] - c[a]) % p == q:
            arcs.append((a, b))
        if (c[a] - c[b]) % p == q:
            arcs.append((b, a))
    return arcs


def homs(n, E, p, q, fix0=True):
    """All homomorphisms to K_{p/q} (vertex 0 coloured 0 if fix0) by exhaustive product (small n only)."""
    res = []
    rng = [range(p)] * n
    if fix0:
        rng = [range(1)] + [range(p)] * (n - 1)
    for c in itertools.product(*rng):
        if all(q <= (c[b] - c[a]) % p <= p - q for a, b in E):
            res.append(c)
    return res


def sat_hom(n, E, p, q):
    """Is there a homomorphism to K_{p/q}?  SAT encoding with exactly-one colour per vertex."""
    x = lambda v, k: v * p + k + 1
    s = Solver(name='minisat22')
    for v in range(n):
        s.add_clause([x(v, k) for k in range(p)])
        for k1, k2 in itertools.combinations(range(p), 2):
            s.add_clause([-x(v, k1), -x(v, k2)])
    for a, b in E:
        for k1 in range(p):
            for k2 in range(p):
                if not (q <= (k2 - k1) % p <= p - q):
                    s.add_clause([-x(a, k1), -x(b, k2)])
    r = s.solve(); s.delete()
    return r


def kpq(p, q):
    G = nx.Graph(); G.add_nodes_from(range(p))
    G.add_edges_from((i, j) for i in range(p) for j in range(i + 1, p) if q <= (j - i) % p <= p - q)
    return G


# ---------- claim 4: K_{8/3} is the Wagner graph --------------------------------------------------------------
K83 = kpq(8, 3)
Wag = nx.cycle_graph(8); Wag.add_edges_from((i, i + 4) for i in range(4))
seq = [0, 3, 6, 1, 4, 7, 2, 5]
cyc_ok = all(K83.has_edge(seq[i], seq[(i + 1) % 8]) for i in range(8))
opp_ok = all(K83.has_edge(seq[i], seq[i + 4]) for i in range(4))
print('K_{8/3}: edges', K83.number_of_edges(), '; 0,3,6,1,4,7,2,5 is an 8-cycle:', cyc_ok,
      '; its other edges join opposite vertices:', opp_ok and K83.number_of_edges() == 12,
      '; isomorphic to the Wagner graph:', nx.is_isomorphic(K83, Wag))
fr = sorted({Fr(a, b) for a in range(1, 10) for b in range(1, a + 1)})
below3 = [f for f in fr if f < 3]
print('largest fraction < 3 with numerator <= 9:', max(below3), '; next fraction >= 3 with numerator <= 9:',
      min(f for f in fr if f >= 3))

# ---------- the two witnesses ----------------------------------------------------------------------------------
results = {}
for name in ['witness_q7.json.gz', 'witness_q31.json.gz']:
    W = load(name); n = len(W['points']); p, q = W['p'], W['q']
    E = [tuple(e) for e in W['edges']]
    Es = {ek(*e) for e in E}
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
    print(f'\n== {name}: n = {n}, edges = {len(E)}, degrees = {[G.degree(v) for v in range(n)]}')
    # colouring stored
    col = W['colouring']
    print('stored colouring proper:', all((col[a] - col[b]) % 3 != 0 for a, b in E))
    # all proper 3-colourings, all 3^9 maps
    allc = [c for c in itertools.product(range(3), repeat=n) if all(c[a] != c[b] for a, b in E)]
    cyc = W['cycles']
    for cy in cyc:
        assert len(set(cy)) == len(cy) and all(ek(cy[i], cy[(i + 1) % len(cy)]) in Es for i in range(len(cy)))
    print('listed cycles are simple cycles along edges; lengths', sorted(len(c) for c in cyc))
    n_tc = n_listed = 0
    for c in allc:
        cy_, _ = has_directed_cycle(n, tight_arcs(c, E, 3, 1))
        n_tc += cy_
        n_listed += any(all((c[cy[(i + 1) % len(cy)]] - c[cy[i]]) % 3 == 1 for i in range(len(cy))) for cy in cyc)
    print(f'proper 3-colourings: {len(allc)}; with a tight cycle: {n_tc}; with a tight listed cycle: {n_listed}')
    # homomorphism to K_{8/3}: exhaustive over 8^8 is too slow in Python; use SAT and a backtracking search
    print('hom to K_{8/3} (SAT):', sat_hom(n, E, 8, 3), '; hom to K_3 (SAT):', sat_hom(n, E, 3, 1))
    # chi_c directly: for all fractions p'/q' < 3 with p' <= 30, no homomorphism
    bad = [f for f in sorted({Fr(a, b) for a in range(2, 31) for b in range(1, a)}) if 2 <= f < 3
           and sat_hom(n, E, f.numerator, f.denominator)]
    print('fractions 2 <= p/q < 3 with p <= 30 admitting a homomorphism:', bad)
    # criticality certificates
    crit = W['critical_colourings']
    allok = True; maxfr = []
    for v in range(n):
        c = crit[v]
        okc = (len(c) == n and c[v] == -1 and all(0 <= c[u] < 3 for u in range(n) if u != v)
               and all(c[a] != c[b] for a, b in E if v not in (a, b)))
        arcs = tight_arcs(c, E, 3, 1, skip=(v,))
        # Kahn on the graph without v: relabel
        others = [u for u in range(n) if u != v]; idx = {u: i for i, u in enumerate(others)}
        cyc_, order = has_directed_cycle(n - 1, [(idx[a], idx[b]) for a, b in arcs])
        pos = {others[o]: k for k, o in enumerate(order)}
        N = n - 1
        f = {u: N * c[u] + pos[u] for u in others}
        P2, Q2 = 3 * N, N + 1     # K_{3N/(N+1)} = K_{24/9} = K_{8/3}
        hom_ok = all(Q2 <= (f[b] - f[a]) % P2 <= P2 - Q2 for a, b in E if v not in (a, b))
        allok &= okc and not cyc_ and hom_ok
        # and directly: H - v maps to K_{8/3} (SAT, independent of the certificate)
        Ev = [(idx[a], idx[b]) for a, b in E if v not in (a, b)]
        maxfr.append(sat_hom(n - 1, Ev, 8, 3))
    print('critical certificates: all proper, acyclic tight digraphs, N*c+pos is a hom to K_{24/9}=K_{8/3}:', allok)
    print('each H - v maps to K_{8/3} (SAT, direct):', all(maxfr))
    results[name] = (G, allc)

# ---------- claim 2: structure -------------------------------------------------------------------------------
G7, allc7 = results['witness_q7.json.gz']
G31, allc31 = results['witness_q31.json.gz']
P = {f'P{i}': i for i in range(8)}; P['S'] = 8
# Wagner with chord P3P7 subdivided by S
WS = nx.Graph()
WS.add_edges_from((i, (i + 1) % 8) for i in range(8))
WS.add_edges_from([(0, 4), (1, 5), (2, 6), (3, 8), (8, 7)])
print('\nH7 equals (identically) the 8-cycle P0..P7 + chords P0P4, P1P5, P2P6 + path P3 S P7:',
      set(map(lambda e: ek(*e), G7.edges)) == set(map(lambda e: ek(*e), WS.edges)))
for ch in range(4):
    Wc = nx.Graph(); Wc.add_edges_from((i, (i + 1) % 8) for i in range(8))
    Wc.add_edges_from((i, i + 4) for i in range(4) if i != ch); Wc.add_edges_from([(ch, 8), (8, ch + 4)])
    print(f'  Wagner with chord {ch}-{ch + 4} subdivided isomorphic to H7:', nx.is_isomorphic(Wc, G7))


def M_graph(v, m):
    """hexagon v[0..5], m[j] joined to v[j], v[j+3]"""
    G = nx.Graph()
    G.add_edges_from((v[i], v[(i + 1) % 6]) for i in range(6))
    G.add_edges_from((m[j], v[j]) for j in range(3)); G.add_edges_from((m[j], v[j + 3]) for j in range(3))
    return G


v7 = [P['P0'], P['P4'], P['P3'], P['P2'], P['P6'], P['P7']]; m7 = [P['P1'], P['P5'], P['S']]
M7 = M_graph(v7, m7)
H7b = M7.copy(); H7b.add_edge(m7[0], m7[1])
print('H7 = M + m0m1 with v = P0,P4,P3,P2,P6,P7 and m = P1,P5,S (identical edge sets):',
      {ek(*e) for e in H7b.edges} == {ek(*e) for e in G7.edges})
v31 = [0, 1, 8, 7, 2, 6]; m31 = [3, 4, 5]
M31 = M_graph(v31, m31)
Q31 = M31.copy(); Q31.add_edges_from([(3, 4), (4, 5)])
print('witness_q31 = M + m0m1 + m1m2 with v = 0,1,8,7,2,6, m = 3,4,5 (identical edge sets):',
      {ek(*e) for e in Q31.edges} == {ek(*e) for e in G31.edges})
print('|Aut| of M, H7, Q31 graph:', [sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(H, H).isomorphisms_iter())
                                    for H in (M7, G7, G31)])
Mstd = M_graph(list(range(6)), [6, 7, 8])
print('M + any one midpoint edge isomorphic to H7:', all(nx.is_isomorphic(nx.compose(Mstd, nx.Graph([e])), G7)
                                                        for e in [(6, 7), (7, 8), (6, 8)]))
print('M + any two midpoint edges isomorphic to the q31 graph:',
      all(nx.is_isomorphic(nx.compose(Mstd, nx.Graph(es)), G31) for es in [[(6, 7), (7, 8)], [(6, 7), (6, 8)], [(6, 8), (7, 8)]]))
print('M is K_{3,3} with a perfect matching subdivided; hexagon + long diagonals is K_{3,3}:',
      nx.is_isomorphic(nx.compose(nx.cycle_graph(6), nx.Graph([(0, 3), (1, 4), (2, 5)])), nx.complete_bipartite_graph(3, 3)))

# ---------- claim 3: the hand proof on M = H7 - P1P5 ----------------------------------------------------------
v, m = v7, m7
V = lambda k: v[k % 6]
Mm = lambda k: m[k % 3]
Z = {k: [V(k), V(k + 1), V(k + 2), V(k + 3), Mm(k)] for k in range(7)}
C = {k: [V(k), Mm(k), V(k + 3), V(k + 4), Mm(k + 1), V(k + 1)] for k in range(6)}
HEX = [V(k) for k in range(6)]


def chain(walk):
    ch = {}
    for i in range(len(walk)):
        a, b = walk[i], walk[(i + 1) % len(walk)]
        if a < b:
            ch[(a, b)] = ch.get((a, b), 0) + 1
        else:
            ch[(b, a)] = ch.get((b, a), 0) - 1
    return {k: x for k, x in ch.items() if x}


def add(c1, c2, s=1):
    r = dict(c1)
    for k, x in c2.items():
        r[k] = r.get(k, 0) + s * x
    return {k: x for k, x in r.items() if x}


Mset = {ek(*e) for e in M7.edges}
print('\nhand proof: every Z_k is a 5-cycle of M:',
      all(len(set(Z[k])) == 5 and all(ek(Z[k][i], Z[k][(i + 1) % 5]) in Mset for i in range(5)) for k in range(6)))
print('every C_k is a 6-cycle of M:',
      all(len(set(C[k])) == 6 and all(ek(C[k][i], C[k][(i + 1) % 6]) in Mset for i in range(6)) for k in range(6)))
print('Z_{k+1} - Z_k = C_k as 1-chains, k = 0..5:', all(add(chain(Z[k + 1]), chain(Z[k]), -1) == chain(C[k]) for k in range(6)))
print('Z_6 = Z_0:', Z[6] == Z[0], '; Z_0 + Z_3 = hexagon:', add(chain(Z[0]), chain(Z[3])) == chain(HEX))
print('C_{k+3} = -C_k (same 6-cycle reversed):', all(chain(C[k + 3]) == {e: -x for e, x in chain(C[k]).items()} for k in range(3)))
EM = list(Mset)
allM = [c for c in itertools.product(range(3), repeat=9) if all(c[a] != c[b] for a, b in EM)]


def dsum(c, walk):
    s = 0
    for i in range(len(walk)):
        a, b = walk[i], walk[(i + 1) % len(walk)]
        s += 1 if (c[b] - c[a]) % 3 == 1 else -1
    return s


okz = okc = okid = okconc = True; how = {'C': 0, 'hex': 0}
for c in allM:
    zs = [dsum(c, Z[k]) for k in range(6)]; cs = [dsum(c, C[k]) for k in range(6)]; hx = dsum(c, HEX)
    okz &= all(z in (3, -3) for z in zs); okc &= all(x in (0, 6, -6) for x in cs)
    okid &= all(zs[(k + 1) % 6] - zs[k] == cs[k] for k in range(6)) and zs[0] + zs[3] == hx
    if all(x == 0 for x in cs):
        okconc &= hx in (6, -6); how['hex'] += 1
    else:
        how['C'] += 1
print(f'M: {len(allM)} proper 3-colourings; delta(Z_k) = +-3: {okz}; delta(C_k) in {{0, +-6}}: {okc}; '
      f'identities hold numerically: {okid}; if no C_k tight then hexagon delta = +-6: {okconc}; split {how}')
tc = sum(has_directed_cycle(9, tight_arcs(c, EM, 3, 1))[0] for c in allM)
print('M: proper 3-colourings with a tight cycle (Kahn):', tc, 'of', len(allM))
print('M maps to K_{8/3} (SAT):', sat_hom(9, EM, 8, 3), '; M is 3-colourable:', len(allM) > 0)
