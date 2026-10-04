"""Regenerate the 9-vertex graphs with 12 and 13 edges over Q(sqrt31) from sampled seeds 47, 2058, 914, 330, and
verify them exactly (own arithmetic): induced unit-distance graphs, isomorphic to M / H7, chi_c = 3."""
import json, subprocess, os, itertools, gzip
from fractions import Fraction as Fr
import networkx as nx
os.chdir('sample31')
U = [tuple(u) for u in json.load(open('../seed31_cands.json'))['units']]
classes = json.load(open('../cycles31_classes.json'))
d, D = 31, 80
def unit(p, q):
    a, b, c, e = (x - y for x, y in zip(p, q))
    return a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0
M = nx.cycle_graph(6)
for j in range(3):
    M.add_edges_from([(j, 6 + j), (6 + j, j + 3)])
H7 = M.copy(); H7.add_edge(6, 7)
def homs_exist(G, p, q):
    from pysat.solvers import Solver
    n = G.number_of_nodes(); x = lambda v, k: v * p + k + 1
    with Solver(name='minisat22') as s:
        for v in range(n):
            s.add_clause([x(v, k) for k in range(p)])
        for a, b in G.edges:
            for k1 in range(p):
                for k2 in range(p):
                    if not (q <= (k2 - k1) % p <= p - q):
                        s.add_clause([-x(a, k1), -x(b, k2)])
        return s.solve()
for k in [47, 2058, 914, 330]:
    a, b, c, dd = map(tuple, classes[k])
    pts = sorted({(0, 0, 0, 0)} | set(U) | {b, c})
    json.dump({'d': 31, 'D': 80, 'points': [list(p) for p in pts]}, open(f'seed_{k}.json', 'w'))
    subprocess.run(['python3', '-B', 'grow.py', f'seed_{k}.json', 'A', f'G_{k}', '2000', '200', '3', '1'], capture_output=True, check=True)
    subprocess.run(['python3', '-B', 'minimise.py', f'G_{k}.json', f'W_{k}.json'], capture_output=True, check=True)
    subprocess.run(['cp', '../repo_copy/critical.py', 'critical.py'], check=True)
    subprocess.run(['python3', '-B', 'critical.py', f'W_{k}.json', f'C_{k}.json'], capture_output=True, check=True)
    W = json.load(open(f'W_{k}.json')); P = [tuple(p) for p in W['points']]; n = len(P)
    units = {(i, j) for i, j in itertools.combinations(range(n), 2) if unit(P[i], P[j])}
    E = {tuple(sorted(e)) for e in W['edges']}
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
    print(f'seed {k}: {n} points, distinct {len(set(P)) == n}, {len(E)} edges, induced (edges == all unit pairs): {units == E}; '
          f'iso M: {nx.is_isomorphic(G, M)}, iso H7: {nx.is_isomorphic(G, H7)}; hom to K_8/3: {homs_exist(G, 8, 3)}; '
          f'3-colourable: {homs_exist(G, 3, 1)}')
    for p in P:
        print('    ', p)
    C = json.load(open(f'C_{k}.json'))
    W['critical_colourings'] = C['critical_colourings']
    json.dump(W, open(f'W_{k}_cert.json', 'w'))
    r = subprocess.run(['python3', '-B', '../repo_copy/check_small.py', f'W_{k}_cert.json'], capture_output=True, text=True)
    print('   check_small:', r.stdout.strip().replace('\n', ' | ')[:400], r.stderr.strip()[:200])
