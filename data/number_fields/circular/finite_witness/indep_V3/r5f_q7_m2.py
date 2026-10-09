"""Regenerate the 14-edge nine-point graph over Q(sqrt7) from sampled seed 21222 and check it exactly."""
import json, subprocess, os, itertools, gzip
import networkx as nx
os.chdir('sample7')
g = json.load(open('../repo_copy/q7_seed.json'))
D, d = 160, 7
pts = [tuple(p) for p in g['points']]
U = [p for p in pts if p[0]**2 + d*p[1]**2 + p[2]**2 + d*p[3]**2 == D*D and p[0]*p[1] + p[2]*p[3] == 0]
Us = set(U); add = lambda p, q: tuple(x + y for x, y in zip(p, q)); zero = (0, 0, 0, 0)
nb = {}
for a in U:
    for u in U:
        b = add(a, u)
        if b != zero and b not in Us:
            nb.setdefault(b, set()).add(a)
cycles = []
for b, As in nb.items():
    for u in U:
        c = add(b, u)
        if c in nb and c != b:
            for a in As:
                for dd in nb[c]:
                    if a != dd:
                        cycles.append((a, b, c, dd))
a, b, c, dd = cycles[21222]
json.dump({'d': 7, 'D': 160, 'points': [list(p) for p in sorted({zero} | Us | {b, c})]}, open('seed_21222.json', 'w'))
subprocess.run(['python3', '-B', 'grow.py', 'seed_21222.json', 'A', 'G_21222', '2000', '200', '3', '1'], capture_output=True, check=True)
subprocess.run(['python3', '-B', 'minimise.py', 'G_21222.json', 'W_21222.json'], capture_output=True, check=True)
subprocess.run(['cp', '../repo_copy/critical.py', 'critical.py'], check=True)
subprocess.run(['python3', '-B', 'critical.py', 'W_21222.json', 'C_21222.json'], capture_output=True, check=True)
W = json.load(open('W_21222.json')); C = json.load(open('C_21222.json')); W['critical_colourings'] = C['critical_colourings']
json.dump(W, open('W_21222_cert.json', 'w'))
P = [tuple(p) for p in W['points']]; n = len(P)
unit = lambda p, q: (lambda a, b, c, e: a*a + d*b*b + c*c + d*e*e == D*D and a*b + c*e == 0)(*(x - y for x, y in zip(p, q)))
units = {(i, j) for i, j in itertools.combinations(range(n), 2) if unit(P[i], P[j])}
E = {tuple(sorted(e)) for e in W['edges']}
G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
S = json.load(gzip.open('../repo_copy/witness_q31.json.gz', 'rt')); G31 = nx.Graph(); G31.add_edges_from(map(tuple, S['edges']))
print('Q(sqrt7) seed 21222:', n, 'points; distinct', len(set(P)) == n, '; edges', len(E), '; induced:', units == E,
      '; isomorphic to M + m0m1 + m1m2 (the q31 graph):', nx.is_isomorphic(G, G31))
for p in P:
    print('   ', p, ' = ((%d + %d sqrt7)/160, (%d + %d sqrt7)/160)' % p)
r = subprocess.run(['python3', '-B', '../repo_copy/check_small.py', 'W_21222_cert.json'], capture_output=True, text=True)
print(r.stdout, r.stderr)
