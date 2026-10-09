"""Random 5-cycle seeds over Q(sqrt7), denominator 160 (two extra points), through grow.py (3,1) and minimise.py."""
import json, random, subprocess, os, re, collections
import networkx as nx
os.chdir('sample7')
for f in ['grow.py', 'minimise.py']:
    subprocess.run(['cp', f'../repo_copy/{f}', f], check=True)
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
print('directed 5-cycles through 0 with two extra points, denominator 160:', len(cycles), flush=True)
M = nx.cycle_graph(6)
for j in range(3):
    M.add_edges_from([(j, 6 + j), (6 + j, j + 3)])
rng = random.Random(7); stats = collections.Counter()
for k in rng.sample(range(len(cycles)), 25):
    a, b, c, dd = cycles[k]
    seed = sorted({zero} | Us | {b, c})
    json.dump({'d': 7, 'D': 160, 'points': [list(p) for p in seed]}, open('seed.json', 'w'))
    out = subprocess.run(['python3', '-B', 'grow.py', 'seed.json', 'A', 'G', '2000', '200', '3', '1'], capture_output=True, text=True).stdout
    m = re.search(r'UNSAT with (\d+) points', out)
    if not m:
        print(k, 'no UNSAT'); continue
    subprocess.run(['python3', '-B', 'minimise.py', 'G.json', 'W.json'], capture_output=True, text=True)
    W = json.load(open('W.json')); G = nx.Graph(); G.add_nodes_from(range(len(W['points']))); G.add_edges_from(map(tuple, W['edges']))
    r = (int(m.group(1)), G.number_of_nodes(), G.number_of_edges(), nx.is_isomorphic(G, M))
    print(k, 'grown %d -> %d vertices, %d edges, iso M: %s' % r, flush=True)
    stats[r[1:]] += 1
    if r[3]:
        subprocess.run(['cp', 'W.json', f'W_M_{k}.json'])
print(stats)
