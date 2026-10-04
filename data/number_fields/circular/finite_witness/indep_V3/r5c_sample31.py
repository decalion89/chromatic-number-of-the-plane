"""Run grow.py (3,1) + minimise.py on a random sample of the 5-cycle seeds over Q(sqrt31) (denominator 80)."""
import json, random, subprocess, re, gzip, os, collections
import networkx as nx
os.chdir('sample31')
for f in ['grow.py', 'minimise.py']:
    subprocess.run(['cp', f'../repo_copy/{f}', f], check=True)
U = [tuple(u) for u in json.load(open('../seed31_cands.json'))['units']]
classes = json.load(open('../cycles31_classes.json'))
S = json.load(gzip.open('../repo_copy/witness_q31.json.gz', 'rt'))
GS = nx.Graph(); GS.add_edges_from(map(tuple, S['edges']))
SP = {tuple(p) for p in S['points']}
rng = random.Random(2026)
sample = rng.sample(range(len(classes)), 40)
stats = collections.Counter()
for k in sample:
    a, b, c, d = map(tuple, classes[k])
    pts = sorted({(0, 0, 0, 0)} | set(U) | {b, c})
    json.dump({'d': 31, 'D': 80, 'points': [list(p) for p in pts]}, open('seed.json', 'w'))
    out = subprocess.run(['python3', '-B', 'grow.py', 'seed.json', 'A', 'G', '2000', '200', '3', '1'],
                         capture_output=True, text=True).stdout
    m = re.search(r'UNSAT with (\d+) points', out)
    if not m:
        print(k, 'no UNSAT:', out.strip().splitlines()[-1]); stats['no UNSAT'] += 1; continue
    nG = int(m.group(1))
    subprocess.run(['python3', '-B', 'minimise.py', 'G.json', 'W.json'], capture_output=True, text=True)
    W = json.load(open('W.json'))
    G = nx.Graph(); G.add_nodes_from(range(len(W['points']))); G.add_edges_from(map(tuple, W['edges']))
    iso = nx.is_isomorphic(G, GS) if len(W['points']) == 9 else False
    # same points up to translation?
    WP = [tuple(p) for p in W['points']]
    trans = any({tuple(x - y for x, y in zip(w, t)) for w in WP} == SP for t in {tuple(x - y for x, y in zip(w, s)) for w in WP for s in SP})
    print(k, 'grown', nG, '-> witness', len(W['points']), 'vertices', G.number_of_edges(), 'edges; iso to q31 graph:', iso,
          '; stored points up to translation:', trans, flush=True)
    stats[(nG, len(W['points']), G.number_of_edges(), iso, trans)] += 1
print(stats)
