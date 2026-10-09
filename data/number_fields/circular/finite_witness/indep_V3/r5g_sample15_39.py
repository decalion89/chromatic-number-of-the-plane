"""Exploration: the same growth (3,1) + deletion over Q(sqrt15) and Q(sqrt39), denominator 80, from 0, the unit
vectors and two points closing a random 5-cycle through 0.  Reports the sizes of the witnesses found."""
import json, random, subprocess, os, re, math, collections
import networkx as nx
os.chdir('sample15')
for f in ['grow.py', 'minimise.py']:
    subprocess.run(['cp', f'../repo_copy/{f}', f], check=True)
def units(d, D):
    out = set()
    B = math.isqrt(D * D // d)
    for b in range(-B, B + 1):
        for e in range(-B, B + 1):
            r = D * D - d * (b * b + e * e)
            if r < 0: continue
            for a in range(-D, D + 1):
                c2 = r - a * a
                if c2 < 0: continue
                c = math.isqrt(c2)
                if c * c != c2: continue
                for cc in {c, -c}:
                    if a * b + cc * e == 0: out.add((a, b, cc, e))
    return sorted(out)
add = lambda p, q: tuple(x + y for x, y in zip(p, q)); zero = (0, 0, 0, 0)
for d, D in [(15, 80), (39, 80)]:
    U = units(d, D); Us = set(U)
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
    print(f'Q(sqrt{d}), D = {D}: {len(U)} unit vectors, {len(cycles)} directed 5-cycles through 0 with two extra points', flush=True)
    if not cycles:
        continue
    rng = random.Random(d); stats = collections.Counter()
    for k in rng.sample(range(len(cycles)), min(12, len(cycles))):
        a, b, c, dd = cycles[k]
        json.dump({'d': d, 'D': D, 'points': [list(p) for p in sorted({zero} | Us | {b, c})]}, open('seed.json', 'w'))
        out = subprocess.run(['python3', '-B', 'grow.py', 'seed.json', 'A', 'G', '3000', '200', '3', '1'], capture_output=True, text=True).stdout
        m = re.search(r'UNSAT with (\d+) points', out)
        if not m:
            print('  ', k, 'no UNSAT:', out.strip().splitlines()[-1][:120], flush=True); continue
        subprocess.run(['python3', '-B', 'minimise.py', 'G.json', 'W.json'], capture_output=True, text=True)
        W = json.load(open('W.json'))
        print(f'   seed {k}: grown {m.group(1)} -> witness {len(W["points"])} vertices, {len(W["edges"])} edges', flush=True)
        stats[len(W['points'])] += 1
        if len(W['points']) == 9:
            subprocess.run(['cp', 'W.json', f'W9_q{d}_{k}.json'])
    print('   sizes:', dict(stats), flush=True)
