"""Referee check 10: compare the reproduced W7.json / C7.json with the stored witness_q7.json.gz."""
import gzip, json, sys

stored = json.load(gzip.open(sys.argv[1], 'rt'))
W = json.load(open(sys.argv[2])); C = json.load(open(sys.argv[3]))
D = stored['denominator']
print('reproduced keys:', sorted(W), '; certificate file keys:', sorted(C) if isinstance(C, dict) else type(C))
print('d, D, p, q equal:', (W['d'], W['denominator'], W.get('p'), W.get('q')) == (stored['d'], D, stored['p'], stored['q']))
for t in [(-D, 0, 0, 0), (D, 0, 0, 0)]:
    # found point f corresponds to stored point s with f = s + t
    S = [tuple(p) for p in stored['points']]
    F = [tuple(p) for p in W['points']]
    m = {}
    for i, f in enumerate(F):
        s = tuple(a - b for a, b in zip(f, t))
        if s in S:
            m[i] = S.index(s)
    print(f'translation found = stored + {tuple(x / D for x in t[:1])}: matched {len(m)} of 9')
    if len(m) == 9:
        break
print('map found -> stored:', m)
Ef = {tuple(sorted((m[a], m[b]))) for a, b in W['edges']}
Es = {tuple(sorted(e)) for e in stored['edges']}
print('edges correspond:', Ef == Es)
print('colouring corresponds:', all(W['colouring'][i] == stored['colouring'][m[i]] for i in range(9)))


def canon_cycle(c):
    k = c.index(min(c)); return tuple(c[k:] + c[:k])


Cf = sorted(canon_cycle([m[v] for v in c]) for c in W['cycles'])
Cs = sorted(canon_cycle(list(c)) for c in stored['cycles'])
print('cycles correspond (as sets of rotated directed cycles):', Cf == Cs)
print('cycles correspond in the same order:', [canon_cycle([m[v] for v in c]) for c in W['cycles']] ==
      [canon_cycle(list(c)) for c in stored['cycles']])
cert = C['critical_colourings'] if isinstance(C, dict) and 'critical_colourings' in C else C
if isinstance(C, dict):
    for k in C:
        if k not in ('critical_colourings',):
            print('   certificate file field', k, ':', str(C[k])[:200])
inv = {v: k for k, v in m.items()}
ok = True
for i in range(9):          # certificate for found vertex i <-> stored certificate for vertex m[i]
    cf = cert[i]; cs = stored['critical_colourings'][m[i]]
    ok &= all(cf[j] == cs[m[j]] for j in range(9))
print('criticality certificates correspond:', ok)
