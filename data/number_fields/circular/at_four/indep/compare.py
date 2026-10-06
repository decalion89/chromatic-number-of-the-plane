"""Compare the per-S statistics of the C checker (canonical S up to automorphisms) with those of the
independent Python checker (all S) on every group run by both."""
import re, glob, os
from fractions import Fraction as Fr
def parse(path):
    out = {}
    for line in open(path):
        if not line.startswith('S='):
            continue
        sp = line.split()
        S = frozenset(tuple(int(x) for x in g.split(',')) for g in re.findall(r'\(([^)]*)\)', sp[0]))
        kv = dict(t.split('=', 1) for t in sp[1:] if '=' in t and not t.startswith('hom'))
        rec = {'kappa': Fr(kv['kappa'])}
        if 'not4col' in line:
            rec['col'] = 0
        else:
            for k in ('col', 'acyc', 'noT4', 'noSq', 'ABeq'):
                rec[k] = int(kv[k])
            rec['hom'] = int(re.search(r'hom\d+/\d+=(\d)', line).group(1))
            rec['ok'] = line.rstrip().endswith('ok')
        out[S] = rec
    return out
tot = mism = 0
for pf in sorted(glob.glob('out_py/*.txt')):
    cf = 'out_c/' + os.path.basename(pf)
    if not os.path.exists(cf):
        continue
    P, C = parse(pf), parse(cf)
    complete = any(l.startswith('# summary') for l in open(pf))
    n_here = 0
    for S, rc in C.items():
        tot += 1; n_here += 1
        rp = P.get(S)
        if not complete and rp is None:
            tot -= 1; n_here -= 1
            continue
        if rc['col'] == 0:
            if rp is not None:
                mism += 1; print('MISMATCH (4-colourability)', pf, sorted(S))
            continue
        if rp is None or any(rp[k] != rc[k] for k in ('kappa', 'col', 'acyc', 'noT4', 'noSq', 'ABeq', 'hom')) or not (rp['ok'] and rc['ok']):
            mism += 1; print('MISMATCH', pf, sorted(S), rc, rp)
    print(f'{os.path.basename(pf)}: {n_here} canonical S compared, python lines {len(P)}' + ('' if complete else ' (python run incomplete)'))
print(f'TOTAL canonical S compared: {tot}, mismatches: {mism}')
