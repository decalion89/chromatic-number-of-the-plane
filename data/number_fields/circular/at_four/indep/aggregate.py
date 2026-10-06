"""Aggregate the C checker's per-S output: how often each claim was exercised."""
import re, glob, os
from fractions import Fraction as Fr
agg = dict(groups=0, S4=0, k_lt=0, k_eq=0, k_gt=0, col=0, ABneq=0, t4_not_sq=0, noT4_not_acyc=0,
           sq_free_at_eq=0, S_eq_with_sqfree=0, S_eq_with_noT4=0, abneq_charfail=0, notcol=0)
for f in sorted(glob.glob('out_c/*.txt')):
    agg['groups'] += 1
    for line in open(f):
        if not line.startswith('S='):
            continue
        kv = dict(t.split('=', 1) for t in line.split()[1:] if '=' in t and not t.startswith('hom'))
        k = Fr(kv['kappa'])
        if 'not4col' in line:
            agg['notcol'] += 1
            continue
        agg['S4'] += 1
        col, acyc, noT4, noSq, ABeq = (int(kv[x]) for x in ('col', 'acyc', 'noT4', 'noSq', 'ABeq'))
        agg['col'] += col
        agg['ABneq'] += col - ABeq
        agg['t4_not_sq'] += int(kv['T4butABeq'])
        agg['noT4_not_acyc'] += noT4 - acyc
        agg['abneq_charfail'] += int(kv['ABneq_charfail'])
        if k < Fr(1, 4): agg['k_lt'] += 1
        elif k == Fr(1, 4):
            agg['k_eq'] += 1
            agg['sq_free_at_eq'] += noSq
            agg['S_eq_with_sqfree'] += noSq > 0
            agg['S_eq_with_noT4'] += noT4 > 0
        else: agg['k_gt'] += 1
for k, v in agg.items():
    print(f'{k:>18}: {v}')
