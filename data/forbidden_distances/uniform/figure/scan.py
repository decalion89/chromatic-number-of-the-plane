"""Floating-point values for the figure only (nothing here is a proof).
- lowscan.tsv: 1/LP(d), where LP(d) is the value of the linear program for the single distance d with unit-triangle
  cliques, for d in [1.01, 6.6] (step 0.005 on [1.6, 2.2], 0.02 elsewhere; envelope row at T = 200);
- tailscan.tsv: 1/value of the program of Lemma 18 (one function for all d >= D) with unit triangles and Moser
  spindles, for D on a logarithmic grid of [5, 2*10^5] (step 10^0.1) and at 10^6 and 10^7.
Then lpcurve.dat and tailcurve.dat (up to 2*10^5), the two curves drawn in uniform_fig.tex. Columns: d (or D), 1/value.
usage: python3 scan.py [low] [tail]   (both by default; about an hour each)"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
todo = sys.argv[1:] or ['low', 'tail']
if 'low' in todo:
    from lpcp import lp_cp
    with open('lowscan.tsv', 'w') as out:
        ds = sorted(set(np.round(np.concatenate([np.arange(1.01, 1.6, 0.02), np.arange(1.6, 2.2001, 0.005),
                                                 np.arange(2.22, 6.6001, 0.02)]), 4)))
        for d in ds:
            v = lp_cp([d], T=200.0, cliques=True)[0]; out.write(f"{d:.3f}\t{1/v:.5f}\n"); out.flush()
if 'tail' in todo:
    from tail import float_lp
    with open('tailscan.tsv', 'w') as out:
        for x in list(np.arange(np.log10(5), np.log10(2e5) + 1e-9, 0.1)) + [6.0, 7.0]:
            v = float_lp(10 ** x, True, True)[0]; out.write(f"{10 ** x:.6g}\t{1/v:.5f}\n"); out.flush()
for src, dst in (('lowscan.tsv', 'lpcurve.dat'), ('tailscan.tsv', 'tailcurve.dat')):
    if os.path.exists(src):
        D = np.loadtxt(src)[:, :2]; D = D[np.argsort(D[:, 0])]
        with open(dst, 'w') as f:
            f.write('d inv\n')
            for d, v in D[D[:, 0] <= 2e5]: f.write(f'{d:.6g} {v:.5f}\n')
