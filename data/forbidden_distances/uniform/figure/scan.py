"""Floating-point values of the two-distance linear program, for the figure only (nothing here is a proof).
1/LP(d) for d in [1.005, 6.6] (step 0.005) -> lowscan.tsv, and on a logarithmic grid of [6.6, 10^4] -> logscan.tsv
(columns: d, 1/LP, and for lowscan the coefficients z0, z1, z2). Then lpcurve.dat, the curve drawn in uniform_fig.tex.
usage: python3 scan.py   (about an hour; most of it for d > 1000)"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lpcp import lp_cp
with open('lowscan.tsv', 'w') as out:
    for d in np.arange(1.005, 6.6001, 0.005):
        v, z = lp_cp([d]); out.write(f"{d:.3f}\t{1/v:.5f}\t{z[0]:.6f}\t{z[1]:.6f}\t{z[2]:.6f}\n")
with open('logscan.tsv', 'w') as out:
    for x in list(np.arange(np.log10(6.6), 3.0, 0.01)) + list(np.arange(3.0, 4.0001, 0.05)):
        v, z = lp_cp([10 ** x]); out.write(f"{10 ** x:.6g}\t{1/v:.5f}\n")
lo = np.loadtxt('lowscan.tsv')[:, :2]; hi = np.loadtxt('logscan.tsv')
mu = 0.4027593576                                   # -min J0; at d = 1 the program has one distance: (1 + mu)/mu
D = np.vstack([[[1.0, (1 + mu) / mu]], lo, hi]); D = D[np.argsort(D[:, 0])]
keep = [i for i in range(len(D)) if D[i, 0] <= 2.2 or D[i, 0] > 6.6 or i % 3 == 0]
with open('lpcurve.dat', 'w') as f:
    f.write('d inv\n')
    for i in keep: f.write(f'{D[i,0]:.6g} {D[i,1]:.5f}\n')
