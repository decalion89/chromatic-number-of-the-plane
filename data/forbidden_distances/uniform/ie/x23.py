"""The 23 points X23 of Ambrus, Csiszárik, Matolcsi, Varga and Zsámboki, The density of planar sets avoiding unit
distances (arXiv:2207.14179; Math. Program. 207 (2024)), transcribed from their verification notebook
verification23.nb, as exact coordinates a + b sqrt3 + c sqrt11 + d sqrt33 (rationals a, b, c, d) of the real and
imaginary parts, and as floats. The unit-distance graph has 47 edges; the 253 distances take 27 values."""
from fractions import Fraction as F
import numpy as np
# (re: [rational, sqrt3, sqrt11, sqrt33], im: [rational, sqrt3, sqrt11, sqrt33])
P = [
 ([0, 0, 0, 0], [0, 0, 0, 0]),
 ([1, 0, 0, 0], [0, 0, 0, 0]),
 ([F(1, 2), 0, 0, 0], [0, F(1, 2), 0, 0]),
 ([F(3, 2), 0, 0, 0], [0, F(1, 2), 0, 0]),
 ([F(5, 6), 0, 0, 0], [0, 0, F(1, 6), 0]),
 ([F(5, 12), 0, 0, F(-1, 12)], [0, F(5, 12), F(1, 12), 0]),
 ([F(5, 4), 0, 0, F(-1, 12)], [0, F(5, 12), F(1, 4), 0]),
 ([F(13, 12), 0, 0, F(1, 12)], [0, F(1, 12), F(-1, 12), 0]),
 ([F(1, 4), 0, 0, F(1, 12)], [0, F(1, 12), F(-1, 4), 0]),
 ([F(13, 12), 0, 0, F(-1, 12)], [0, F(-1, 12), F(-1, 12), 0]),
 ([F(2, 3), 0, 0, 0], [0, F(1, 2), F(-1, 6), 0]),
 ([F(7, 12), 0, 0, F(-1, 12)], [0, F(5, 12), F(-1, 12), 0]),
 ([F(5, 6), 0, 0, F(-1, 6)], [0, F(-1, 6), F(1, 6), 0]),
 ([F(7, 12), 0, 0, F(1, 12)], [0, F(7, 12), F(-1, 12), 0]),
 ([F(1, 12), 0, 0, F(-1, 12)], [0, F(-1, 12), F(-1, 12), 0]),
 ([F(2, 3), 0, 0, F(-1, 6)], [0, F(1, 3), F(-1, 6), 0]),
 ([F(19, 12), 0, 0, F(-1, 12)], [0, F(5, 12), F(-1, 12), 0]),
 ([F(1, 3), 0, 0, 0], [0, F(1, 2), F(1, 6), 0]),
 ([F(3, 4), 0, 0, F(-1, 12)], [0, F(-1, 12), F(1, 4), 0]),
 ([1, 0, 0, F(-1, 6)], [0, F(-1, 6), 0, 0]),
 ([F(7, 6), 0, 0, 0], [0, 0, F(-1, 6), 0]),
 ([F(4, 3), 0, 0, F(-1, 6)], [0, F(1, 3), F(1, 6), 0]),
 ([F(19, 12), 0, 0, F(-1, 4)], [0, F(1, 4), F(-1, 12), 0]),
]
B = np.array([1.0, np.sqrt(3), np.sqrt(11), np.sqrt(33)])
def fl(v): return float(np.dot([float(c) for c in v], B))
X = np.array([complex(fl(re), fl(im)) for re, im in P])
