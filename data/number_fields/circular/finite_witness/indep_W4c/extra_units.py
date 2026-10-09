#!/usr/bin/env python3
"""Informational: list the unit-distance non-edge differences (from the
geometry checker's JSON output) and try to write each, exactly, as a word
in the rotations of the construction with larger exponents:
  Q(sqrt3,sqrt11): R60^j RA^k RG^l (1,0),  |j|,|k|,|l| <= 4
  Q(sqrt2,sqrt3):  zeta^j w^l,              j in Z/24, |l| <= 6
Every identification found is verified in exact arithmetic.

Usage: python3 extra_units.py WITNESS.json.gz GEOM_RESULT.json
"""
import gzip
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from biquad import Field, vec_to_int8  # noqa: E402
from check_geom import elt, cpow  # noqa: E402
from fractions import Fraction

with gzip.open(sys.argv[1], "rt") as fh:
    d = json.load(fh)
res = json.load(open(sys.argv[2]))
extra = res["stats"]["extra_unit_differences_up_to_sign"]
if isinstance(extra, str):
    extra = eval(extra)  # stored via default=str only if not JSON-native
D = d["denominator"]
if d["field"] == "Q(sqrt3, sqrt11)":
    F = Field(3, 11)
    R60 = (elt(Fraction(1, 2)), elt(0, Fraction(1, 2)))
    RA = (elt(Fraction(5, 6)), elt(0, 0, Fraction(1, 6)))
    RG = (elt(Fraction(11, 14)), elt(0, Fraction(5, 14)))
    e1 = (elt(1), elt(0))
    table = {}
    for j in range(-4, 5):
        for k in range(-4, 5):
            for l in range(-4, 5):
                v = F.cmul(F.cmul(cpow(F, R60, j), cpow(F, RA, k)), F.cmul(cpow(F, RG, l), e1))
                try:
                    t = vec_to_int8(v, D)
                except ValueError:
                    continue
                table.setdefault(t, (j, k, l))
    label = "R60^j RA^k RG^l"
else:
    F = Field(2, 3)
    zeta = (elt(0, Fraction(1, 4), 0, Fraction(1, 4)), elt(0, Fraction(-1, 4), 0, Fraction(1, 4)))
    w = (elt(Fraction(1, 3)), elt(0, Fraction(2, 3)))
    table = {}
    for j in range(24):
        for l in range(-6, 7):
            v = F.cmul(cpow(F, zeta, j), cpow(F, w, l))
            try:
                t = vec_to_int8(v, D)
            except ValueError:
                continue
            table.setdefault(t, (j, l))
    label = "zeta^j w^l"
sa, sb = math.sqrt(F.A), math.sqrt(F.B)
fl = lambda t: (t[0] + t[1] * sa + t[2] * sb + t[3] * sa * sb) / D
print("unit-distance non-edge differences (up to sign), with multiplicity:")
tot = 0
for vec, cnt in extra:
    t = tuple(vec)
    nt = tuple(-x for x in t)
    ang = math.degrees(math.atan2(fl(t[4:]), fl(t[:4])))
    word = table.get(t)
    sgn = "+"
    if word is None and nt in table:
        word, sgn = table[nt], "-"
    tot += cnt
    if word:
        desc = "= %s%s with exponents %s" % (sgn, label, word)
    else:
        desc = "not of the form +-%s in the searched exponent range" % label
    print("  %-38s x%-3d angle %9.4f deg   %s" % (list(t), cnt, ang, desc))
print("total unit non-edge pairs:", tot)
