"""Referee E: method A (ref_clip.run) with several functional orders; the sets of index vectors
(as maps functional -> index) must coincide."""
import random
from fractions import Fraction as Fr
from ref_clip import run
res = {}
for th in (Fr(7, 25), Fr(201, 700)):
    sets = []
    for seed in (None, 1, 2, 3):
        order = list(range(16))
        if seed is not None:
            random.Random(seed).shuffle(order)
        funcs, ford, leaves, nodes = run(th, order=order)
        sets.append(frozenset(frozenset(zip(ford, idx)) for _, idx in leaves))
        print("theta=%s order seed %s: %d leaves, %d nodes" % (th, seed, len(leaves), nodes))
    print("  all orders give the same set of index vectors:", all(s == sets[0] for s in sets))
