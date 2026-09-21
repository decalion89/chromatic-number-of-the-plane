"""Which shapes of hub disjunction close, and which are free?

The lemma's strength is entirely decided by one finite question.  A disjunction
over a set W of ring points, written as exponents of the spindle rotation, is
defeated exactly when the shifted copies W, W+1, ..., W+T admit an independent
transversal in the integers -- a choice of one exponent from each, no two
consecutive.  A transversal is a colouring, so shapes that admit one are free
and the lemma says nothing about them.

Width two settles into "even closes at m+1 copies, odd never".  But a wider W
is a weaker hypothesis and therefore an easier thing to find in a graph, so the
question worth answering is how far the width can go before every shape becomes
free.  Enumerated exhaustively over every shape inside a window.
"""
import sys
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from itertools import combinations
from hn.homcol import independent_transversal_exists as ite

CAP = 26           # copies to try before calling a shape free
for L in (12,):
    print(f"window 0..{L}, copies up to {CAP}\n", flush=True)
    for size in (1, 2, 3, 4):
        killers = []
        total = 0
        for rest in combinations(range(1, L + 1), size - 1):
            W = (0,) + rest
            total += 1
            need = next((c for c in range(1, CAP + 1) if not ite(W, c)), None)
            if need:
                killers.append((W, need))
        print(f"|W|={size}: {len(killers)} of {total} shapes close", flush=True)
        for W, need in killers[:14]:
            print(f"    {W} closes at {need} copies", flush=True)
        if len(killers) > 14:
            print(f"    ... and {len(killers)-14} more", flush=True)
        print(flush=True)
