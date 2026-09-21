"""Not "at least one" but "at least how many"?

The witness says: in every proper 5-colouring, at least one of its 161 pairs
is monochromatic.  That is the weakest form of the statement.  The real
quantity is the MINIMUM, over colourings, of how many of the 161 are
monochromatic at once -- and it is computable exactly, not estimated.

A pair is monochromatic exactly when its two endpoints agree, so introduce one
indicator per pair, tie it to the five same-colour cases, and ask for a
colouring with at most t indicators true.  Binary search on t and the smallest
feasible value is the answer.

A larger minimum is a stronger hypothesis and a stronger hypothesis is easier
to consume: "one of 161" gives a rotated copy 160 ways to escape, while "ten
of 161" forces ten simultaneous coincidences and leaves far less room.  If the
minimum turns out to be one, the hypothesis really is as weak as it looks.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.field import Field
from hn.geometry import Point
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = 5
t0 = time.time()
doc = json.load(open("/home/user/darwin-50/research/hadwiger-nelson/"
                     "data/witness_five.json"))
K = Field(tuple(doc["field"]["generators"]))
P = [Point(K.element([Fr(c) for c in x]), K.element([Fr(c) for c in y]))
     for x, y in doc["points"]]
n = len(P)
E = [tuple(e) for e in doc["unit_edges"]]
FORB = [(a, b) for a, b, _ in doc["forbidden_pairs"]]
print(f"witness: {n} points, {len(E)} unit edges, {len(FORB)} pairs"
      f"  [{time.time()-t0:.0f}s]", flush=True)

pool = IDPool(start_from=n * k + 1)
cnf = CNF()
for v in range(n):
    cnf.append([1 + v * k + c for c in range(k)])
for a, b in E:
    for c in range(k):
        cnf.append([-(1 + a * k + c), -(1 + b * k + c)])
# m[p] is true when pair p is monochromatic
ind = []
for p, (a, b) in enumerate(FORB):
    mv = pool.id(("mono", p))
    ind.append(mv)
    # same colour in colour c  =>  m
    for c in range(k):
        cnf.append([-(1 + a * k + c), -(1 + b * k + c), mv])
    # m  =>  they agree somewhere (so m cannot be set for free)
    cnf.append([-mv] + [pool.id(("ag", p, c)) for c in range(k)])
    for c in range(k):
        ag = pool.id(("ag", p, c))
        cnf.append([-ag, 1 + a * k + c])
        cnf.append([-ag, 1 + b * k + c])

for t in range(0, 25):
    card = CardEnc.atmost(lits=ind, bound=t, vpool=pool,
                          encoding=EncType.seqcounter)
    sv = Solver(name="cd15", bootstrap_with=cnf.clauses + card.clauses)
    ok = sv.solve()
    sv.delete()
    print(f"   at most {t:2d} monochromatic -> "
          f"{'possible' if ok else 'impossible'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if ok:
        print(f"\n>>> the minimum is {t}: every proper 5-colouring of the "
              f"witness leaves at least {t} of its {len(FORB)} pairs "
              f"monochromatic", flush=True)
        break
print("DONE", flush=True)
