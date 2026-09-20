"""The decision loop, fixed: shrink the hitting set before asking for an escape.

The first version asked a cardinality-bounded SAT call for ANY 63-set hitting
the family, then asked for a colouring whose colour-0 class escapes it.  That
loop does not converge, and the diagnosis is concrete: over 400 rounds the
escaping classes were large (282 to 314 vertices) but 152 vertices were common
to EVERY ONE of them.  A single vertex hit the whole family, the minimum
hitting set read 2 -- below the floor rho >= 5 -- and the loop was spinning on
near-identical colourings rather than sitting on a boundary.

The cause is slack.  An arbitrary 63-set hits the family through whichever
vertices it likes, so the next class is barely constrained and comes back
almost unchanged.  Shrinking S to a MINIMAL hitting set first removes the
slack: every vertex left in S is needed by some class, so a colouring escaping
S has to avoid all of them, and the class it returns is genuinely new.

Same two questions, same cheap side -- no B-set hits F means rho > B -- but
now each round costs one shrink and buys a class that moves the family.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver
from hn import degrey


def shrink_hitting(S, fam):
    """Drop any vertex of S that no class needs, cheapest first."""
    cover = {v: 0 for v in S}
    sset = set(S)
    for cl in fam:
        hits = [v for v in cl if v in sset]
        if len(hits) == 1:
            cover[hits[0]] += 1
    keep = [v for v in S if cover[v] > 0]
    # re-check: a class may now be uncovered, so add back greedily
    kset = set(keep)
    missing = [cl for cl in fam if not (set(cl) & kset)]
    for cl in missing:
        pick = next((v for v in cl if v in sset), cl[0])
        keep.append(pick)
        kset.add(pick)
    return sorted(set(keep))


def run(g, k, budget, rounds=500000, report=200):
    n = g.n

    def x(v, c):
        return 1 + v * k + c

    cls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in g.edges():
        for c in range(k):
            cls.append([-x(a, c), -x(b, c)])
    colour = Solver(name="cd19", bootstrap_with=cls)
    pool = IDPool(start_from=n + 1)
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=budget,
                          vpool=pool, encoding=EncType.seqcounter)
    hit = Solver(name="cd19", bootstrap_with=card.clauses)
    fam, t0 = [], time.time()
    try:
        for rnd in range(rounds):
            if not hit.solve():
                print(f"  round {rnd}: no {budget}-set hits {len(fam)} "
                      f"classes -- rho > {budget}  [{time.time()-t0:.0f}s]",
                      flush=True)
                return False
            m = hit.get_model()
            S = shrink_hitting([v for v in range(n) if m[v] > 0], fam)
            if not colour.solve(assumptions=[-x(v, 0) for v in S]):
                print(f"  round {rnd}: S of {len(S)} is FORCING -- "
                      f"rho <= {len(S)}  [{time.time()-t0:.0f}s]", flush=True)
                return True
            mm = set(colour.get_model())
            c0 = [v for v in range(n) if x(v, 0) in mm]
            fam.append(c0)
            hit.add_clause([v + 1 for v in c0])
            if rnd and rnd % report == 0:
                inter = set(fam[0])
                for cl in fam[1:]:
                    inter &= set(cl)
                print(f"  round {rnd}: |F| = {len(fam)}, |S| = {len(S)}, "
                      f"class {len(c0)}, common to all {len(inter)}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    finally:
        colour.delete()
        hit.delete()
    return None


if __name__ == "__main__":
    from hn.graph import build_graph
    print("control -- Sa at four colours, budget 5 (floor is 4):", flush=True)
    print("  verdict:", run(build_graph(degrey.build_Sa()), 4, 5, report=500),
          flush=True)
    print("\nde Grey's G at five colours, budget 63:", flush=True)
    print("  verdict:", run(degrey.build_G(), 5, 63, report=200))
