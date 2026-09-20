"""Ask the target question directly: is rho(G,5) at most 63?

Computing rho exactly is expensive in one direction only -- proving a set
forcing is UNSAT -- so asking for the VALUE pays that price repeatedly.  The
question that actually matters is a decision, and it has a cheap side.

rho <= 63 says some 63 vertices meet every realisable colour class.  Keep a
finite family F of classes and alternate:

  (a) is there a set of at most 63 vertices hitting every class in F?
      A cardinality-constrained SAT call.  If NO, then no 63-set hits even
      this subfamily, so none hits all of them: rho > 63, DECIDED, cheaply.

  (b) given such a set S, is there a proper 5-colouring whose colour-0 class
      misses S entirely?  If YES, that class joins F and the loop turns.
      If NO, S is forcing and rho <= 63, decided the expensive way.

Every round of (a) that succeeds and (b) that answers costs seconds. The loop
ends at a real verdict either way, and the cheap end is the one that would
rule de Grey's G out for a core of three.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver


def decide(graph, k, budget, rounds=100000, report=20):
    n = graph.n

    def x(v, c):
        return 1 + v * k + c

    colcls = [[x(v, c) for c in range(k)] for v in range(n)]
    for a, b in graph.edges():
        for c in range(k):
            colcls.append([-x(a, c), -x(b, c)])
    colour = Solver(name="cd19", bootstrap_with=colcls)

    fam, t0 = [], time.time()
    # hitting-set side: variable v+1 means "vertex v is in S"
    pool = IDPool(start_from=n + 1)
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=budget,
                          vpool=pool, encoding=EncType.seqcounter)
    hit = Solver(name="cd19", bootstrap_with=card.clauses)
    try:
        for rnd in range(rounds):
            if not hit.solve():
                print(f"  round {rnd}: no {budget}-set hits the {len(fam)} "
                      f"classes collected -- rho > {budget}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
                return False, len(fam)
            m = hit.get_model()
            S = [v for v in range(n) if m[v] > 0]
            if rnd % report == 0:
                print(f"  round {rnd}: |F| = {len(fam)}, "
                      f"candidate set of {len(S)}  [{time.time()-t0:.0f}s]",
                      flush=True)
            if not colour.solve(assumptions=[-x(v, 0) for v in S]):
                print(f"  round {rnd}: no colouring escapes it -- S of "
                      f"{len(S)} is FORCING, rho <= {budget}  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
                return True, S
            mm = set(colour.get_model())
            cls = [v for v in range(n) if x(v, 0) in mm]
            fam.append(cls)
            hit.add_clause([v + 1 for v in cls])
    finally:
        colour.delete()
        hit.delete()
    return None, len(fam)


if __name__ == "__main__":
    from hn import degrey
    from hn.graph import build_graph

    print("control -- Sa at four colours, where rho = 7:", flush=True)
    sa = build_graph(degrey.build_Sa())
    print("  budget 7:", decide(sa, 4, 7, report=50)[0], flush=True)
    print("  budget 6:", decide(sa, 4, 6, report=50)[0], flush=True)

    print("\nde Grey's G at five colours, budget 63:", flush=True)
    ok, info = decide(degrey.build_G(), 5, 63, report=20)
    print(f"  verdict: rho <= 63 is {ok}")
