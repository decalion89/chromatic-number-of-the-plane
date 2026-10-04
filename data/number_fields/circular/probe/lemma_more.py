from lemma_lp import *
import sys
for k in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    N = 5 ** k
    sC = min(C_threshold(k, sg)[0] for sg in (1, -1))
    eQ = None
    for a in (1, 2):
        for b in (1, 2):
            q = (F(N * a, 3), F(N * b, 3))
            for sg in (1, -1):
                o = Q_threshold(k, sg, q)
                eQ = o[0] if eQ is None else min(eQ, o[0])
    print(f"k={k}: r1_C = {HALF - sC} = {float(HALF - sC):.7f}; r1_Q = {F(1,3) - eQ} = {float(F(1,3) - eQ):.7f}")
    sys.stdout.flush()
