import random, itertools
from fractions import Fraction as Fr
from kappa_exact import Group, symmetrize, kappa, tight_set, positive_relation
from test_kappa_relations import generates
G=Group(2,1)
rng=random.Random(7)
found={}
cands=[]
for a in itertools.product(range(-3,4),repeat=6):
    pass
tries=0
while tries<400 and len(found)<12:
    tries+=1
    k=rng.choice([3,4])
    Sh=[(rng.randint(-3,3),rng.randint(-3,3),0) for _ in range(k)]
    S=symmetrize(G,Sh)
    if any(s[0]==0 and s[1]==0 for s in S): continue
    if not generates(G,S): continue
    kap,opt=kappa(G,S)
    if Fr(1,4)<kap<Fr(1,3):
        chi=1/kap
        if chi not in found:
            found[chi]=Sh
            T=tight_set(G,S,opt[0][0],opt[0][1],kap)
            print(chi, Sh, 'opt', [tuple(map(str,o[0])) for o in opt][:4], 'T',T, 'rel',positive_relation(G,T), flush=True)
