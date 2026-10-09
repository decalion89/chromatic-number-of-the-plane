from fractions import Fraction as Fr
from kappa_exact import Group, symmetrize, kappa, half, gval
import itertools, random
G=Group(1,1)
def brute1(D, den=2000):
    best=0
    for k in range(den):
        a=k/den
        best=max(best,min(min((a*d)%1,1-(a*d)%1) for d in D))
    return best
for D in [(1,4),(2,3),(1,2),(1,2,3),(1,2,5),(3,4),(2,5),(1,6),(2,3,7),(1,5,8),(3,5,7)]:
    S=symmetrize(G,[(d,0) for d in D])
    k,opt=kappa(G,S)
    a,b=D[0],D[1]
    pred = Fr(a+b-1,2*(a+b)) if len(D)==2 and (a+b)%2==1 else None
    print(D,k,float(k),round(brute1(D),4),'pred2:',pred,'opt:',[(str(o[0][0]),o[1]) for o in opt][:6])
# Z^2 check against brute force grid
G2=Group(2,1)
rng=random.Random(5)
for trial in range(6):
    Sh=[(rng.randint(-2,2),rng.randint(-2,2),0) for _ in range(3)]
    if any(s[0]==0 and s[1]==0 for s in Sh): continue
    S=symmetrize(G2,Sh)
    try:
        k,opt=kappa(G2,S)
    except Exception as e:
        print('err',e); continue
    den=240
    best=0
    for i in range(den):
        for j in range(den):
            a=(i/den,j/den)
            v=min(min((a[0]*s[0]+a[1]*s[1])%1, 1-(a[0]*s[0]+a[1]*s[1])%1) for s in Sh)
            best=max(best,v)
    print(Sh,k,float(k),round(best,4))
