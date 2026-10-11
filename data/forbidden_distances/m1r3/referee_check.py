"""The internal referee's own check of an M1R3 certificate, written separately from both checkers (kept as it ran):
brute force over all 2^n subsets for the independent sets, exact min Q, every congruence (and how many need a
reflection), Phi - m on a dense floating-point grid (not a proof), and the positive root of c x^2 + m x - psiT.
usage: python3 -I referee_check.py CERT.json   (needs 2^n bytes of memory and a few minutes for n = 27)"""
import json, sys, itertools
from fractions import Fraction as F
import numpy as np, mpmath as mp
C=json.load(open(sys.argv[1]))
P=[tuple(F(c) for c in p) for p in C['points']]; n=len(P)
d2=lambda i,j: sum((a-b)**2 for a,b in zip(P[i],P[j]))/2
adj=[[i!=j and d2(i,j)==1 for j in range(n)] for i in range(n)]
# brute force all 2^n subsets with numpy
N=1<<n; ok=np.ones(N,dtype=bool); idx=np.arange(N,dtype=np.int64)
for i in range(n):
  for j in range(i+1,n):
    if adj[i][j]: ok &= ~(((idx>>i)&1).astype(bool) & ((idx>>j)&1).astype(bool))
S=idx[ok]; print('independent sets (brute force):',len(S))
psiT=F(C['psiT']); lam=[F(x) for x in C['lam']]; U={(a-1,b-1):F(v) for a,b,v in C['u']}
WC=[([a-1 for a in I],[b-1 for b in J],F(v)) for I,J,v in C['wc']]
qmin=None
for s in S.tolist():
  mem=[i for i in range(n) if s>>i&1]; ms=set(mem)
  q=psiT+sum(lam[i] for i in mem)+sum(v for (a,b),v in U.items() if a in ms and b in ms)
  q+=sum(v*((all(a in ms for a in I))-(all(b in ms for b in J))) for I,J,v in WC)
  if qmin is None or q<qmin: qmin=q;arg=mem
print('min Q exact', float(qmin), arg)
# congruence: independent check incl. chirality
import numpy.linalg as la
s2=np.sqrt(0.5); X=np.array([[float(c) for c in p] for p in P])*s2
nchiral=0; bad=0
for I,J,v in WC:
  if any(d2(I[a],I[b])!=d2(J[a],J[b]) for a in range(len(I)) for b in range(a+1,len(I))): bad+=1
  A=X[I]-X[I].mean(0); B=X[J]-X[J].mean(0)
  if len(I)>=4 and la.matrix_rank(A,1e-9)==3:
    R=la.lstsq(A,B,rcond=None)[0]
    if la.det(R)<0: nchiral+=1
print('wc bad',bad,'wc entries needing a reflection (non-coplanar, det<0):',nchiral, 'sizes', sorted(set(len(I) for I,_,_ in WC)))
# numeric Phi
psi0=F(C['psi0']); m=F(C['m'])
G={}
for (a,b),v in U.items(): G[d2(a,b)]=G.get(d2(a,b),0)+v
mp.mp.dps=30
def phi(t):
  t=mp.mpf(t); s=lambda x: mp.sinc(x)
  return mp.mpf(psi0.numerator)/psi0.denominator*s(t)-mp.mpf(float(sum(lam)))-sum(mp.mpf(float(v))*s(mp.sqrt(float(q))*t) for q,v in G.items())
ts=np.linspace(1e-6,100,200001)
vals=np.array([float(psi0)*np.sinc(t/np.pi)-float(sum(lam))-sum(float(v)*np.sinc(np.sqrt(float(q))*t/np.pi) for q,v in G.items()) for t in ts])
k=vals.argmin(); print('grid min Phi',vals[k],'at t=',ts[k],' m=',float(m),' margin',vals[k]-float(m))
c=psi0-sum(lam)-sum(U.values())-m
print('Phi(0)=',float(psi0-sum(lam)-sum(U.values())),'c=',float(c))
r=(-mp.mpf(float(m))+mp.sqrt(mp.mpf(float(m))**2+4*mp.mpf(float(c))*mp.mpf(float(psiT))))/(2*mp.mpf(float(c)))
print('root',r)
