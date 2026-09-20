"""Which rotations of K meet Sa at all?

Measured on de Grey's own: Sa and Sb share exactly ONE point, the origin, and
Y carries exactly SIX edges joining the two exclusive halves.  Six edges and a
shared vertex are the whole of what turns a graph with no forced pair into one
that has one -- Sa alone forces nothing, Y forces (2,0),(-2,0).

That is a cheap filter and a sharp one.  A rotation about the origin that
produces no cross edge glues two lumps at a point and can force nothing new;
only the ones that bite are worth a pair scan.  So sweep every unit of K in
hand -- 3030 of them, each a rotation -- count the shared points and the cross
edges by float, and keep the ones that bite.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "kneck2.py").read()
exec(src[:src.index("pairsum = {}")])
from hn.degrey import S_POINTS

SQ3 = ksub(kmul(krat(2), Z6), K1)


def to_k(xs, ys):
    p, q = Fr(xs.get(1, 0)), Fr(xs.get(33, 0))
    r, t = Fr(ys.get(3, 0)), Fr(ys.get(11, 0))
    return kadd(kadd(krat(p), kmul(krat(r), SQ3)),
                kmul(ksub(krat(t), kmul(krat(q), SQ3)), S))


Sa, seenp = [], set()
for z in [to_k(x, y) for x, y in S_POINTS]:
    for base in (z, kconj(z)):
        q = base
        for _ in range(6):
            if q not in seenp:
                seenp.add(q)
                Sa.append(q)
            q = kmul(Z6, q)
zsA = [zof(p) for p in Sa]
saset = set(Sa)
print(f"Sa over K: {len(Sa)} points  [{time.time()-t0:.0f}s]", flush=True)

cellA = {}
for i, z in enumerate(zsA):
    cellA.setdefault((int(z.real // 1), int(z.imag // 1)), []).append(i)

hits = []
for ri, u in enumerate(steps):
    img = [kmul(u, q) for q in Sa]
    imset = set(img)
    shared = len(imset & saset)
    zi = [zof(q) for q in img]
    cross = 0
    for j, z in enumerate(zi):
        if img[j] in saset:
            continue
        cx, cy = int(z.real // 1), int(z.imag // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cellA.get((cx + da, cy + db), ()):
                    if abs(abs(z - zsA[i]) - 1) < 1e-9:
                        cross += 1
    if cross:
        hits.append((cross, shared, ri))
        print(f"  rotation {ri}: {shared} shared, {cross} cross edges  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    if ri % 300 == 0:
        print(f"  ... {ri}/{len(steps)}, {len(hits)} biting so far  "
              f"[{time.time()-t0:.0f}s]", flush=True)
hits.sort(reverse=True)
print(f"\n{len(hits)} rotations of K bite Sa; best "
      f"{hits[:10]}  [{time.time()-t0:.0f}s]", flush=True)
import pickle
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/kcross.pkl", "wb") as fh:
    pickle.dump([(c, s, flat(steps[r])) for c, s, r in hits], fh)
