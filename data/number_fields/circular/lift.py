"""Follow characters of (1/N)Z[i] coherently: a character at N = 5^k (point c mod N Z[i]) lifts to N' = 5N as one of
the 25 points c + N*lam (lam mod 5).  Keep the lifts whose kappa_{k+1} (min over G_{N'} of ||Re(conj(c) gamma)||,
locally maximised) stays above THR; report the survivors and their best kappa at each level."""
import re, sys, numpy as np
THR = float(sys.argv[1]); KMAX = int(sys.argv[2])
cents = []
for line in open(sys.argv[3] if len(sys.argv) > 3 else "types_k3.txt"):
    m = re.search(r"kappa ([0-9.]+)\s+c = \(([0-9.]+), ([0-9.]+)\)", line)
    if m and float(m.group(1)) >= THR: cents.append(complex(float(m.group(2)), float(m.group(3))))
rho = complex(3, 4) / 5
def G(k):
    return np.array([e * rho ** j for j in range(-k, k + 1) for e in (1, 1j, -1, -1j)])
def kap(C, Gk):
    kv = np.full(C.shape, 0.5)
    for g in Gk:
        v = (np.conj(C) * g).real
        np.minimum(kv, np.abs(v - np.round(v)), out=kv)
    return kv
s = np.arange(-0.15, 0.1501, 0.005)
X, Y = np.meshgrid(s, s, indexing="ij"); OFF = X + 1j * Y
def refine(c, Gk):
    best = (-1, c)
    for _ in range(3):
        K = kap(best[1] + OFF, Gk)
        a = np.unravel_index(np.argmax(K), K.shape)
        cand = (float(K[a]), best[1] + OFF[a])
        if cand[0] <= best[0] + 1e-12: break
        best = cand
    return best
live = [(c, 3) for c in cents]
print(f"k = 3: {len(live)} starting points with kappa >= {THR}")
for k in range(4, KMAX + 1):
    N0 = 5 ** (k - 1); Gk = G(k); nxt = []
    for c, _ in live:
        for a in range(5):
            for b in range(5):
                kv, cc = refine(c + N0 * complex(a, b), Gk)
                if kv >= THR:
                    nxt.append((cc, kv))
    # deduplicate modulo N Z[i]
    N = 5 ** k; uniq = []
    for cc, kv in sorted(nxt, key=lambda t: -t[1]):
        z = complex(cc.real % N, cc.imag % N)
        if all(abs(((z - u).real + N / 2) % N - N / 2) + abs(((z - u).imag + N / 2) % N - N / 2) > 0.5 for u, _ in uniq):
            uniq.append((z, kv))
    live = [(z, k) for z, _ in uniq]
    print(f"k = {k}: {len(uniq)} survivors; kappa:", " ".join(f"{kv:.4f}" for _, kv in uniq[:16]))
    for z, kv in uniq[:16]:
        print(f"    c/N = ({z.real / N:.5f}, {z.imag / N:.5f})  kappa {kv:.4f}")
    if not live: break
