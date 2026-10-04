"""Every triangle-free graph with at most 8 vertices has a homomorphism to K_{8/3}, so its circular chromatic number
is at most 8/3 < 3; hence a triangle-free graph with chi_c = 3 (such as witness_q7, the Wagner graph with one chord
subdivided) has at least 9 vertices.

usage: python3 small_triangle_free.py [N]   (N = 8 by default; N <= 8, since K_{8/3} is only the right target there)

It is enough to check the maximal triangle-free graphs on exactly N vertices: pad a smaller graph with isolated
vertices, then add edges while no triangle appears; a homomorphism of the larger graph restricts to the smaller one.
All triangle-free graphs on the labelled vertices 0..N-1 are generated, each once: vertex v is joined to an
independent set of {0, ..., v-1}.  For the maximal ones (every non-adjacent pair has a common neighbour) a
homomorphism to K_{8/3} (vertices Z/8, i ~ j iff 3 <= (i - j) mod 8 <= 5) is searched by backtracking.
Output for N = 8: 4682270 triangle-free graphs, 15247 maximal, 15120 of them not bipartite, all mapping to K_{8/3}
(under a minute).  The counts 7, 41, 388, 5789 for N = 3, ..., 6 agree with a direct enumeration of all
2^(N(N-1)/2) graphs on N labelled vertices.
"""
import itertools, sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
assert 1 <= N <= 8
P, Q = 8, 3
OK = [[Q <= (a - b) % P <= P - Q for b in range(P)] for a in range(P)]


def hom(adj):
    """True if the graph (adjacency bitmasks) has a homomorphism to K_{8/3}."""
    order = sorted(range(N), key=lambda v: -bin(adj[v]).count('1'))
    nb = [[w for w in range(N) if adj[v] >> w & 1] for v in range(N)]
    c = [-1] * N

    def go(t):
        if t == N:
            return True
        v = order[t]
        for k in ((0,) if t == 0 else range(P)):
            if all(c[w] < 0 or OK[k][c[w]] for w in nb[v]):
                c[v] = k
                if go(t + 1):
                    return True
                c[v] = -1
        return False
    return go(0)


def bipartite(adj):
    col = [-1] * N
    for s in range(N):
        if col[s] >= 0:
            continue
        col[s] = 0; st = [s]
        while st:
            v = st.pop()
            for w in range(N):
                if adj[v] >> w & 1:
                    if col[w] < 0:
                        col[w] = 1 - col[v]; st.append(w)
                    elif col[w] == col[v]:
                        return False
    return True


count = maximal = nonbip = 0
fails = []
adj = [0] * N


def rec(v):
    global count, maximal, nonbip
    if v == N:
        count += 1
        for i in range(N):
            for j in range(i + 1, N):
                if not (adj[i] >> j & 1) and not (adj[i] & adj[j]):
                    return
        maximal += 1
        if not bipartite(adj):
            nonbip += 1
            if not hom(adj):
                fails.append(list(adj))
        return
    for r in range(v + 1):
        for S in itertools.combinations(range(v), r):
            m = 0
            for s in S:
                if adj[s] & m:
                    break
                m |= 1 << s
            else:
                for s in S:
                    adj[s] |= 1 << v
                adj[v] = m
                rec(v + 1)
                for s in S:
                    adj[s] &= ~(1 << v)
                adj[v] = 0


rec(0)
print(f'triangle-free graphs on {N} labelled vertices: {count}; maximal: {maximal}; not bipartite: {nonbip}; '
      f'without a homomorphism to K_8/3: {len(fails)}')
assert not fails, fails[:3]
