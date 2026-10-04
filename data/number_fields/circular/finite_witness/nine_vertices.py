"""The triangle-free graphs with 9 vertices and circular chromatic number 3.

usage: python3 nine_vertices.py      (needs networkx, for isomorphism tests; about a minute)

Generates all triangle-free graphs on n = 1, ..., 9 vertices up to isomorphism: each graph on n vertices extends one
on n - 1 vertices by a vertex joined to an independent set, and isomorphic copies are removed (Weisfeiler-Lehman hash,
then an exact isomorphism test).  The counts 1, 2, 3, 7, 14, 38, 107, 410, 1897 are those of OEIS A006785.  A
triangle-free graph with at most 10 vertices has chi <= 3 (the smallest triangle-free 4-chromatic graph has 11), so for
these graphs chi_c = 3 exactly when
  (a) there is no homomorphism to K_{8/3} (by backtracking; chi_c is a fraction with numerator at most 9, and 8/3 is
      the largest one below 3), and, independently,
  (b) every proper 3-colouring has a tight cycle (a directed cycle along which the colour increases by 1 mod 3),
      by enumerating the 3-colourings with c(0) = 0 (Lemma 20 of papers/three-colours).
The program checks that (a) and (b) agree on every graph, and reports the graphs with chi_c = 3: on 9 vertices there
are exactly three, the hexagon with its three long diagonals subdivided (12 edges), and that graph with one or two
edges added between the midpoints of the diagonals (13 edges: witness_q7, the Wagner graph with one chord
subdivided; 14 edges: witness_q31).
"""
import itertools, sys
import networkx as nx


def _req(ok, *msg):
    """An explicit check (not assert, so that python -O cannot skip it)."""
    if not ok:
        print('REJECTED:', *msg, file=sys.stderr)
        sys.exit(1)

P, Q = 8, 3
OK = [[Q <= (a - b) % P <= P - Q for b in range(P)] for a in range(P)]


def maps_to_k83(G):
    n = G.number_of_nodes()
    order = sorted(G.nodes, key=lambda v: -G.degree(v))
    c = {}

    def go(t):
        if t == n:
            return True
        v = order[t]
        for k in ((0,) if t == 0 else range(P)):
            if all(w not in c or OK[k][c[w]] for w in G[v]):
                c[v] = k
                if go(t + 1):
                    return True
                del c[v]
        return False
    return go(0)


def every_3colouring_has_tight_cycle(G):
    V = sorted(G.nodes); n = len(V); idx = {v: i for i, v in enumerate(V)}
    E = [(idx[a], idx[b]) for a, b in G.edges]
    adj = [[] for _ in range(n)]
    for a, b in E:
        adj[a].append(b); adj[b].append(a)
    seen = False
    for rest in itertools.product(range(3), repeat=n - 1):
        c = (0,) + rest
        if any(c[a] == c[b] for a, b in E):
            continue
        seen = True
        out = [[u for u in adj[x] if (c[u] - c[x]) % 3 == 1] for x in range(n)]
        state = [0] * n
        cyc = False
        for s in range(n):
            if state[s] or cyc:
                continue
            stack = [(s, iter(out[s]))]; state[s] = 1
            while stack and not cyc:
                x, it = stack[-1]
                y = next(it, None)
                if y is None:
                    stack.pop(); state[x] = 2
                elif state[y] == 1:
                    cyc = True
                elif state[y] == 0:
                    state[y] = 1; stack.append((y, iter(out[y])))
        if not cyc:
            return False
    return seen


def dedup(graphs):
    buckets, out = {}, []
    for G in graphs:
        lst = buckets.setdefault(nx.weisfeiler_lehman_graph_hash(G, iterations=3), [])
        if not any(nx.is_isomorphic(G, H) for H in lst):
            lst.append(G); out.append(G)
    return out


level = [nx.empty_graph(1)]
counts = [1]
for n in range(2, 10):
    cand = []
    for G in level:
        V = list(G.nodes)
        for r in range(len(V) + 1):
            for S in itertools.combinations(V, r):
                if any(G.has_edge(a, b) for a, b in itertools.combinations(S, 2)):
                    continue
                H = G.copy(); H.add_node(n - 1); H.add_edges_from((s, n - 1) for s in S)
                cand.append(H)
    level = dedup(cand)
    counts.append(len(level))
    a = [not maps_to_k83(G) for G in level]
    b = [every_3colouring_has_tight_cycle(G) for G in level]
    _req(a == b, 'the two tests disagree on', n, 'vertices')
    print(f'n = {n}: {len(level)} triangle-free graphs up to isomorphism, {sum(a)} with chi_c = 3 (both tests agree)',
          flush=True)
_req(counts == [1, 2, 3, 7, 14, 38, 107, 410, 1897], 'unexpected counts', counts)
three = sorted((G for G, x in zip(level, a) if x), key=lambda G: G.number_of_edges())
print('the triangle-free graphs with 9 vertices and chi_c = 3:', [G.number_of_edges() for G in three], 'edges')
# the hexagon v0..v5 with the diagonals v_j v_{j+3} subdivided by m_j (vertices 6, 7, 8)
M = nx.cycle_graph(6)
for j in range(3):
    M.add_edges_from([(j, 6 + j), (6 + j, j + 3)])
M1 = M.copy(); M1.add_edge(6, 7)
M2 = M1.copy(); M2.add_edge(7, 8)
_req(len(three) == 3 and all(nx.is_isomorphic(G, H) for G, H in zip(three, [M, M1, M2])), 'unexpected graphs')
print('they are the hexagon with its long diagonals subdivided, and it with one or two edges between the midpoints')
