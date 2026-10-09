/* Referee check 6: every triangle-free graph on n <= 8 labelled vertices maps to K_{8/3}.
   Brute force over ALL 2^(n(n-1)/2) labelled graphs (no reduction to maximal graphs, no vertex-by-vertex
   generation): count triangle-free graphs, maximal ones, non-bipartite maximal ones, and test a homomorphism to
   K_{8/3} on EVERY triangle-free graph.  Also tests 3-colourability (K_3) of every triangle-free graph.
   usage: tf_small n                                                                                         */
#include <stdio.h>
#include <stdlib.h>

static int n;
static unsigned adj[16];
static int col[16], ord_[16], root_[16];

/* homomorphism to K_{P/Q} by backtracking along a BFS order; first vertex of each component coloured 0 */
static int P, Q;
static unsigned allowed[16]; /* allowed[k] = colours j with Q <= (j-k) mod P <= P-Q */

static void bfs_order(void) {
    int seen = 0, t = 0;
    for (int s = 0; s < n; s++) {
        if (seen >> s & 1) continue;
        int head = t;
        ord_[t] = s; root_[t] = 1; t++; seen |= 1 << s;
        while (head < t) {
            int v = ord_[head++];
            for (int w = 0; w < n; w++)
                if ((adj[v] >> w & 1) && !(seen >> w & 1)) { ord_[t] = w; root_[t] = 0; t++; seen |= 1 << w; }
        }
    }
}

static int bt(int t) {
    if (t == n) return 1;
    int v = ord_[t];
    unsigned dom;
    if (root_[t]) dom = 1u;           /* colour 0 */
    else {
        dom = (1u << P) - 1;
        for (int s = 0; s < t; s++) { int w = ord_[s]; if (adj[v] >> w & 1) dom &= allowed[col[w]]; }
    }
    for (int k = 0; k < P; k++) if (dom >> k & 1) { col[v] = k; if (bt(t + 1)) return 1; }
    col[v] = -1;
    return 0;
}

static int hom(int p, int q) {
    P = p; Q = q;
    for (int k = 0; k < P; k++) {
        allowed[k] = 0;
        for (int j = 0; j < P; j++) { int dd = ((j - k) % P + P) % P; if (Q <= dd && dd <= P - Q) allowed[k] |= 1u << j; }
    }
    for (int v = 0; v < n; v++) col[v] = -1;
    bfs_order();
    return bt(0);
}

static int bipartite(void) { return hom(2, 1); }

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    int ne = n * (n - 1) / 2, ei[64], ej[64], k = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { ei[k] = i; ej[k] = j; k++; }
    long long tf = 0, maximal = 0, nonbip_max = 0, nohom83 = 0, no3col = 0, nohom52 = 0, maxnohom52 = 0;
    long long byedges[64] = {0};
    for (unsigned long long m = 0; m < (1ull << ne); m++) {
        for (int v = 0; v < n; v++) adj[v] = 0;
        for (int e = 0; e < ne; e++) if (m >> e & 1) { adj[ei[e]] |= 1u << ej[e]; adj[ej[e]] |= 1u << ei[e]; }
        int ok = 1;
        for (int e = 0; e < ne && ok; e++) if ((m >> e & 1) && (adj[ei[e]] & adj[ej[e]])) ok = 0;
        if (!ok) continue;
        tf++; byedges[__builtin_popcountll(m)]++;
        if (!hom(8, 3)) nohom83++;
        if (!hom(3, 1)) no3col++;
        if (!hom(5, 2)) nohom52++;
        int mx = 1;
        for (int e = 0; e < ne && mx; e++) if (!(m >> e & 1) && !(adj[ei[e]] & adj[ej[e]])) mx = 0;
        if (mx) { maximal++; if (!bipartite()) nonbip_max++; if (!hom(5, 2)) maxnohom52++; }
    }
    printf("n=%d: labelled graphs %llu, triangle-free %lld, maximal triangle-free %lld, maximal non-bipartite %lld\n",
           n, 1ull << ne, tf, maximal, nonbip_max);
    printf("  triangle-free without a homomorphism to K_8/3: %lld; not 3-colourable: %lld; "
           "without a homomorphism to K_5/2: %lld (maximal ones: %lld)\n", nohom83, no3col, nohom52, maxnohom52);
    return 0;
}
