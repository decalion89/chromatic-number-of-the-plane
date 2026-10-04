/* Referee check 7: triangle-free graphs on n labelled vertices, counted up to isomorphism by Burnside's lemma
   (no isomorphism test, no canonical form):  #classes = (1/n!) sum over permutations s of Fix(s),
   where Fix(s) = number of s-invariant labelled triangle-free graphs; the same for the isomorphism-invariant
   property "no homomorphism to K_{8/3}".  For each cycle type one representative s is used, weighted by the size
   of its conjugacy class.  s-invariant graphs are unions of edge orbits; they are enumerated by backtracking over
   the orbits (include / exclude), rejecting an edge whose ends already have a common neighbour.  For s = identity
   this is an edge-by-edge enumeration of all labelled triangle-free graphs; there the graphs without a
   homomorphism to K_{8/3} are also counted by number of edges, tested for 3-colourability, and the maximal
   triangle-free graphs are counted.
   usage: tf_burnside n                                                                                       */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n;
static unsigned adj[16];
static int col[16], ord_[16], root_[16];
static int P, Q;
static unsigned allowed[16];

static void set_target(int p, int q) {
    P = p; Q = q;
    for (int k = 0; k < P; k++) {
        allowed[k] = 0;
        for (int j = 0; j < P; j++) { int dd = ((j - k) % P + P) % P; if (Q <= dd && dd <= P - Q) allowed[k] |= 1u << j; }
    }
}
static void bfs_order(void) {
    int seen = 0, t = 0;
    for (int s = 0; s < n; s++) {
        if (seen >> s & 1) continue;
        int head = t; ord_[t] = s; root_[t] = 1; t++; seen |= 1 << s;
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
    if (root_[t]) dom = 1u;
    else {
        dom = (1u << P) - 1;
        for (int s = 0; s < t; s++) { int w = ord_[s]; if (adj[v] >> w & 1) dom &= allowed[col[w]]; }
    }
    for (int k = 0; k < P; k++) if (dom >> k & 1) { col[v] = k; if (bt(t + 1)) return 1; }
    return 0;
}
static int hom(int p, int q) { set_target(p, q); bfs_order(); return bt(0); }

/* edge orbits of the current permutation */
static int norb, orbsz[64], orbi[64][64], orbj[64][64];
static long long fix_all, fix_nohom;
static int identity_mode;
static long long id_byedges[64], id_nohom_byedges[64], id_nohom_not3col, id_maximal, id_maximal_nohom;
static int nedges;
static int ei_all[64], ej_all[64], ne_all;

static void leaf(void) {
    fix_all++;
    int nh = !hom(8, 3);
    if (nh) fix_nohom++;
    if (identity_mode) {
        id_byedges[nedges]++;
        if (nh) { id_nohom_byedges[nedges]++; if (!hom(3, 1)) id_nohom_not3col++; }
        int mx = 1;
        for (int e = 0; e < ne_all && mx; e++) {
            int a = ei_all[e], b = ej_all[e];
            if (!(adj[a] >> b & 1) && !(adj[a] & adj[b])) mx = 0;
        }
        if (mx) { id_maximal++; if (nh) id_maximal_nohom++; }
    }
}

static void rec(int o) {
    if (o == norb) { leaf(); return; }
    rec(o + 1);                                  /* exclude orbit o */
    /* include orbit o: add its edges one by one, each must not close a triangle */
    int added = 0, ok = 1;
    for (int k = 0; k < orbsz[o]; k++) {
        int a = orbi[o][k], b = orbj[o][k];
        if (adj[a] & adj[b]) { ok = 0; break; }
        adj[a] |= 1u << b; adj[b] |= 1u << a; added++;
    }
    if (ok) { nedges += orbsz[o]; rec(o + 1); nedges -= orbsz[o]; }
    for (int k = 0; k < added; k++) { int a = orbi[o][k], b = orbj[o][k]; adj[a] &= ~(1u << b); adj[b] &= ~(1u << a); }
}

static void orbits(const int *sg) {
    int done[16][16]; memset(done, 0, sizeof done);
    norb = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        if (done[i][j]) continue;
        int a = i, b = j, k = 0;
        while (!done[a][b]) {
            done[a][b] = done[b][a] = 1;
            orbi[norb][k] = a; orbj[norb][k] = b; k++;
            a = sg[a]; b = sg[b];
        }
        orbsz[norb++] = k;
    }
}

static long long fact(int k) { long long r = 1; for (int i = 2; i <= k; i++) r *= i; return r; }

static int parts[16], np;
static long long total_w_all, total_w_nohom, total_class;

static void do_partition(void) {
    int sg[16], v = 0;
    for (int c = 0; c < np; c++) {
        int L = parts[c];
        for (int t = 0; t < L; t++) sg[v + t] = v + (t + 1) % L;
        v += L;
    }
    long long denom = 1; int mult[16] = {0};
    for (int c = 0; c < np; c++) { denom *= parts[c]; mult[parts[c]]++; }
    for (int L = 1; L <= n; L++) denom *= fact(mult[L]);
    long long cls = fact(n) / denom;
    orbits(sg);
    identity_mode = (np == n);
    fix_all = fix_nohom = 0; nedges = 0;
    for (int u = 0; u < n; u++) adj[u] = 0;
    rec(0);
    printf("  cycle type");
    for (int c = 0; c < np; c++) printf(" %d", parts[c]);
    printf(": class size %lld, edge orbits %d, Fix = %lld, Fix(no hom to K_8/3) = %lld\n", cls, norb, fix_all, fix_nohom);
    fflush(stdout);
    total_w_all += cls * fix_all; total_w_nohom += cls * fix_nohom; total_class += cls;
}

static void gen(int rem, int maxp) {
    if (rem == 0) { do_partition(); return; }
    for (int p = (rem < maxp ? rem : maxp); p >= 1; p--) { parts[np++] = p; gen(rem - p, p); np--; }
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    ne_all = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { ei_all[ne_all] = i; ej_all[ne_all] = j; ne_all++; }
    printf("n = %d\n", n);
    gen(n, n);
    long long nf = fact(n);
    printf("sum of class sizes = %lld (n! = %lld)\n", total_class, nf);
    printf("weighted sum %lld, / n! = %lld remainder %lld  => triangle-free graphs up to isomorphism\n",
           total_w_all, total_w_all / nf, total_w_all % nf);
    printf("weighted sum %lld, / n! = %lld remainder %lld  => ... without a homomorphism to K_8/3\n",
           total_w_nohom, total_w_nohom / nf, total_w_nohom % nf);
    printf("identity: labelled triangle-free graphs by number of edges:");
    for (int e = 0; e <= ne_all; e++) if (id_byedges[e]) printf(" %d:%lld", e, id_byedges[e]);
    printf("\nidentity: labelled triangle-free graphs without a homomorphism to K_8/3, by number of edges:");
    for (int e = 0; e <= ne_all; e++) if (id_nohom_byedges[e]) printf(" %d:%lld", e, id_nohom_byedges[e]);
    printf("\nidentity: of these, not 3-colourable: %lld; maximal triangle-free labelled graphs: %lld, of which "
           "without a homomorphism to K_8/3: %lld\n", id_nohom_not3col, id_maximal, id_maximal_nohom);
    return 0;
}
