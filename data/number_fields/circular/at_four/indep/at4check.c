/*
 * Exhaustive check of Lemma F13, Corollary F14 and Corollary F13' (draft q3_at4.md)
 * on small finite abelian groups Gamma = Z/m_1 x ... x Z/m_d.
 *
 * For every symmetric generating S (0 not in S), up to automorphisms of Gamma, with Cay(Gamma,S)
 * 4-colourable, enumerate ALL proper 4-colourings c with c(0)=0 (colour rotations c+k give the rest)
 * and check, with integer arithmetic only:
 *   (F13a) if A(x,s,t) != B(x,s,t) then {A,B}={2,6} and the square x,x+s,x+s+t,x+t is a tight cycle
 *          (in the orientation given by A=2);
 *   (F13b) a colouring without tight squares (hence without tight 4-cycles) has A=B everywhere,
 *          L(s)+L(-s)=4n where L(s)=sum_x l(x,s), the map s -> L(s)/(4n) mod 1 is the restriction of a
 *          character, and the closed walks (ord(s) steps s) and (s,t,-(s+t)) have x-independent
 *          winding equal to the averaged sum, divisible by 4;
 *   (F14)  a colouring with acyclic tight digraph has L(s) in (n,3n) for all s, and kappa(S) > 1/4;
 *          (ii) kappa <= 1/4 => no colouring has acyclic tight digraph;
 *          (iii) kappa > 1/4 => some colouring has; the colouring floor(4 xi) of an optimal xi has;
 *   (F13') (some colouring has no tight square) <=> kappa(S) >= 1/4, and for EVERY character xi with
 *          xi(S) in [1/4,3/4] the colouring floor(4 xi) is proper and has no tight square;
 *   (hom)  independent test: Cay maps to K_{P/Q}, P=4Q-1, Q=floor((n+1)/4) (the largest fraction below 4
 *          with numerator <= n), iff kappa > 1/4.
 * Usage: at4check m_1 ... m_d
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 32
#define MAXS 32

static int d, m[4], stride[4], n, M;
static int addt[MAXN][MAXN], negt[MAXN], ordt[MAXN];
static int coord[MAXN][4];

static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

static void build_group(void) {
    n = 1; M = 1;
    for (int i = 0; i < d; i++) { stride[i] = n; n *= m[i]; M = M / gcd(M, m[i]) * m[i]; }
    for (int x = 0; x < n; x++) { int r = x; for (int i = 0; i < d; i++) { coord[x][i] = r % m[i]; r /= m[i]; } }
    for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) {
        int z = 0; for (int i = 0; i < d; i++) z += ((coord[x][i] + coord[y][i]) % m[i]) * stride[i];
        addt[x][y] = z;
    }
    for (int x = 0; x < n; x++) { int z = 0; for (int i = 0; i < d; i++) z += ((m[i] - coord[x][i]) % m[i]) * stride[i]; negt[x] = z; }
    /* order: invariant y = k*x, the loop stops at the least k >= 1 with k*x = 0 */
    for (int x = 0; x < n; x++) { int k = 1, y = x; while (y != 0) { y = addt[y][x]; k++; } ordt[x] = k; }
    for (int x = 1; x < n; x++) { /* verify */ int y = 0; for (int k = 0; k < ordt[x]; k++) y = addt[y][x];
        if (y != 0) { fprintf(stderr, "order bug\n"); exit(1); }
        y = 0; for (int k = 1; k < ordt[x]; k++) { y = addt[y][x]; if (y == 0) { fprintf(stderr, "order bug2\n"); exit(1);} } }
}

/* character value v_j(x) in [0,M): xi_j(x) = v_j(x)/M */
static int charval(int j, int x) {
    long v = 0; for (int i = 0; i < d; i++) v += (long)coord[j][i] * coord[x][i] * (M / m[i]);
    return (int)(v % M);
}
static int charv[MAXN][MAXN];

/* automorphisms */
static int naut; static unsigned char (*aut)[MAXN];

static void build_auts(void) {
    int cap = 1; for (int i = 0; i < d; i++) cap *= n;
    aut = malloc((size_t)cap * MAXN);
    naut = 0;
    int g[4];
    long total = 1; for (int i = 0; i < d; i++) total *= n;
    for (long t = 0; t < total; t++) {
        long r = t; int ok = 1;
        for (int i = 0; i < d; i++) { g[i] = (int)(r % n); r /= n; if (m[i] % ordt[g[i]] != 0) ok = 0; }
        if (!ok) continue;
        int img[MAXN]; int seen[MAXN]; memset(seen, 0, sizeof seen);
        for (int x = 0; x < n; x++) {
            int y = 0;
            for (int i = 0; i < d; i++) for (int k = 0; k < coord[x][i]; k++) y = addt[y][g[i]];
            img[x] = y; if (seen[y]) { ok = 0; break; } seen[y] = 1;
        }
        if (!ok) continue;
        for (int x = 0; x < n; x++) aut[naut][x] = (unsigned char)img[x];
        naut++;
    }
}

/* classes {x,-x} */
static int ncls, clsrep[MAXN], clsof[MAXN];

/* current S */
static int ns, S[MAXS], sidx[MAXN], negk[MAXS];

static int generates(void) {
    int seen[MAXN] = {0}, q[MAXN], h = 0, tl = 0; seen[0] = 1; q[tl++] = 0;
    while (h < tl) { int x = q[h++]; for (int k = 0; k < ns; k++) { int y = addt[x][S[k]]; if (!seen[y]) { seen[y] = 1; q[tl++] = y; } } }
    return tl == n;
}

/* colouring enumeration */
static int vorder[MAXN], col[MAXN];
static long long ncol;

/* per-S statistics */
static long long n_acyc, n_not4, n_nosq, n_abeq, v_f13a, v_not4_abneq, v_nosq_abneq, v_char, v_walk,
    n_abneq_charfail, v_acyc, n_t4_abeq;
static int kbest, jbest; /* kappa = kbest/M */

static int lt[MAXN][MAXS];

static int is_char_restriction(const long long *L) {
    /* exists character j with L[k]/(4n) == v_j(S[k])/M mod 1 for all k, i.e. L[k]*M == 4n*v (mod 4nM) */
    long long mod = 4LL * n * M;
    for (int j = 0; j < n; j++) {
        int ok = 1;
        for (int k = 0; k < ns && ok; k++) {
            long long lhs = (L[k] * M) % mod, rhs = (4LL * n * charv[j][S[k]]) % mod;
            if (lhs != rhs) ok = 0;
        }
        if (ok) return 1;
    }
    return 0;
}

static int acyclic_tight(void) {
    int indeg[MAXN] = {0};
    for (int x = 0; x < n; x++) for (int k = 0; k < ns; k++) if (lt[x][k] == 1) indeg[addt[x][S[k]]]++;
    int q[MAXN], h = 0, tl = 0;
    for (int x = 0; x < n; x++) if (!indeg[x]) q[tl++] = x;
    while (h < tl) { int x = q[h++]; for (int k = 0; k < ns; k++) if (lt[x][k] == 1) { int y = addt[x][S[k]]; if (--indeg[y] == 0) q[tl++] = y; } }
    return tl == n;
}

static void check_colouring(void) {
    ncol++;
    for (int x = 0; x < n; x++) for (int k = 0; k < ns; k++) {
        int v = (col[addt[x][S[k]]] - col[x]) & 3;
        if (v == 0) { fprintf(stderr, "improper\n"); exit(1); }
        lt[x][k] = v;
    }
    /* tight closed walks of length 4 (automatically 4-cycles: colours c,c+1,c+2,c+3) and tight squares */
    int t4 = 0, sq = 0;
    for (int x = 0; x < n && !(t4 && sq); x++) for (int k1 = 0; k1 < ns; k1++) if (lt[x][k1] == 1) {
        int y = addt[x][S[k1]];
        for (int k2 = 0; k2 < ns; k2++) if (lt[y][k2] == 1) {
            int z = addt[y][S[k2]];
            for (int k3 = 0; k3 < ns; k3++) if (lt[z][k3] == 1) {
                int w = addt[z][S[k3]];
                int dlt = addt[x][negt[w]]; /* x - w */
                int k4 = sidx[dlt];
                if (k4 >= 0 && lt[w][k4] == 1) { t4 = 1; if (k3 == negk[k1]) sq = 1; }
            }
        }
    }
    /* A = B check */
    int abneq = 0, viol = 0;
    for (int x = 0; x < n; x++) for (int k1 = 0; k1 < ns; k1++) for (int k2 = 0; k2 < ns; k2++) {
        int s = S[k1], t = S[k2];
        int xs = addt[x][s], xt = addt[x][t], xst = addt[xs][t];
        int A = lt[x][k1] + lt[xs][k2], B = lt[x][k2] + lt[xt][k1];
        if (A != B) {
            abneq = 1;
            int good = 0;
            if (A == 2 && B == 6) good = (lt[x][k1] == 1 && lt[xs][k2] == 1 && lt[xst][negk[k1]] == 1 && lt[xt][negk[k2]] == 1);
            else if (A == 6 && B == 2) good = (lt[x][k2] == 1 && lt[xt][k1] == 1 && lt[xst][negk[k2]] == 1 && lt[xs][negk[k1]] == 1);
            if (!good) viol = 1;
        }
    }
    if (viol) v_f13a++;
    if (!t4) n_not4++;
    if (!sq) n_nosq++;
    if (!abneq) n_abeq++;
    if (!t4 && abneq) v_not4_abneq++;
    if (!sq && abneq) v_nosq_abneq++;
    if (t4 && !abneq) n_t4_abeq++;
    long long L[MAXS];
    for (int k = 0; k < ns; k++) { L[k] = 0; for (int x = 0; x < n; x++) L[k] += lt[x][k]; }
    if (!abneq) {
        int bad = 0;
        for (int k = 0; k < ns; k++) if (L[k] + L[negk[k]] != 4LL * n) bad = 1;
        if (!is_char_restriction(L)) bad = 1;
        if (bad) v_char++;
        /* closed walks */
        int wbad = 0;
        for (int k = 0; k < ns; k++) {
            int s = S[k], o = ordt[s];
            long long lam0 = -1;
            for (int x = 0; x < n; x++) {
                long long lam = 0; int y = x;
                for (int i = 0; i < o; i++) { lam += lt[y][k]; y = addt[y][s]; }
                if (y != x) { fprintf(stderr, "walk not closed\n"); exit(1); }
                if (lam0 < 0) lam0 = lam; else if (lam != lam0) wbad = 1;
                if (lam % 4) wbad = 1;
            }
            if ((long long)n * lam0 != (long long)o * L[k]) wbad = 1;
        }
        for (int k1 = 0; k1 < ns; k1++) for (int k2 = 0; k2 < ns; k2++) {
            int u = addt[S[k1]][S[k2]];
            if (u == 0 || sidx[u] < 0) continue;
            int k3 = sidx[negt[u]];
            long long lam0 = -1;
            for (int x = 0; x < n; x++) {
                int y1 = addt[x][S[k1]], y2 = addt[y1][S[k2]];
                long long lam = lt[x][k1] + lt[y1][k2] + lt[y2][k3];
                if (addt[y2][S[k3]] != x) { fprintf(stderr, "walk3 not closed\n"); exit(1); }
                if (lam0 < 0) lam0 = lam; else if (lam != lam0) wbad = 1;
                if (lam % 4) wbad = 1;
            }
            if ((long long)n * lam0 != L[k1] + L[k2] + L[k3]) wbad = 1;
        }
        if (wbad) v_walk++;
    } else {
        if (!is_char_restriction(L)) n_abneq_charfail++;
    }
    int acyc = acyclic_tight();
    if (acyc) {
        n_acyc++;
        int bad = 0;
        if (t4 || sq || abneq) bad = 1;
        for (int k = 0; k < ns; k++) if (!(L[k] > n && L[k] < 3LL * n)) bad = 1;
        if (!(4 * kbest > M)) bad = 1;
        if (bad) v_acyc++;
    }
}

static int adjlist[MAXN][MAXS], adjn[MAXN], posv[MAXN];

static void rec(int i) {
    if (i == n) { check_colouring(); return; }
    int v = vorder[i];
    for (int c = 0; c < 4; c++) {
        if (i == 0 && c != 0) break;
        int ok = 1;
        for (int a = 0; a < adjn[v] && ok; a++) { int u = adjlist[v][a]; if (posv[u] < i && col[u] == c) ok = 0; }
        if (!ok) continue;
        col[v] = c; rec(i + 1);
    }
}

/* homomorphism to K_{P/Q} by backtracking (independent of the above) */
static int P_, Q_, hc[MAXN];
static int hrec(int i) {
    if (i == n) return 1;
    int v = vorder[i];
    for (int c = 0; c < P_; c++) {
        if (i == 0 && c != 0) break;
        int ok = 1;
        for (int a = 0; a < adjn[v] && ok; a++) { int u = adjlist[v][a]; if (posv[u] < i) { int dd = ((c - hc[u]) % P_ + P_) % P_; if (dd < Q_ || dd > P_ - Q_) ok = 0; } }
        if (!ok) continue;
        hc[v] = c; if (hrec(i + 1)) return 1;
    }
    return 0;
}

static void print_elem(int x) {
    printf("(");
    for (int i = 0; i < d; i++) printf("%s%d", i ? "," : "", coord[x][i]);
    printf(")");
}

int main(int argc, char **argv) {
    d = argc - 1; if (d < 1 || d > 4) { fprintf(stderr, "usage\n"); return 1; }
    for (int i = 0; i < d; i++) m[i] = atoi(argv[i + 1]);
    build_group();
    if (n > MAXN) { fprintf(stderr, "too big\n"); return 1; }
    for (int j = 0; j < n; j++) for (int x = 0; x < n; x++) charv[j][x] = charval(j, x);
    build_auts();
    ncls = 0;
    for (int x = 1; x < n; x++) { int r = x < negt[x] ? x : negt[x]; if (r == x) { clsrep[ncls] = x; clsof[x] = ncls; ncls++; } }
    for (int x = 1; x < n; x++) clsof[x] = clsof[x < negt[x] ? x : negt[x]];
    int nmask = 1 << ncls;
    unsigned char *visited = calloc(nmask, 1);
    /* class permutation of each automorphism */
    int (*cperm)[MAXN] = malloc(sizeof(int[MAXN]) * naut);
    for (int a = 0; a < naut; a++) for (int c = 0; c < ncls; c++) cperm[a][c] = clsof[aut[a][clsrep[c]]];
    printf("# group Z/%d", m[0]); for (int i = 1; i < d; i++) printf(" x Z/%d", m[i]);
    printf("  n=%d  exponent M=%d  |Aut|=%d  classes=%d\n", n, M, naut, ncls);
    long long tot_S = 0, tot_col = 0, tot_bad = 0, nS4 = 0;
    for (int mask = 1; mask < nmask; mask++) {
        if (visited[mask]) continue;
        for (int a = 0; a < naut; a++) { int im = 0; for (int c = 0; c < ncls; c++) if (mask >> c & 1) im |= 1 << cperm[a][c]; visited[im] = 1; }
        ns = 0; for (int x = 0; x < n; x++) sidx[x] = -1;
        for (int x = 1; x < n; x++) if (mask >> clsof[x] & 1) { sidx[x] = ns; S[ns++] = x; }
        for (int k = 0; k < ns; k++) negk[k] = sidx[negt[S[k]]];
        if (!generates()) continue;
        tot_S++;
        /* kappa */
        kbest = -1; jbest = -1;
        for (int j = 0; j < n; j++) {
            int mn = M;
            for (int k = 0; k < ns; k++) { int v = charv[j][S[k]]; int dd = v < M - v ? v : M - v; if (dd < mn) mn = dd; }
            if (mn > kbest) { kbest = mn; jbest = j; }
        }
        /* adjacency, BFS order */
        for (int x = 0; x < n; x++) { adjn[x] = 0; for (int k = 0; k < ns; k++) { int y = addt[x][S[k]]; int dup = 0; for (int a = 0; a < adjn[x]; a++) if (adjlist[x][a] == y) dup = 1; if (!dup) adjlist[x][adjn[x]++] = y; } }
        { int seen[MAXN] = {0}, h = 0, tl = 0; seen[0] = 1; vorder[tl++] = 0;
          while (h < tl) { int x = vorder[h++]; for (int a = 0; a < adjn[x]; a++) { int y = adjlist[x][a]; if (!seen[y]) { seen[y] = 1; vorder[tl++] = y; } } }
          for (int i = 0; i < n; i++) posv[vorder[i]] = i; }
        ncol = n_acyc = n_not4 = n_nosq = n_abeq = v_f13a = v_not4_abneq = v_nosq_abneq = v_char = v_walk = n_abneq_charfail = v_acyc = n_t4_abeq = 0;
        rec(0);
        if (ncol == 0) {
            printf("S=");
            for (int k = 0; k < ns; k++) print_elem(S[k]);
            printf(" |S|=%d kappa=%d/%d not4col\n", ns, kbest, M);
            continue;
        }
        nS4++;
        tot_col += ncol;
        /* (ii)/(iii) */
        int bad = 0; char why[512]; why[0] = 0;
        if (v_f13a) { bad = 1; strcat(why, " F13a"); }
        if (v_not4_abneq) { bad = 1; strcat(why, " noT4-but-AneB"); }
        if (v_nosq_abneq) { bad = 1; strcat(why, " noSq-but-AneB"); }
        if (v_char) { bad = 1; strcat(why, " char"); }
        if (v_walk) { bad = 1; strcat(why, " walk"); }
        if (v_acyc) { bad = 1; strcat(why, " acyc"); }
        int kap_gt = 4 * kbest > M, kap_ge = 4 * kbest >= M;
        if (!kap_gt && n_acyc) { bad = 1; strcat(why, " (ii)"); }
        if (kap_gt && !n_acyc) { bad = 1; strcat(why, " (iii)"); }
        if (kap_gt) {
            /* floor(4 xi) colouring of the optimal character: proper and acyclic */
            for (int x = 0; x < n; x++) col[x] = (4 * charv[jbest][x]) / M;
            for (int x = 0; x < n; x++) for (int k = 0; k < ns; k++) { int v = (col[addt[x][S[k]]] - col[x]) & 3; if (!v) { bad = 1; strcat(why, " floorimproper"); } lt[x][k] = v; }
            if (!acyclic_tight()) { bad = 1; strcat(why, " floor-not-acyclic"); }
        }
        /* F13': some colouring without tight square iff kappa >= 1/4 */
        if ((n_nosq > 0) != kap_ge) { bad = 1; strcat(why, " F13'"); }
        int nchar_ge = 0;
        for (int j = 0; j < n; j++) {
            int mn = M;
            for (int k = 0; k < ns; k++) { int v = charv[j][S[k]]; int dd = v < M - v ? v : M - v; if (dd < mn) mn = dd; }
            if (4 * mn < M) continue;
            nchar_ge++;
            for (int x = 0; x < n; x++) col[x] = (4 * charv[j][x]) / M;
            int improper = 0, sqf = 0;
            for (int x = 0; x < n; x++) for (int k = 0; k < ns; k++) { int v = (col[addt[x][S[k]]] - col[x]) & 3; if (!v) improper = 1; lt[x][k] = v; }
            if (!improper)
                for (int x = 0; x < n; x++) for (int k1 = 0; k1 < ns; k1++) for (int k2 = 0; k2 < ns; k2++) {
                    int xs = addt[x][S[k1]], xst = addt[xs][S[k2]], xt = addt[x][S[k2]];
                    if (lt[x][k1] == 1 && lt[xs][k2] == 1 && lt[xst][negk[k1]] == 1 && lt[xt][negk[k2]] == 1) sqf = 1;
                }
            if (improper || sqf) { bad = 1; strcat(why, " F13'-if"); break; }
        }
        /* independent hom test */
        Q_ = (n + 1) / 4; if (Q_ < 1) Q_ = 1; P_ = 4 * Q_ - 1;
        int hom = hrec(0);
        if (hom != kap_gt) { bad = 1; strcat(why, " hom"); }
        if (bad) tot_bad++;
        printf("S=");
        for (int k = 0; k < ns; k++) print_elem(S[k]);
        printf(" |S|=%d kappa=%d/%d col=%lld acyc=%lld noT4=%lld noSq=%lld ABeq=%lld T4butABeq=%lld ABneq_charfail=%lld chars_ge=%d hom%d/%d=%d %s%s\n",
               ns, kbest, M, ncol, n_acyc, n_not4, n_nosq, n_abeq, n_t4_abeq, n_abneq_charfail, nchar_ge, P_, Q_, hom,
               bad ? "BAD:" : "ok", why);
        fflush(stdout);
    }
    printf("# summary: generating S classes=%lld, 4-colourable=%lld, colourings(c(0)=0)=%lld, BAD=%lld\n", tot_S, nS4, tot_col, tot_bad);
    return 0;
}
