/*
 * SA reliability test: takes a specific D11 on stdin, searches for D12.
 * Input format: space-separated D11 elements on one line.
 *
 * Usage: echo "1 2 3 4 6 10 12 13 15 19 20 23 24 28 30 31 33 37 39 40 41 42" | ./sa_orbit_test <p> [trials] [iter]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_M 200

static int P, N_PARAM, M, N_TOTAL;
static int RED_THRESH, BLUE_THRESH;
static int D11_SIZE, D12_SIZE;

static int D11[MAX_M], D12[MAX_M], D12T[MAX_M];
static int delta_11[MAX_M], delta_12[MAX_M], delta_12T[MAX_M];
static int mod_add[MAX_M][MAX_M], mod_sub[MAX_M][MAX_M], mod_neg[MAX_M];

static int d12_mem[MAX_M], d12_nmem[MAX_M];
static int d12_cnt, d12_ncnt;
static int d12_mem_idx[MAX_M], d12_nmem_idx[MAX_M];

static unsigned long long rng_s[4];

static inline unsigned long long rotl(const unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline unsigned long long rng_next(void) {
    const unsigned long long result = rotl(rng_s[1] * 5, 7) * 9;
    const unsigned long long t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0]; rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2]; rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t; rng_s[3] = rotl(rng_s[3], 45);
    return result;
}

static void rng_seed(unsigned long long seed) {
    rng_s[0] = seed ^ 0x9E3779B97F4A7C15ULL;
    rng_s[1] = seed ^ 0x6C62272E07BB0142ULL;
    rng_s[2] = seed ^ 0xBF58476D1CE4E5B9ULL;
    rng_s[3] = seed ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) rng_next();
}

static inline int rng_int(int n) { return (int)(rng_next() % (unsigned long long)n); }
static inline double rng_double(void) { return (rng_next() >> 11) * 0x1.0p-53; }

static void init_mod_tables(void) {
    for (int a = 0; a < M; a++) {
        mod_neg[a] = (M - a) % M;
        for (int b = 0; b < M; b++) {
            mod_add[a][b] = (a + b) % M;
            mod_sub[a][b] = (a - b + M) % M;
        }
    }
}

static void compute_delta_full(const int *ind, int *delta) {
    for (int d = 0; d < M; d++) {
        int count = 0;
        for (int a = 0; a < M; a++)
            if (ind[a] && ind[mod_sub[a][d]]) count++;
        delta[d] = count;
    }
}

static inline void update_delta_inc(int *delta, const int *S, int rem, int add_el) {
    for (int d = 1; d < M; d++) {
        int change = S[mod_sub[add_el][d]] + S[mod_add[add_el][d]]
                   - S[mod_sub[rem][d]] - S[mod_add[rem][d]];
        if (mod_sub[rem][d] == add_el) change++;
        if (mod_add[rem][d] == add_el) change++;
        delta[d] += change;
    }
}

static void rebuild_d12_lists(void) {
    d12_cnt = 0; d12_ncnt = 0;
    for (int i = 0; i < M; i++) {
        if (D12[i]) { d12_mem_idx[i] = d12_cnt; d12_mem[d12_cnt++] = i; d12_nmem_idx[i] = -1; }
        else { d12_nmem_idx[i] = d12_ncnt; d12_nmem[d12_ncnt++] = i; d12_mem_idx[i] = -1; }
    }
}

static void d12_swap_to_mem(int e) {
    int ni = d12_nmem_idx[e]; int last = d12_nmem[d12_ncnt-1];
    d12_nmem[ni] = last; d12_nmem_idx[last] = ni; d12_ncnt--;
    d12_mem_idx[e] = d12_cnt; d12_mem[d12_cnt++] = e; d12_nmem_idx[e] = -1;
}

static void d12_swap_to_nmem(int e) {
    int mi = d12_mem_idx[e]; int last = d12_mem[d12_cnt-1];
    d12_mem[mi] = last; d12_mem_idx[last] = mi; d12_cnt--;
    d12_nmem_idx[e] = d12_ncnt; d12_nmem[d12_ncnt++] = e; d12_mem_idx[e] = -1;
}

static int compute_cost(void) {
    int d11_size = D11_SIZE, d12_size = D12_SIZE;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size, d2 = d22_size + d12_size;
    int cost = 0;
    for (int d = 1; d < M; d++) {
        int cv11 = delta_11[d] + delta_12[d];
        if (D11[d]) { int ex = cv11 - RED_THRESH; if (ex > 0) cost += ex; }
        else { int bc = (N_TOTAL-2) - 2*d1 + cv11; int ex = bc - BLUE_THRESH; if (ex > 0) cost += ex; }
        int d22d = delta_11[d] + (M-2-2*d11_size) + 2*D11[d];
        int cv22 = d22d + delta_12T[d];
        if (!D11[d]) { int ex = cv22 - RED_THRESH; if (ex > 0) cost += ex; }
        else { int bc = (N_TOTAL-2) - 2*d2 + cv22; int ex = bc - BLUE_THRESH; if (ex > 0) cost += ex; }
    }
    return cost;
}

static int verify_v1v2(void) {
    int d11_size = D11_SIZE, d12_size = D12_SIZE;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size, d2 = d22_size + d12_size;
    int D22[MAX_M];
    for (int i = 0; i < M; i++) D22[i] = (i > 0 && !D11[i]) ? 1 : 0;
    int cost = 0;
    for (int d = 0; d < M; d++) {
        int sigma = 0; for (int a = 0; a < M; a++) if (D11[a] && D12[mod_sub[d][a]]) sigma++;
        int delta = 0; for (int a = 0; a < M; a++) if (D12[a] && D22[mod_sub[a][d]]) delta++;
        int common = sigma + delta;
        if (D12[d]) { int ex = common - RED_THRESH; if (ex > 0) cost += ex; }
        else { int bc = (N_TOTAL-2) - d1 - d2 + common; int ex = bc - BLUE_THRESH; if (ex > 0) cost += ex; }
    }
    return cost;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: echo 'D11 elems' | %s <p> [trials] [iter]\n", argv[0]); return 1; }

    P = atoi(argv[1]);
    int sa_trials = argc > 2 ? atoi(argv[2]) : 20;
    int sa_iter = argc > 3 ? atoi(argv[3]) : 2000000;

    N_PARAM = (P + 1) / 2; M = P; N_TOTAL = 4 * N_PARAM - 2;
    RED_THRESH = N_PARAM - 2; BLUE_THRESH = N_PARAM - 1;
    D11_SIZE = (P + 1) / 2; D12_SIZE = (P - 1) / 2;

    init_mod_tables();

    /* Read D11 from stdin */
    memset(D11, 0, sizeof(int) * M);
    int actual_d11 = 0;
    int x;
    while (scanf("%d", &x) == 1) {
        if (x >= 0 && x < M) { D11[x] = 1; actual_d11++; }
    }
    printf("p=%d, |D11|=%d (expected %d)\n", P, actual_d11, D11_SIZE);
    compute_delta_full(D11, delta_11);

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    for (int trial = 0; trial < sa_trials; trial++) {
        rng_seed((unsigned long long)trial * 999983ULL + 17);

        memset(D12, 0, sizeof(int) * M);
        memset(D12T, 0, sizeof(int) * M);
        D12[0] = 1; D12T[0] = 1;

        int elems[MAX_M]; int ne = 0;
        for (int i = 1; i < M; i++) elems[ne++] = i;
        for (int i = ne-1; i > 0; i--) { int j = rng_int(i+1); int t = elems[i]; elems[i] = elems[j]; elems[j] = t; }
        for (int i = 0; i < D12_SIZE - 1; i++) { D12[elems[i]] = 1; D12T[mod_neg[elems[i]]] = 1; }
        rebuild_d12_lists();
        compute_delta_full(D12, delta_12);
        compute_delta_full(D12T, delta_12T);

        int cost = compute_cost();
        int best = cost;

        double T = 8.0, T_min = 0.1;
        double alpha = exp(log(T_min / T) / (double)sa_iter);

        for (int it = 0; it < sa_iter; it++) {
            if (cost == 0) {
                int v12 = verify_v1v2();
                if (v12 == 0) {
                    clock_gettime(CLOCK_MONOTONIC, &ts_now);
                    double el = (ts_now.tv_sec - ts_start.tv_sec) + (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
                    printf("  FOUND at trial=%d iter=%d (%.2fs)\n", trial, it, el);
                    printf("  D12=[");
                    for (int i = 0, f = 1; i < M; i++) if (D12[i]) { if (!f) printf(","); printf("%d", i); f = 0; }
                    printf("]\n");
                    return 0;
                }
            }

            if (d12_cnt <= 1 || d12_ncnt == 0) continue;
            int ri = rng_int(d12_cnt); int rem = d12_mem[ri];
            if (rem == 0) continue;
            int ai = rng_int(d12_ncnt); int add_el = d12_nmem[ai];

            D12[rem] = 0; D12[add_el] = 1;
            int rem_t = mod_neg[rem], add_t = mod_neg[add_el];
            D12T[rem_t] = 0; D12T[add_t] = 1;
            update_delta_inc(delta_12, D12, rem, add_el);
            update_delta_inc(delta_12T, D12T, rem_t, add_t);

            int nc = compute_cost();
            int dc = nc - cost;
            if (dc <= 0 || rng_double() < exp(-(double)dc / fmax(T, 0.01))) {
                cost = nc; if (cost < best) best = cost;
                d12_swap_to_nmem(rem); d12_swap_to_mem(add_el);
            } else {
                D12[rem] = 1; D12[add_el] = 0;
                D12T[rem_t] = 1; D12T[add_t] = 0;
                update_delta_inc(delta_12, D12, add_el, rem);
                update_delta_inc(delta_12T, D12T, add_t, rem_t);
            }
            T *= alpha;
        }

        clock_gettime(CLOCK_MONOTONIC, &ts_now);
        double el = (ts_now.tv_sec - ts_start.tv_sec) + (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
        printf("  trial %d: best_cost=%d (%.2fs)\n", trial, best, el);
    }

    printf("NOT FOUND after %d trials\n", sa_trials);
    return 1;
}
