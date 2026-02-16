/*
 * SA orbit scanner: sample random symmetric D11, search for valid D12.
 * Fixes D11 and only optimizes D12 with incremental delta updates.
 *
 * Usage: ./sa_orbit_scan <p> [n_orbits] [sa_trials] [sa_iter]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_M 600

static int P, N_PARAM, M, N_TOTAL;
static int RED_THRESH, BLUE_THRESH;
static int D11_SIZE, D12_SIZE;

/* State arrays */
static int D11[MAX_M];
static int D12[MAX_M];
static int D12T[MAX_M];

/* Autocorrelations */
static int delta_11[MAX_M];  /* Fixed for each orbit */
static int delta_12[MAX_M];
static int delta_12T[MAX_M];

/* Precomputed modular arithmetic */
static int mod_add[MAX_M][MAX_M];
static int mod_sub[MAX_M][MAX_M];
static int mod_neg[MAX_M];

/* Symmetric pairs */
static int pairs[MAX_M][2];
static int num_pairs;

/* D12 member/nonmember lists */
static int d12_mem[MAX_M], d12_nmem[MAX_M];
static int d12_cnt, d12_ncnt;
static int d12_mem_idx[MAX_M], d12_nmem_idx[MAX_M];

/* RNG (xoshiro256**) */
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

static inline int rng_int(int n) {
    return (int)(rng_next() % (unsigned long long)n);
}

static inline double rng_double(void) {
    return (rng_next() >> 11) * 0x1.0p-53;
}

static void init_mod_tables(void) {
    for (int a = 0; a < M; a++) {
        mod_neg[a] = (M - a) % M;
        for (int b = 0; b < M; b++) {
            mod_add[a][b] = (a + b) % M;
            mod_sub[a][b] = (a - b + M) % M;
        }
    }
}

static void build_pairs(void) {
    num_pairs = 0;
    for (int x = 1; x <= (M - 1) / 2; x++) {
        pairs[num_pairs][0] = x;
        pairs[num_pairs][1] = M - x;
        num_pairs++;
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
        if (D12[i]) {
            d12_mem_idx[i] = d12_cnt;
            d12_mem[d12_cnt++] = i;
            d12_nmem_idx[i] = -1;
        } else {
            d12_nmem_idx[i] = d12_ncnt;
            d12_nmem[d12_ncnt++] = i;
            d12_mem_idx[i] = -1;
        }
    }
}

static void d12_swap_to_mem(int elem) {
    int ni = d12_nmem_idx[elem];
    int last = d12_nmem[d12_ncnt - 1];
    d12_nmem[ni] = last; d12_nmem_idx[last] = ni;
    d12_ncnt--;
    d12_mem_idx[elem] = d12_cnt;
    d12_mem[d12_cnt++] = elem;
    d12_nmem_idx[elem] = -1;
}

static void d12_swap_to_nmem(int elem) {
    int mi = d12_mem_idx[elem];
    int last = d12_mem[d12_cnt - 1];
    d12_mem[mi] = last; d12_mem_idx[last] = mi;
    d12_cnt--;
    d12_nmem_idx[elem] = d12_ncnt;
    d12_nmem[d12_ncnt++] = elem;
    d12_mem_idx[elem] = -1;
}

/* Compute V1V1+V2V2 cost (D11 is fixed, delta_11 precomputed) */
static int compute_cost(void) {
    int d11_size = D11_SIZE;
    int d12_size = D12_SIZE;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;
    int cost = 0;

    for (int d = 1; d < M; d++) {
        int cv11 = delta_11[d] + delta_12[d];
        if (D11[d]) {
            int excess = cv11 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }

        int d22d = delta_11[d] + (M - 2 - 2 * d11_size) + 2 * D11[d];
        int cv22 = d22d + delta_12T[d];
        if (!D11[d]) {
            int excess = cv22 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d2 + cv22;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }
    return cost;
}

/* Full V1V2 verification */
static int verify_v1v2(void) {
    int d11_size = D11_SIZE;
    int d12_size = D12_SIZE;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int D22[MAX_M];
    for (int i = 0; i < M; i++)
        D22[i] = (i > 0 && !D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < M; d++) {
        int sigma = 0;
        for (int a = 0; a < M; a++)
            if (D11[a] && D12[mod_sub[d][a]]) sigma++;
        int delta = 0;
        for (int a = 0; a < M; a++)
            if (D12[a] && D22[mod_sub[a][d]]) delta++;
        int common = sigma + delta;

        if (D12[d]) {
            int excess = common - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - d1 - d2 + common;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }
    return cost;
}

/* Generate random symmetric D11 */
static void gen_random_d11(void) {
    memset(D11, 0, sizeof(int) * M);
    int half_select = D11_SIZE / 2;

    int perm[MAX_M];
    for (int i = 0; i < num_pairs; i++) perm[i] = i;
    for (int i = num_pairs - 1; i > 0; i--) {
        int j = rng_int(i + 1);
        int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }
    for (int i = 0; i < half_select; i++) {
        D11[pairs[perm[i]][0]] = 1;
        D11[pairs[perm[i]][1]] = 1;
    }
}

/* SA search for D12 with D11 fixed. Returns 1 if solution found, 0 otherwise. */
static int sa_search_d12(int trial_seed, int max_iter) {
    rng_seed((unsigned long long)trial_seed * 999983ULL + 17);

    /* Random D12: always include 0, then choose D12_SIZE-1 from {1,...,M-1} */
    memset(D12, 0, sizeof(int) * M);
    memset(D12T, 0, sizeof(int) * M);
    D12[0] = 1;
    D12T[0] = 1;

    int elems[MAX_M];
    int ne = 0;
    for (int i = 1; i < M; i++) elems[ne++] = i;
    for (int i = ne - 1; i > 0; i--) {
        int j = rng_int(i + 1);
        int tmp = elems[i]; elems[i] = elems[j]; elems[j] = tmp;
    }
    for (int i = 0; i < D12_SIZE - 1; i++) {
        D12[elems[i]] = 1;
        D12T[mod_neg[elems[i]]] = 1;
    }
    rebuild_d12_lists();

    compute_delta_full(D12, delta_12);
    compute_delta_full(D12T, delta_12T);

    int current_cost = compute_cost();
    int local_best = current_cost;

    double T = 8.0;
    double T_min = 0.0001;
    double alpha = exp(log(T_min / T) / (double)max_iter);

    for (int it = 0; it < max_iter; it++) {
        if (current_cost == 0) {
            int v12 = verify_v1v2();
            if (v12 == 0) return 1;
            /* V1V2 fails, continue searching */
        }

        int r = rng_int(100);

        if (r < 75) {
            /* Single swap in D12 (keep 0 fixed) */
            if (d12_cnt <= 1 || d12_ncnt == 0) continue;

            /* Pick random member != 0 */
            int ri = rng_int(d12_cnt);
            int rem = d12_mem[ri];
            if (rem == 0) continue;  /* Don't remove 0 */

            int ai = rng_int(d12_ncnt);
            int add_el = d12_nmem[ai];

            D12[rem] = 0; D12[add_el] = 1;
            int rem_t = mod_neg[rem], add_t = mod_neg[add_el];
            D12T[rem_t] = 0; D12T[add_t] = 1;

            update_delta_inc(delta_12, D12, rem, add_el);
            update_delta_inc(delta_12T, D12T, rem_t, add_t);

            int new_cost = compute_cost();
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double() < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                d12_swap_to_nmem(rem);
                d12_swap_to_mem(add_el);
                if (current_cost < local_best) local_best = current_cost;
            } else {
                D12[rem] = 1; D12[add_el] = 0;
                D12T[rem_t] = 1; D12T[add_t] = 0;
                update_delta_inc(delta_12, D12, add_el, rem);
                update_delta_inc(delta_12T, D12T, add_t, rem_t);
            }
        } else {
            /* Double swap */
            if (d12_cnt <= 2 || d12_ncnt < 2) continue;

            int ri1 = rng_int(d12_cnt);
            int rem1 = d12_mem[ri1];
            if (rem1 == 0) continue;
            int ri2 = rng_int(d12_cnt - 1);
            int rem2 = d12_mem[ri2 >= ri1 ? ri2 + 1 : ri2];
            if (rem2 == 0) continue;

            int ai1 = rng_int(d12_ncnt);
            int ai2 = rng_int(d12_ncnt - 1);
            int add1 = d12_nmem[ai1];
            int add2 = d12_nmem[ai2 >= ai1 ? ai2 + 1 : ai2];

            /* Save old deltas */
            int old_d12[MAX_M], old_d12T[MAX_M];
            memcpy(old_d12, delta_12, sizeof(int) * M);
            memcpy(old_d12T, delta_12T, sizeof(int) * M);

            D12[rem1] = 0; D12[rem2] = 0;
            D12[add1] = 1; D12[add2] = 1;
            D12T[mod_neg[rem1]] = 0; D12T[mod_neg[rem2]] = 0;
            D12T[mod_neg[add1]] = 1; D12T[mod_neg[add2]] = 1;

            compute_delta_full(D12, delta_12);
            compute_delta_full(D12T, delta_12T);

            int new_cost = compute_cost();
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double() < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                d12_swap_to_nmem(rem1); d12_swap_to_nmem(rem2);
                d12_swap_to_mem(add1); d12_swap_to_mem(add2);
                if (current_cost < local_best) local_best = current_cost;
            } else {
                D12[rem1] = 1; D12[rem2] = 1;
                D12[add1] = 0; D12[add2] = 0;
                D12T[mod_neg[rem1]] = 1; D12T[mod_neg[rem2]] = 1;
                D12T[mod_neg[add1]] = 0; D12T[mod_neg[add2]] = 0;
                memcpy(delta_12, old_d12, sizeof(int) * M);
                memcpy(delta_12T, old_d12T, sizeof(int) * M);
            }
        }

        T *= alpha;
    }

    /* Final check if we reached cost 0 */
    if (current_cost == 0) {
        int v12 = verify_v1v2();
        if (v12 == 0) return 1;
    }

    return 0;
}


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <p> [n_orbits] [sa_trials] [sa_iter]\n", argv[0]);
        return 1;
    }

    P = atoi(argv[1]);
    int n_orbits = argc > 2 ? atoi(argv[2]) : 500;
    int sa_trials = argc > 3 ? atoi(argv[3]) : 16;
    int sa_iter = argc > 4 ? atoi(argv[4]) : 2000000;

    N_PARAM = (P + 1) / 2;
    M = P;
    N_TOTAL = 4 * N_PARAM - 2;
    RED_THRESH = N_PARAM - 2;
    BLUE_THRESH = N_PARAM - 1;
    D11_SIZE = (P + 1) / 2;
    D12_SIZE = (P - 1) / 2;

    printf("SA orbit scan: p=%d, n=%d, m=%d\n", P, N_PARAM, M);
    printf("Thresholds: red=%d, blue=%d, |D11|=%d, |D12|=%d\n",
           RED_THRESH, BLUE_THRESH, D11_SIZE, D12_SIZE);
    printf("Orbits: %d, SA trials/orbit: %d, SA iter/trial: %d\n",
           n_orbits, sa_trials, sa_iter);
    fflush(stdout);

    init_mod_tables();
    build_pairs();

    /* Master RNG for orbit generation */
    rng_seed(42);

    int n_working = 0;
    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    for (int orbit = 0; orbit < n_orbits; orbit++) {
        /* Generate random symmetric D11 */
        gen_random_d11();
        compute_delta_full(D11, delta_11);

        int found = 0;
        for (int trial = 0; trial < sa_trials; trial++) {
            int seed = orbit * 10000 + trial;
            if (sa_search_d12(seed, sa_iter)) {
                found = 1;
                break;
            }
        }

        if (found) {
            n_working++;
            /* Print the working D11 */
            printf("  WORKING orbit %d: D11=[", orbit);
            for (int i = 0, first = 1; i < M; i++)
                if (D11[i]) { if (!first) printf(","); printf("%d", i); first = 0; }
            printf("] D12=[");
            for (int i = 0, first = 1; i < M; i++)
                if (D12[i]) { if (!first) printf(","); printf("%d", i); first = 0; }
            printf("]\n");
            fflush(stdout);
        }

        if ((orbit + 1) % 10 == 0 || found) {
            clock_gettime(CLOCK_MONOTONIC, &ts_now);
            double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                           (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
            double pw = (double)n_working / (orbit + 1);
            double ppw = P * pw;
            printf("[%d/%d] working=%d p_working=%.4f p*p_working=%.3f (%.0fs)\n",
                   orbit + 1, n_orbits, n_working, pw, ppw, elapsed);
            fflush(stdout);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                   (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
    double pw = (double)n_working / n_orbits;
    double ppw = P * pw;

    printf("\n========================================\n");
    printf("RESULTS: p=%d\n", P);
    printf("  Orbits tested:  %d\n", n_orbits);
    printf("  Working:        %d\n", n_working);
    printf("  p_working:      %.4f\n", pw);
    printf("  p * p_working:  %.4f\n", ppw);
    printf("  Time:           %.1fs\n", elapsed);
    printf("========================================\n");

    return 0;
}
