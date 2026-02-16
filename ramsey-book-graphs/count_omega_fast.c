/*
 * Fast omega counter: count distinct D12-orbits for a given D11 at p=43.
 * Uses incremental delta updates from sa_orbit_scan.c for speed.
 * Reads D11 elements from command line args.
 *
 * Usage: ./count_omega_fast <n_trials> <d11_elem1> <d11_elem2> ...
 * Compile: gcc -O3 -o count_omega_fast count_omega_fast.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#ifndef P_VAL
#define P_VAL 43
#endif

#define P P_VAL
#define D11_SIZE ((P+1)/2)
#define D12_SIZE ((P-1)/2)
#define N_PARAM  ((P+1)/2)
#define N_TOTAL  (4*N_PARAM - 2)
#define RED_THRESH  (N_PARAM - 2)
#define BLUE_THRESH (N_PARAM - 1)

/* State arrays */
static int D11[P];
static int D12[P];
static int D12T[P];

/* Autocorrelations */
static int delta_11[P];
static int delta_12[P];
static int delta_12T[P];

/* Precomputed modular arithmetic */
static int mod_add[P][P];
static int mod_sub[P][P];
static int mod_neg[P];

/* D12 member/nonmember lists */
static int d12_mem[P], d12_nmem[P];
static int d12_cnt, d12_ncnt;
static int d12_mem_idx[P], d12_nmem_idx[P];

/* RNG (xoshiro256**) */
static unsigned long long rng_s[4];

static inline unsigned long long rotl64(const unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline unsigned long long rng_next(void) {
    const unsigned long long result = rotl64(rng_s[1] * 5, 7) * 9;
    const unsigned long long t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0]; rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2]; rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t; rng_s[3] = rotl64(rng_s[3], 45);
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
    for (int a = 0; a < P; a++) {
        mod_neg[a] = (P - a) % P;
        for (int b = 0; b < P; b++) {
            mod_add[a][b] = (a + b) % P;
            mod_sub[a][b] = (a - b + P) % P;
        }
    }
}

static void compute_delta_full(const int *ind, int *delta) {
    for (int d = 0; d < P; d++) {
        int count = 0;
        for (int a = 0; a < P; a++)
            if (ind[a] && ind[mod_sub[a][d]]) count++;
        delta[d] = count;
    }
}

static inline void update_delta_inc(int *delta, const int *S, int rem, int add_el) {
    for (int d = 1; d < P; d++) {
        int change = S[mod_sub[add_el][d]] + S[mod_add[add_el][d]]
                   - S[mod_sub[rem][d]] - S[mod_add[rem][d]];
        if (mod_sub[rem][d] == add_el) change++;
        if (mod_add[rem][d] == add_el) change++;
        delta[d] += change;
    }
}

static void rebuild_d12_lists(void) {
    d12_cnt = 0; d12_ncnt = 0;
    for (int i = 0; i < P; i++) {
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

/* Compute V1V1+V2V2 cost */
static int compute_cost(void) {
    int d11_size = D11_SIZE;
    int d12_size = D12_SIZE;
    int d22_size = P - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;
    int cost = 0;

    for (int d = 1; d < P; d++) {
        int cv11 = delta_11[d] + delta_12[d];
        if (D11[d]) {
            int excess = cv11 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }

        int d22d = delta_11[d] + (P - 2 - 2 * d11_size) + 2 * D11[d];
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
    int d22_size = P - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int D22[P];
    for (int i = 0; i < P; i++)
        D22[i] = (i > 0 && !D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < P; d++) {
        int sigma = 0;
        for (int a = 0; a < P; a++)
            if (D11[a] && D12[mod_sub[d][a]]) sigma++;
        int delta_val = 0;
        for (int a = 0; a < P; a++)
            if (D12[a] && D22[mod_sub[a][d]]) delta_val++;
        int common = sigma + delta_val;

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

/* SA search for D12. Returns 1 if valid solution found. */
static int sa_search_d12(int trial_seed, int max_iter) {
    rng_seed((unsigned long long)trial_seed * 999983ULL + 17);

    /* Random D12: always include 0, then choose D12_SIZE-1 from {1,...,P-1} */
    memset(D12, 0, sizeof(D12));
    memset(D12T, 0, sizeof(D12T));
    D12[0] = 1;
    D12T[0] = 1;

    int elems[P];
    int ne = 0;
    for (int i = 1; i < P; i++) elems[ne++] = i;
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

    double T = 8.0;
    double T_min = 0.0001;
    double alpha = exp(log(T_min / T) / (double)max_iter);

    for (int it = 0; it < max_iter; it++) {
        if (current_cost == 0) {
            int v12 = verify_v1v2();
            if (v12 == 0) return 1;
        }

        int r = rng_int(100);

        if (r < 75) {
            /* Single swap in D12 (keep 0 fixed) */
            if (d12_cnt <= 1 || d12_ncnt == 0) continue;
            int ri = rng_int(d12_cnt);
            int rem = d12_mem[ri];
            if (rem == 0) continue;
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

            int old_d12[P], old_d12T[P];
            memcpy(old_d12, delta_12, sizeof(delta_12));
            memcpy(old_d12T, delta_12T, sizeof(delta_12T));

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
            } else {
                D12[rem1] = 1; D12[rem2] = 1;
                D12[add1] = 0; D12[add2] = 0;
                D12T[mod_neg[rem1]] = 1; D12T[mod_neg[rem2]] = 1;
                D12T[mod_neg[add1]] = 0; D12T[mod_neg[add2]] = 0;
                memcpy(delta_12, old_d12, sizeof(delta_12));
                memcpy(delta_12T, old_d12T, sizeof(delta_12T));
            }
        }

        T *= alpha;
    }

    if (current_cost == 0) {
        int v12 = verify_v1v2();
        if (v12 == 0) return 1;
    }
    return 0;
}

/* Canonicalize D12 under shift + negation */
static void canonicalize(int *d12_ind, int *canon) {
    int elems[D12_SIZE], ne = 0;
    for (int i = 0; i < P; i++) if (d12_ind[i]) elems[ne++] = i;

    int best[D12_SIZE];
    for (int i = 0; i < D12_SIZE; i++) best[i] = P;

    for (int negate = 0; negate < 2; negate++) {
        for (int shift = 0; shift < P; shift++) {
            int transformed[D12_SIZE];
            for (int i = 0; i < ne; i++) {
                int x = elems[i];
                if (negate) x = (P - x) % P;
                transformed[i] = (x + shift) % P;
            }
            /* Sort */
            for (int i = 0; i < ne - 1; i++)
                for (int j = i + 1; j < ne; j++)
                    if (transformed[j] < transformed[i]) {
                        int t = transformed[i];
                        transformed[i] = transformed[j];
                        transformed[j] = t;
                    }
            /* Compare lexicographically */
            int smaller = 0;
            for (int i = 0; i < ne; i++) {
                if (transformed[i] < best[i]) { smaller = 1; break; }
                if (transformed[i] > best[i]) break;
            }
            if (smaller) memcpy(best, transformed, ne * sizeof(int));
        }
    }
    memcpy(canon, best, ne * sizeof(int));
}

/* Hash set for distinct D12 orbits */
static uint64_t hash_canon(int *canon, int len) {
    uint64_t h = 0;
    for (int i = 0; i < len; i++)
        h = h * 131 + canon[i];
    return h;
}

#define HASH_SIZE 65536
typedef struct {
    int canon[D12_SIZE];
    int used;
} HashEntry;
static HashEntry hash_table[HASH_SIZE];
static int n_distinct = 0;

static int hash_insert(int *canon) {
    uint64_t h = hash_canon(canon, D12_SIZE) % HASH_SIZE;
    while (hash_table[h].used) {
        if (memcmp(hash_table[h].canon, canon, D12_SIZE * sizeof(int)) == 0)
            return 0;
        h = (h + 1) % HASH_SIZE;
    }
    memcpy(hash_table[h].canon, canon, D12_SIZE * sizeof(int));
    hash_table[h].used = 1;
    n_distinct++;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2 + D11_SIZE) {
        fprintf(stderr, "Usage: %s <n_trials> <d11_elem1> ... <d11_elem%d>\n",
                argv[0], D11_SIZE);
        return 1;
    }

    int n_trials = atoi(argv[1]);
    int d11_list[D11_SIZE];
    for (int i = 0; i < D11_SIZE; i++)
        d11_list[i] = atoi(argv[i + 2]);

    init_mod_tables();

    /* Set up D11 indicator */
    memset(D11, 0, sizeof(D11));
    for (int i = 0; i < D11_SIZE; i++) D11[d11_list[i]] = 1;
    compute_delta_full(D11, delta_11);

    memset(hash_table, 0, sizeof(hash_table));

    int max_iter = 1500000;
    int found = 0;

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    for (int trial = 0; trial < n_trials; trial++) {
        int seed = trial * 10000 + 42;
        if (sa_search_d12(seed, max_iter)) {
            found++;
            int canon[D12_SIZE];
            canonicalize(D12, canon);
            hash_insert(canon);
        }

        if ((trial + 1) % 50 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts_now);
            double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                           (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
            fprintf(stderr, "[%d/%d] found=%d distinct=%d (%.0fs)\n",
                    trial + 1, n_trials, found, n_distinct, elapsed);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                   (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;

    printf("Trials: %d, Found valid: %d (%.1f%%)\n",
           n_trials, found, 100.0 * found / n_trials);
    printf("Distinct D12-orbits (omega): %d\n", n_distinct);
    printf("N_estimate = omega * %d = %d\n", 2 * P, n_distinct * 2 * P);
    printf("Time: %.1fs\n", elapsed);

    return 0;
}
