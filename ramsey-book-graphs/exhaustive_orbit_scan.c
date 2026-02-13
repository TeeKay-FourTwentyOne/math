/*
 * Exhaustive orbit scan: enumerate ALL multiplicative orbits of symmetric D11,
 * run SA to search for valid D12 for each.
 *
 * Compile: cc -O3 -o exhaustive_orbit_scan exhaustive_orbit_scan.c -lpthread -lm
 * Usage:   ./exhaustive_orbit_scan [options]
 *   -p P        prime p ≡ 3 mod 4 (default 43)
 *   -t THREADS  worker threads (default: nproc - 1)
 *   -n TRIALS   SA trials per orbit (default 10)
 *   -i ITER     SA iterations per trial (default 800000)
 *   -c FILE     checkpoint file (default p<P>_checkpoint.jsonl)
 *   -o FILE     output file (default p<P>_results.json)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <stdint.h>
#include <errno.h>

#define MAX_M 64  /* max prime (limited by uint64_t bitmask) */

/* ── Runtime parameters (set in main) ────────────────────────────── */
static int P, N_PARAM, M, N_TOTAL;
static int RED_THRESH, BLUE_THRESH;
static int D11_SIZE, D12_SIZE;
static int NUM_PAIRS, HALF_SELECT, NUM_MULTS;
static int SA_TRIALS, SA_ITER;

/* ── Shared read-only data ───────────────────────────────────────── */
static int mod_add[MAX_M][MAX_M];
static int mod_sub[MAX_M][MAX_M];
static int mod_neg[MAX_M];

static int pairs[MAX_M][2];       /* pairs[i] = {d, P-d} for d=1..(P-1)/2 */
static int multipliers[MAX_M];    /* representatives of Z_P* / {±1} */

static int     n_orbits;
static uint64_t *orbit_bitmasks;          /* canonical bitmask per orbit */
static int    (*orbit_d11)[MAX_M];        /* indicator array per orbit */

/* ── Synchronization ─────────────────────────────────────────────── */
static volatile int interrupted = 0;
static int next_orbit;                    /* atomic work counter */
static pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
static FILE *ckpt_fp;

/* ── Progress counters (updated under mutex) ─────────────────────── */
static int orbits_done = 0;
static int orbits_working = 0;
static int orbits_resumed = 0;
static int *orbit_done;  /* boolean: already in checkpoint */

/* ── Per-orbit result ────────────────────────────────────────────── */
typedef struct {
    int orbit_id;
    int working;
    int best_cost;
    int trials_run;
    int d12[MAX_M];  /* indicator array, valid only if working */
} orbit_result_t;

/* ── Per-thread SA context ───────────────────────────────────────── */
typedef struct {
    int D11[MAX_M];
    int D12[MAX_M];
    int D12T[MAX_M];
    int delta_11[MAX_M];
    int delta_12[MAX_M];
    int delta_12T[MAX_M];
    int d12_mem[MAX_M], d12_nmem[MAX_M];
    int d12_cnt, d12_ncnt;
    int d12_mem_idx[MAX_M], d12_nmem_idx[MAX_M];
    unsigned long long rng_s[4];
} sa_ctx_t;

/* ════════════════════════════════════════════════════════════════════
 *  RNG (xoshiro256**)
 * ════════════════════════════════════════════════════════════════════ */
static inline unsigned long long rotl(const unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}
static inline unsigned long long rng_next(sa_ctx_t *c) {
    const unsigned long long result = rotl(c->rng_s[1] * 5, 7) * 9;
    const unsigned long long t = c->rng_s[1] << 17;
    c->rng_s[2] ^= c->rng_s[0]; c->rng_s[3] ^= c->rng_s[1];
    c->rng_s[1] ^= c->rng_s[2]; c->rng_s[0] ^= c->rng_s[3];
    c->rng_s[2] ^= t; c->rng_s[3] = rotl(c->rng_s[3], 45);
    return result;
}
static void rng_seed(sa_ctx_t *c, unsigned long long seed) {
    c->rng_s[0] = seed ^ 0x9E3779B97F4A7C15ULL;
    c->rng_s[1] = seed ^ 0x6C62272E07BB0142ULL;
    c->rng_s[2] = seed ^ 0xBF58476D1CE4E5B9ULL;
    c->rng_s[3] = seed ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) rng_next(c);
}
static inline int rng_int(sa_ctx_t *c, int n) {
    return (int)(rng_next(c) % (unsigned long long)n);
}
static inline double rng_double(sa_ctx_t *c) {
    return (rng_next(c) >> 11) * 0x1.0p-53;
}

/* ════════════════════════════════════════════════════════════════════
 *  Initialization
 * ════════════════════════════════════════════════════════════════════ */
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
    NUM_PAIRS = (M - 1) / 2;
    for (int i = 0; i < NUM_PAIRS; i++) {
        pairs[i][0] = i + 1;
        pairs[i][1] = M - (i + 1);
    }
}

static void compute_multipliers(void) {
    /* Find primitive root mod P */
    int g = 0;
    for (int a = 2; a < P; a++) {
        int ord = 1;
        long long x = a;
        while (x != 1) { x = (x * a) % P; ord++; }
        if (ord == P - 1) { g = a; break; }
    }
    if (g == 0) { fprintf(stderr, "No primitive root found for p=%d\n", P); exit(1); }

    /* Representatives of Z_P* / {±1}: g^0, g^1, ..., g^{(P-3)/2} */
    NUM_MULTS = (P - 1) / 2;
    multipliers[0] = 1;
    for (int i = 1; i < NUM_MULTS; i++)
        multipliers[i] = (int)(((long long)multipliers[i-1] * g) % P);

    /* Verify: g^{NUM_MULTS} ≡ -1 mod P */
    int check = (int)(((long long)multipliers[NUM_MULTS-1] * g) % P);
    if (check != P - 1) {
        fprintf(stderr, "Multiplier verification failed: g=%d, g^%d=%d (expected %d)\n",
                g, NUM_MULTS, check, P-1);
        exit(1);
    }
    printf("Primitive root: %d, %d multipliers verified\n", g, NUM_MULTS);
}

/* ════════════════════════════════════════════════════════════════════
 *  Orbit enumeration
 * ════════════════════════════════════════════════════════════════════ */

static uint64_t canonicalize(uint64_t bm) {
    uint64_t canon = bm;
    for (int mi = 1; mi < NUM_MULTS; mi++) {
        int g = multipliers[mi];
        uint64_t new_bm = 0;
        for (int d = 1; d < M; d++) {
            if (bm & (1ULL << d))
                new_bm |= (1ULL << ((g * d) % P));
        }
        if (new_bm < canon) canon = new_bm;
    }
    return canon;
}

static int cmp_u64(const void *a, const void *b) {
    uint64_t va = *(const uint64_t *)a, vb = *(const uint64_t *)b;
    return (va > vb) - (va < vb);
}

static void enumerate_orbits(void) {
    HALF_SELECT = D11_SIZE / 2;

    /* Compute total combinations C(NUM_PAIRS, HALF_SELECT) */
    long long total_combos = 1;
    for (int i = 0; i < HALF_SELECT; i++) {
        total_combos = total_combos * (NUM_PAIRS - i) / (i + 1);
    }
    printf("Enumerating %lld symmetric D11 (%d pairs choose %d)...\n",
           total_combos, NUM_PAIRS, HALF_SELECT);

    uint64_t *all_bm = malloc(total_combos * sizeof(uint64_t));
    if (!all_bm) { fprintf(stderr, "malloc failed for %lld bitmasks\n", total_combos); exit(1); }

    /* Generate all combinations and canonicalize */
    int combo[MAX_M];
    for (int i = 0; i < HALF_SELECT; i++) combo[i] = i;

    long long idx = 0;
    while (1) {
        /* Build bitmask for this combination */
        uint64_t bm = 0;
        for (int i = 0; i < HALF_SELECT; i++) {
            bm |= (1ULL << pairs[combo[i]][0]);
            bm |= (1ULL << pairs[combo[i]][1]);
        }
        all_bm[idx++] = canonicalize(bm);

        /* Next combination */
        int k = HALF_SELECT - 1;
        while (k >= 0 && combo[k] == NUM_PAIRS - HALF_SELECT + k) k--;
        if (k < 0) break;
        combo[k]++;
        for (int j = k + 1; j < HALF_SELECT; j++) combo[j] = combo[j-1] + 1;
    }

    if (idx != total_combos) {
        fprintf(stderr, "Bug: generated %lld combos, expected %lld\n", idx, total_combos);
        exit(1);
    }

    /* Sort and deduplicate */
    qsort(all_bm, total_combos, sizeof(uint64_t), cmp_u64);

    n_orbits = 0;
    for (long long i = 0; i < total_combos; i++) {
        if (i == 0 || all_bm[i] != all_bm[i-1])
            n_orbits++;
    }

    orbit_bitmasks = malloc(n_orbits * sizeof(uint64_t));
    orbit_d11 = malloc(n_orbits * sizeof(int[MAX_M]));
    orbit_done = calloc(n_orbits, sizeof(int));
    if (!orbit_bitmasks || !orbit_d11 || !orbit_done) {
        fprintf(stderr, "malloc failed for orbit data\n"); exit(1);
    }

    int oidx = 0;
    for (long long i = 0; i < total_combos; i++) {
        if (i == 0 || all_bm[i] != all_bm[i-1]) {
            orbit_bitmasks[oidx] = all_bm[i];
            memset(orbit_d11[oidx], 0, sizeof(int) * MAX_M);
            for (int d = 0; d < M; d++)
                orbit_d11[oidx][d] = (all_bm[i] & (1ULL << d)) ? 1 : 0;
            oidx++;
        }
    }

    free(all_bm);

    int expected = (int)(total_combos / NUM_MULTS);
    printf("Found %d orbits (expected %d, ratio %.2f)\n",
           n_orbits, expected, (double)total_combos / n_orbits);
    if (n_orbits != expected) {
        printf("  WARNING: orbit count mismatch (some orbits may have smaller size)\n");
    }
}

/* ════════════════════════════════════════════════════════════════════
 *  SA core functions
 * ════════════════════════════════════════════════════════════════════ */

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

static int compute_cost(sa_ctx_t *c) {
    int d11_size = D11_SIZE;
    int d12_size = D12_SIZE;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;
    int cost = 0;

    for (int d = 1; d < M; d++) {
        /* V1V1 constraints */
        int cv11 = c->delta_11[d] + c->delta_12[d];
        if (c->D11[d]) {
            int excess = cv11 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }

        /* V2V2 constraints */
        int d22d = c->delta_11[d] + (M - 2 - 2 * d11_size) + 2 * c->D11[d];
        int cv22 = d22d + c->delta_12T[d];
        if (!c->D11[d]) {
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

static int verify_v1v2(sa_ctx_t *c) {
    int d11_size = D11_SIZE;
    int d12_size = D12_SIZE;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int D22[MAX_M];
    for (int i = 0; i < M; i++)
        D22[i] = (i > 0 && !c->D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < M; d++) {
        int sigma = 0;
        for (int a = 0; a < M; a++)
            if (c->D11[a] && c->D12[mod_sub[d][a]]) sigma++;
        int delta = 0;
        for (int a = 0; a < M; a++)
            if (c->D12[a] && D22[mod_sub[a][d]]) delta++;
        int common = sigma + delta;

        if (c->D12[d]) {
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

/* ════════════════════════════════════════════════════════════════════
 *  D12 member-list management
 * ════════════════════════════════════════════════════════════════════ */

static void rebuild_d12_lists(sa_ctx_t *c) {
    c->d12_cnt = 0; c->d12_ncnt = 0;
    for (int i = 0; i < M; i++) {
        if (c->D12[i]) {
            c->d12_mem_idx[i] = c->d12_cnt;
            c->d12_mem[c->d12_cnt++] = i;
            c->d12_nmem_idx[i] = -1;
        } else {
            c->d12_nmem_idx[i] = c->d12_ncnt;
            c->d12_nmem[c->d12_ncnt++] = i;
            c->d12_mem_idx[i] = -1;
        }
    }
}

static void d12_swap_to_mem(sa_ctx_t *c, int elem) {
    int ni = c->d12_nmem_idx[elem];
    int last = c->d12_nmem[c->d12_ncnt - 1];
    c->d12_nmem[ni] = last; c->d12_nmem_idx[last] = ni;
    c->d12_ncnt--;
    c->d12_mem_idx[elem] = c->d12_cnt;
    c->d12_mem[c->d12_cnt++] = elem;
    c->d12_nmem_idx[elem] = -1;
}

static void d12_swap_to_nmem(sa_ctx_t *c, int elem) {
    int mi = c->d12_mem_idx[elem];
    int last = c->d12_mem[c->d12_cnt - 1];
    c->d12_mem[mi] = last; c->d12_mem_idx[last] = mi;
    c->d12_cnt--;
    c->d12_nmem_idx[elem] = c->d12_ncnt;
    c->d12_nmem[c->d12_ncnt++] = elem;
    c->d12_mem_idx[elem] = -1;
}

/* ════════════════════════════════════════════════════════════════════
 *  SA search for D12 (fixed D11)
 * ════════════════════════════════════════════════════════════════════ */

static int sa_search_d12(sa_ctx_t *c, int trial_seed, int max_iter) {
    rng_seed(c, (unsigned long long)trial_seed * 999983ULL + 17);

    /* Random D12: always include 0, then choose D12_SIZE-1 from {1,...,M-1} */
    memset(c->D12, 0, sizeof(int) * MAX_M);
    memset(c->D12T, 0, sizeof(int) * MAX_M);
    c->D12[0] = 1;
    c->D12T[0] = 1;

    int elems[MAX_M];
    int ne = 0;
    for (int i = 1; i < M; i++) elems[ne++] = i;
    for (int i = ne - 1; i > 0; i--) {
        int j = rng_int(c, i + 1);
        int tmp = elems[i]; elems[i] = elems[j]; elems[j] = tmp;
    }
    for (int i = 0; i < D12_SIZE - 1; i++) {
        c->D12[elems[i]] = 1;
        c->D12T[mod_neg[elems[i]]] = 1;
    }
    rebuild_d12_lists(c);

    compute_delta_full(c->D12, c->delta_12);
    compute_delta_full(c->D12T, c->delta_12T);

    int current_cost = compute_cost(c);
    int local_best = current_cost;

    double T = 8.0;
    double T_min = 0.0001;
    double alpha = exp(log(T_min / T) / (double)max_iter);

    for (int it = 0; it < max_iter && !interrupted; it++) {
        if (current_cost == 0) {
            int v12 = verify_v1v2(c);
            if (v12 == 0) return 1;
        }

        int r = rng_int(c, 100);

        if (r < 75) {
            /* ── Single swap ── */
            if (c->d12_cnt <= 1 || c->d12_ncnt == 0) goto cool;
            int ri = rng_int(c, c->d12_cnt);
            int rem = c->d12_mem[ri];
            if (rem == 0) goto cool;
            int ai = rng_int(c, c->d12_ncnt);
            int add_el = c->d12_nmem[ai];

            c->D12[rem] = 0; c->D12[add_el] = 1;
            int rem_t = mod_neg[rem], add_t = mod_neg[add_el];
            c->D12T[rem_t] = 0; c->D12T[add_t] = 1;

            update_delta_inc(c->delta_12, c->D12, rem, add_el);
            update_delta_inc(c->delta_12T, c->D12T, rem_t, add_t);

            int new_cost = compute_cost(c);
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double(c) < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                d12_swap_to_nmem(c, rem);
                d12_swap_to_mem(c, add_el);
                if (current_cost < local_best) local_best = current_cost;
            } else {
                c->D12[rem] = 1; c->D12[add_el] = 0;
                c->D12T[rem_t] = 1; c->D12T[add_t] = 0;
                update_delta_inc(c->delta_12, c->D12, add_el, rem);
                update_delta_inc(c->delta_12T, c->D12T, add_t, rem_t);
            }
        } else {
            /* ── Double swap (incremental forward, memcpy restore) ── */
            if (c->d12_cnt <= 2 || c->d12_ncnt < 2) goto cool;

            int ri1 = rng_int(c, c->d12_cnt);
            int rem1 = c->d12_mem[ri1];
            if (rem1 == 0) goto cool;
            int ri2 = rng_int(c, c->d12_cnt - 1);
            int rem2 = c->d12_mem[ri2 >= ri1 ? ri2 + 1 : ri2];
            if (rem2 == 0) goto cool;

            int ai1 = rng_int(c, c->d12_ncnt);
            int ai2 = rng_int(c, c->d12_ncnt - 1);
            int add1 = c->d12_nmem[ai1];
            int add2 = c->d12_nmem[ai2 >= ai1 ? ai2 + 1 : ai2];

            /* Save deltas for restore on reject */
            int saved_d12[MAX_M], saved_d12T[MAX_M];
            memcpy(saved_d12, c->delta_12, sizeof(int) * M);
            memcpy(saved_d12T, c->delta_12T, sizeof(int) * M);

            /* Forward: two sequential incremental updates */
            int rem1_t = mod_neg[rem1], add1_t = mod_neg[add1];
            int rem2_t = mod_neg[rem2], add2_t = mod_neg[add2];

            c->D12[rem1] = 0; c->D12[add1] = 1;
            c->D12T[rem1_t] = 0; c->D12T[add1_t] = 1;
            update_delta_inc(c->delta_12, c->D12, rem1, add1);
            update_delta_inc(c->delta_12T, c->D12T, rem1_t, add1_t);

            c->D12[rem2] = 0; c->D12[add2] = 1;
            c->D12T[rem2_t] = 0; c->D12T[add2_t] = 1;
            update_delta_inc(c->delta_12, c->D12, rem2, add2);
            update_delta_inc(c->delta_12T, c->D12T, rem2_t, add2_t);

            int new_cost = compute_cost(c);
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double(c) < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                d12_swap_to_nmem(c, rem1); d12_swap_to_nmem(c, rem2);
                d12_swap_to_mem(c, add1); d12_swap_to_mem(c, add2);
                if (current_cost < local_best) local_best = current_cost;
            } else {
                /* Restore */
                c->D12[rem1] = 1; c->D12[rem2] = 1;
                c->D12[add1] = 0; c->D12[add2] = 0;
                c->D12T[rem1_t] = 1; c->D12T[rem2_t] = 1;
                c->D12T[add1_t] = 0; c->D12T[add2_t] = 0;
                memcpy(c->delta_12, saved_d12, sizeof(int) * M);
                memcpy(c->delta_12T, saved_d12T, sizeof(int) * M);
            }
        }
cool:
        T *= alpha;
    }

    if (current_cost == 0) {
        int v12 = verify_v1v2(c);
        if (v12 == 0) return 1;
    }
    return 0;
}

/* ════════════════════════════════════════════════════════════════════
 *  Process one orbit
 * ════════════════════════════════════════════════════════════════════ */

static orbit_result_t process_orbit(sa_ctx_t *c, int orbit_id) {
    orbit_result_t res;
    res.orbit_id = orbit_id;
    res.working = 0;
    res.best_cost = 9999;
    res.trials_run = 0;
    memset(res.d12, 0, sizeof(res.d12));

    /* Load D11 for this orbit */
    memcpy(c->D11, orbit_d11[orbit_id], sizeof(int) * MAX_M);
    compute_delta_full(c->D11, c->delta_11);

    for (int trial = 0; trial < SA_TRIALS && !interrupted; trial++) {
        int seed = orbit_id * 10000 + trial;
        int found = sa_search_d12(c, seed, SA_ITER);
        res.trials_run = trial + 1;

        if (found) {
            res.working = 1;
            res.best_cost = 0;
            memcpy(res.d12, c->D12, sizeof(int) * MAX_M);
            break;
        }

        /* Track best cost from this trial's final state */
        int cost = compute_cost(c);
        if (cost < res.best_cost) res.best_cost = cost;
    }

    return res;
}

/* ════════════════════════════════════════════════════════════════════
 *  Checkpoint I/O
 * ════════════════════════════════════════════════════════════════════ */

static void write_checkpoint_line(orbit_result_t *r) {
    char d11_buf[1024], d12_buf[1024];
    int pos;

    /* Format D11 */
    pos = 0;
    d11_buf[pos++] = '[';
    for (int i = 0, first = 1; i < M; i++) {
        if (orbit_d11[r->orbit_id][i]) {
            if (!first) d11_buf[pos++] = ',';
            pos += snprintf(d11_buf + pos, sizeof(d11_buf) - pos, "%d", i);
            first = 0;
        }
    }
    d11_buf[pos++] = ']'; d11_buf[pos] = '\0';

    /* Format D12 (if working) */
    if (r->working) {
        pos = 0;
        d12_buf[pos++] = '[';
        for (int i = 0, first = 1; i < M; i++) {
            if (r->d12[i]) {
                if (!first) d12_buf[pos++] = ',';
                pos += snprintf(d12_buf + pos, sizeof(d12_buf) - pos, "%d", i);
                first = 0;
            }
        }
        d12_buf[pos++] = ']'; d12_buf[pos] = '\0';
    } else {
        strcpy(d12_buf, "null");
    }

    fprintf(ckpt_fp,
        "{\"id\":%d,\"d11\":%s,\"w\":%d,\"cost\":%d,\"d12\":%s,\"t\":%d}\n",
        r->orbit_id, d11_buf, r->working, r->best_cost, d12_buf, r->trials_run);
    fflush(ckpt_fp);
}

static int load_checkpoint(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) return 0;

    int loaded = 0, working = 0;
    char line[4096];
    while (fgets(line, sizeof(line), f)) {
        int id = -1, w = 0;
        if (sscanf(line, "{\"id\":%d,", &id) == 1 && id >= 0 && id < n_orbits) {
            /* Parse working flag */
            char *wp = strstr(line, "\"w\":");
            if (wp) w = atoi(wp + 4);
            orbit_done[id] = 1;
            loaded++;
            if (w) working++;
        }
    }
    fclose(f);
    return loaded;
}

static void write_summary(const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) { fprintf(stderr, "Cannot write %s\n", filename); return; }

    double pw = (orbits_done > 0) ? (double)orbits_working / n_orbits : 0.0;

    fprintf(f, "{\n");
    fprintf(f, "  \"p\": %d,\n", P);
    fprintf(f, "  \"total_orbits\": %d,\n", n_orbits);
    fprintf(f, "  \"orbits_tested\": %d,\n", orbits_done);
    fprintf(f, "  \"working_orbits\": %d,\n", orbits_working);
    fprintf(f, "  \"p_working\": %.8f,\n", pw);
    fprintf(f, "  \"p_times_p_working\": %.4f,\n", P * pw);
    fprintf(f, "  \"sa_config\": {\n");
    fprintf(f, "    \"trials_per_orbit\": %d,\n", SA_TRIALS);
    fprintf(f, "    \"iterations_per_trial\": %d,\n", SA_ITER);
    fprintf(f, "    \"start_temp\": 8.0,\n");
    fprintf(f, "    \"cooling_rate_note\": \"geometric from 8.0 to 0.0001\"\n");
    fprintf(f, "  }\n");
    fprintf(f, "}\n");
    fclose(f);
}

/* ════════════════════════════════════════════════════════════════════
 *  Worker thread
 * ════════════════════════════════════════════════════════════════════ */

static void *worker(void *arg) {
    (void)arg;
    sa_ctx_t ctx;
    memset(&ctx, 0, sizeof(ctx));

    while (!interrupted) {
        int oid = __atomic_fetch_add(&next_orbit, 1, __ATOMIC_SEQ_CST);
        if (oid >= n_orbits) break;
        if (orbit_done[oid]) continue;

        orbit_result_t res = process_orbit(&ctx, oid);

        pthread_mutex_lock(&output_mutex);
        write_checkpoint_line(&res);
        orbits_done++;
        if (res.working) {
            orbits_working++;
            /* Print working orbit immediately */
            printf("  WORKING orbit %d: D11=[", res.orbit_id);
            for (int i = 0, first = 1; i < M; i++)
                if (orbit_d11[res.orbit_id][i]) {
                    if (!first) printf(",");
                    printf("%d", i);
                    first = 0;
                }
            printf("]\n");
        }
        if (orbits_done % 100 == 0 || res.working) {
            double pw = (double)orbits_working / n_orbits;
            printf("[%d/%d] done=%d working=%d p_working=%.5f p*pw=%.3f\n",
                   orbits_done, n_orbits - orbits_resumed,
                   orbits_done, orbits_working, pw, P * pw);
            fflush(stdout);
        }
        pthread_mutex_unlock(&output_mutex);
    }
    return NULL;
}

/* ════════════════════════════════════════════════════════════════════
 *  Signal handler
 * ════════════════════════════════════════════════════════════════════ */
static void sigint_handler(int sig) {
    (void)sig;
    interrupted = 1;
}

/* ════════════════════════════════════════════════════════════════════
 *  Main
 * ════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv) {
    P = 43;
    int num_threads = 0;
    SA_TRIALS = 10;
    SA_ITER = 800000;
    char ckpt_file[256] = "";
    char out_file[256] = "";

    /* Parse arguments */
    int opt;
    while ((opt = getopt(argc, argv, "p:t:n:i:c:o:")) != -1) {
        switch (opt) {
            case 'p': P = atoi(optarg); break;
            case 't': num_threads = atoi(optarg); break;
            case 'n': SA_TRIALS = atoi(optarg); break;
            case 'i': SA_ITER = atoi(optarg); break;
            case 'c': strncpy(ckpt_file, optarg, sizeof(ckpt_file)-1); break;
            case 'o': strncpy(out_file, optarg, sizeof(out_file)-1); break;
            default:
                fprintf(stderr,
                    "Usage: %s [-p P] [-t threads] [-n trials] [-i iter] [-c ckpt] [-o out]\n",
                    argv[0]);
                return 1;
        }
    }

    /* Validate */
    if (P < 7 || P >= MAX_M) {
        fprintf(stderr, "Prime must be in [7, %d)\n", MAX_M);
        return 1;
    }
    if (P % 4 != 3) {
        fprintf(stderr, "Prime must be ≡ 3 mod 4 (got %d ≡ %d mod 4)\n", P, P % 4);
        return 1;
    }

    /* Default filenames */
    if (ckpt_file[0] == '\0') snprintf(ckpt_file, sizeof(ckpt_file), "p%d_checkpoint.jsonl", P);
    if (out_file[0] == '\0') snprintf(out_file, sizeof(out_file), "p%d_results.json", P);

    /* Derive parameters */
    M = P;
    N_PARAM = (P + 1) / 2;
    N_TOTAL = 4 * N_PARAM - 2;
    RED_THRESH = N_PARAM - 2;
    BLUE_THRESH = N_PARAM - 1;
    D11_SIZE = (P + 1) / 2;
    D12_SIZE = (P - 1) / 2;

    /* Auto-detect threads */
    if (num_threads <= 0) {
        num_threads = (int)sysconf(_SC_NPROCESSORS_ONLN) - 1;
        if (num_threads < 1) num_threads = 1;
    }

    printf("═══════════════════════════════════════════════════════\n");
    printf(" Exhaustive Orbit Scan: p=%d\n", P);
    printf("═══════════════════════════════════════════════════════\n");
    printf("  n=%d, m=%d, |D11|=%d, |D12|=%d\n", N_PARAM, M, D11_SIZE, D12_SIZE);
    printf("  Thresholds: red=%d, blue=%d\n", RED_THRESH, BLUE_THRESH);
    printf("  SA: %d trials × %d iterations\n", SA_TRIALS, SA_ITER);
    printf("  Threads: %d\n", num_threads);
    printf("  Checkpoint: %s\n", ckpt_file);
    printf("  Output: %s\n", out_file);
    printf("═══════════════════════════════════════════════════════\n");
    fflush(stdout);

    /* Setup */
    init_mod_tables();
    build_pairs();
    compute_multipliers();

    struct timespec ts_start;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    enumerate_orbits();

    struct timespec ts_enum;
    clock_gettime(CLOCK_MONOTONIC, &ts_enum);
    double enum_time = (ts_enum.tv_sec - ts_start.tv_sec) +
                       (ts_enum.tv_nsec - ts_start.tv_nsec) * 1e-9;
    printf("Orbit enumeration: %.2fs\n", enum_time);

    /* Load checkpoint */
    int resumed = load_checkpoint(ckpt_file);
    if (resumed > 0) {
        /* Count working from checkpoint */
        /* Reparse to count working */
        FILE *f = fopen(ckpt_file, "r");
        char line[4096];
        int working_count = 0;
        if (f) {
            while (fgets(line, sizeof(line), f)) {
                char *wp = strstr(line, "\"w\":");
                if (wp && atoi(wp + 4) == 1) working_count++;
            }
            fclose(f);
        }
        orbits_done = resumed;
        orbits_working = working_count;
        orbits_resumed = resumed;
        printf("Resumed from checkpoint: %d orbits done (%d working), %d remaining\n",
               resumed, working_count, n_orbits - resumed);
    }

    /* Open checkpoint for append */
    ckpt_fp = fopen(ckpt_file, "a");
    if (!ckpt_fp) {
        fprintf(stderr, "Cannot open checkpoint file: %s\n", ckpt_file);
        return 1;
    }

    /* Estimate runtime */
    int remaining = n_orbits - resumed;
    double est_per_orbit = (double)SA_TRIALS * SA_ITER * M * 8e-9;  /* rough estimate */
    double est_total = est_per_orbit * remaining / num_threads;
    printf("\nEstimated time: %.0f seconds (%.1f hours) for %d remaining orbits\n",
           est_total, est_total / 3600, remaining);
    printf("Starting SA search...\n\n");
    fflush(stdout);

    /* Install signal handler */
    signal(SIGINT, sigint_handler);

    /* Launch workers */
    next_orbit = 0;
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], NULL, worker, NULL);
    }

    /* Wait for completion */
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    free(threads);

    fclose(ckpt_fp);

    /* Final timing */
    struct timespec ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double total_time = (ts_end.tv_sec - ts_start.tv_sec) +
                        (ts_end.tv_nsec - ts_start.tv_nsec) * 1e-9;

    double pw = (double)orbits_working / n_orbits;

    printf("\n═══════════════════════════════════════════════════════\n");
    if (interrupted) printf("  INTERRUPTED (partial results saved)\n");
    printf("  RESULTS: p=%d\n", P);
    printf("  Total orbits:      %d\n", n_orbits);
    printf("  Orbits tested:     %d\n", orbits_done);
    printf("  Working orbits:    %d\n", orbits_working);
    printf("  p_working:         %.6f\n", pw);
    printf("  p × p_working:     %.4f\n", P * pw);
    printf("  Time:              %.1fs (%.2f hours)\n", total_time, total_time / 3600);
    printf("  Checkpoint:        %s\n", ckpt_file);
    printf("═══════════════════════════════════════════════════════\n");

    /* Write summary */
    write_summary(out_file);
    printf("Summary written to %s\n", out_file);

    /* Cleanup */
    free(orbit_bitmasks);
    free(orbit_d11);
    free(orbit_done);

    return 0;
}
