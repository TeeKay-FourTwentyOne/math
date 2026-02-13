/*
 * p43_exhaustive_scan.c
 *
 * Exhaustive scan of all 16,796 multiplicative orbits of symmetric D11
 * at p=43 for the R(B_{n-1}, B_n) = 4n-1 proof.
 *
 * For each orbit, runs simulated annealing to search for a valid D12.
 * Uses pthreads for parallelism and saves results incrementally as JSONL.
 *
 * Compile: cc -O3 -o p43_exhaustive_scan p43_exhaustive_scan.c -lpthread -lm
 * Usage:   ./p43_exhaustive_scan [options]
 *   --threads N    Worker threads (default: number of CPU cores)
 *   --trials N     SA trials per orbit (default: 10)
 *   --iter N       SA iterations per trial (default: 1500000)
 *   --resume       Skip orbits already in output file
 *   --output FILE  Output JSONL file (default: p43_exhaustive_results.jsonl)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <sys/time.h>

/* ---- Problem constants for p=43 ---- */
#define P          43
#define M          43
#define HALF_P     21       /* (p-1)/2 */
#define D11_SZ     22       /* (p+1)/2 */
#define D12_SZ     21       /* (p-1)/2 */
#define N_PARAM    22       /* n = (p+1)/2 */
#define N_TOTAL    86       /* 4n - 2 */
#define D1         43       /* d1 = D11_SZ + D12_SZ */
#define D2         41       /* d2 = D22_SZ + D12_SZ */
#define RED_THRESH 20       /* n - 2 */
#define BLUE_THRESH 21      /* n - 1 */
#define N_PAIRS    21       /* (p-1)/2 */
#define PAIRS_SEL  11       /* D11_SZ / 2 */

#define EXPECTED_COMBOS 352716
#define EXPECTED_ORBITS 16796
#define MAX_ORBITS      17000

/* ---- Modular arithmetic (shared, read-only after init) ---- */
static int mod_add_t[M][M];
static int mod_sub_t[M][M];
static int mod_neg_t[M];

/* ---- Orbit data (shared, read-only after enumeration) ---- */
static int orbit_d11_sorted[MAX_ORBITS][D11_SZ];
static int orbit_d11_ind[MAX_ORBITS][M];
static int num_orbits = 0;

/* ---- Results ---- */
typedef struct {
    _Atomic int done;
    int working;
    int best_cost;
    int trials_run;
    int d12_elems[D12_SZ];
} orbit_result_t;

static orbit_result_t results[MAX_ORBITS];

/* ---- Shared counters ---- */
static _Atomic int next_orbit_idx;
static _Atomic int completed_count;
static _Atomic int working_count;

/* ---- Output ---- */
static pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
static char output_path[512] = "p43_exhaustive_results.jsonl";
static FILE *output_fp = NULL;

/* ---- Timing ---- */
static struct timeval start_tv;

/* ---- Configuration ---- */
static int cfg_threads = 0;
static int cfg_trials = 10;
static int cfg_iter = 1500000;
static int cfg_resume = 0;

/* ================================================================
 * RNG: xoshiro256** (per-thread via struct)
 * ================================================================ */
typedef struct {
    unsigned long long s[4];
} rng_t;

static inline unsigned long long rotl64(unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline unsigned long long rng_next(rng_t *r) {
    unsigned long long result = rotl64(r->s[1] * 5, 7) * 9;
    unsigned long long t = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl64(r->s[3], 45);
    return result;
}

static void rng_seed(rng_t *r, unsigned long long seed) {
    r->s[0] = seed ^ 0x9E3779B97F4A7C15ULL;
    r->s[1] = seed ^ 0x6C62272E07BB0142ULL;
    r->s[2] = seed ^ 0xBF58476D1CE4E5B9ULL;
    r->s[3] = seed ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) rng_next(r);
}

static inline int rng_int(rng_t *r, int n) {
    return (int)(rng_next(r) % (unsigned long long)n);
}

static inline double rng_double(rng_t *r) {
    return (rng_next(r) >> 11) * 0x1.0p-53;
}

/* ================================================================
 * Initialization
 * ================================================================ */
static void init_mod_tables(void) {
    for (int a = 0; a < M; a++) {
        mod_neg_t[a] = (M - a) % M;
        for (int b = 0; b < M; b++) {
            mod_add_t[a][b] = (a + b) % M;
            mod_sub_t[a][b] = (a - b + M) % M;
        }
    }
}

/* ================================================================
 * Orbit enumeration
 * ================================================================ */
static void isort(int *a, int n) {
    for (int i = 1; i < n; i++) {
        int key = a[i], j = i - 1;
        while (j >= 0 && a[j] > key) { a[j+1] = a[j]; j--; }
        a[j+1] = key;
    }
}

static int cmp_arrays(const int *a, const int *b, int n) {
    for (int i = 0; i < n; i++) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

/* Check if sorted D11 is the lex-min representative of its orbit
 * under Z_43* / {±1} (21 multipliers). */
static int is_canonical(const int *d11_sorted) {
    int img[D11_SZ];
    for (int g = 2; g <= HALF_P; g++) {
        for (int i = 0; i < D11_SZ; i++)
            img[i] = (g * d11_sorted[i]) % M;
        isort(img, D11_SZ);
        if (cmp_arrays(img, d11_sorted, D11_SZ) < 0)
            return 0;
    }
    return 1;
}

static void enumerate_orbits(void) {
    int combo[PAIRS_SEL];
    for (int i = 0; i < PAIRS_SEL; i++) combo[i] = i;

    int total_combos = 0;
    num_orbits = 0;

    while (1) {
        total_combos++;

        /* Build sorted D11 from selected pairs */
        int d11[D11_SZ];
        int idx = 0;
        for (int i = 0; i < PAIRS_SEL; i++) {
            int x = combo[i] + 1;
            d11[idx++] = x;
            d11[idx++] = M - x;
        }
        isort(d11, D11_SZ);

        if (is_canonical(d11)) {
            if (num_orbits >= MAX_ORBITS) {
                fprintf(stderr, "ERROR: too many orbits\n");
                exit(1);
            }
            memcpy(orbit_d11_sorted[num_orbits], d11, sizeof(int) * D11_SZ);
            memset(orbit_d11_ind[num_orbits], 0, sizeof(int) * M);
            for (int i = 0; i < D11_SZ; i++)
                orbit_d11_ind[num_orbits][d11[i]] = 1;
            num_orbits++;
        }

        /* Next combination in lex order */
        int i = PAIRS_SEL - 1;
        while (i >= 0 && combo[i] == N_PAIRS - PAIRS_SEL + i) i--;
        if (i < 0) break;
        combo[i]++;
        for (int j = i + 1; j < PAIRS_SEL; j++)
            combo[j] = combo[j-1] + 1;
    }

    fprintf(stderr, "Enumerated %d combinations, %d canonical orbits\n",
            total_combos, num_orbits);
    if (total_combos != EXPECTED_COMBOS)
        fprintf(stderr, "WARNING: expected %d combos, got %d\n",
                EXPECTED_COMBOS, total_combos);
    if (num_orbits != EXPECTED_ORBITS)
        fprintf(stderr, "WARNING: expected %d orbits, got %d\n",
                EXPECTED_ORBITS, num_orbits);
}

/* ================================================================
 * Per-thread SA state
 * ================================================================ */
typedef struct {
    int D11[M];
    int D12[M];
    int D12T[M];
    int delta_11[M];
    int delta_12[M];
    int delta_12T[M];
    int d12_mem[M], d12_nmem[M];
    int d12_cnt, d12_ncnt;
    int d12_mem_idx[M], d12_nmem_idx[M];
    rng_t rng;
} sa_state_t;

/* ================================================================
 * Delta computation
 * ================================================================ */
static void compute_delta_full(const int *ind, int *delta) {
    for (int d = 0; d < M; d++) {
        int count = 0;
        for (int a = 0; a < M; a++)
            if (ind[a] && ind[mod_sub_t[a][d]]) count++;
        delta[d] = count;
    }
}

static inline void update_delta_inc(int *delta, const int *S, int rem, int add_el) {
    for (int d = 1; d < M; d++) {
        int change = S[mod_sub_t[add_el][d]] + S[mod_add_t[add_el][d]]
                   - S[mod_sub_t[rem][d]] - S[mod_add_t[rem][d]];
        if (mod_sub_t[rem][d] == add_el) change++;
        if (mod_add_t[rem][d] == add_el) change++;
        delta[d] += change;
    }
}

/* ================================================================
 * D12 member list management
 * ================================================================ */
static void rebuild_d12_lists(sa_state_t *st) {
    st->d12_cnt = 0; st->d12_ncnt = 0;
    for (int i = 0; i < M; i++) {
        if (st->D12[i]) {
            st->d12_mem_idx[i] = st->d12_cnt;
            st->d12_mem[st->d12_cnt++] = i;
            st->d12_nmem_idx[i] = -1;
        } else {
            st->d12_nmem_idx[i] = st->d12_ncnt;
            st->d12_nmem[st->d12_ncnt++] = i;
            st->d12_mem_idx[i] = -1;
        }
    }
}

static inline void list_to_mem(sa_state_t *st, int elem) {
    int ni = st->d12_nmem_idx[elem];
    int last = st->d12_nmem[st->d12_ncnt - 1];
    st->d12_nmem[ni] = last; st->d12_nmem_idx[last] = ni;
    st->d12_ncnt--;
    st->d12_mem_idx[elem] = st->d12_cnt;
    st->d12_mem[st->d12_cnt++] = elem;
    st->d12_nmem_idx[elem] = -1;
}

static inline void list_to_nmem(sa_state_t *st, int elem) {
    int mi = st->d12_mem_idx[elem];
    int last = st->d12_mem[st->d12_cnt - 1];
    st->d12_mem[mi] = last; st->d12_mem_idx[last] = mi;
    st->d12_cnt--;
    st->d12_nmem_idx[elem] = st->d12_ncnt;
    st->d12_nmem[st->d12_ncnt++] = elem;
    st->d12_mem_idx[elem] = -1;
}

/* ================================================================
 * Cost function (V1V1 + V2V2 constraints)
 * ================================================================ */
static int compute_cost(const sa_state_t *st) {
    int cost = 0;
    for (int d = 1; d < M; d++) {
        int cv11 = st->delta_11[d] + st->delta_12[d];
        if (st->D11[d]) {
            /* V1V1 red: A(d)+B(d) <= 20 */
            int excess = cv11 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            /* V1V1 blue: A(d)+B(d)-2 <= 21 → A(d)+B(d) <= 23 */
            int bc = (N_TOTAL - 2) - 2 * D1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }

        int d22d = st->delta_11[d] + (M - 2 - 2 * D11_SZ) + 2 * st->D11[d];
        int cv22 = d22d + st->delta_12T[d];
        if (!st->D11[d]) {
            /* V2V2 red: C(d)+B'(d) <= 20 */
            int excess = cv22 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            /* V2V2 blue: C(d)+B'(d)+2 <= 21 → C(d)+B'(d) <= 19 (tightest) */
            int bc = (N_TOTAL - 2) - 2 * D2 + cv22;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }
    return cost;
}

/* ================================================================
 * V1V2 verification (O(M^2), only called at cost=0)
 * ================================================================ */
static int verify_v1v2(const sa_state_t *st) {
    int D22[M];
    for (int i = 0; i < M; i++)
        D22[i] = (i > 0 && !st->D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < M; d++) {
        int sigma = 0;
        for (int a = 0; a < M; a++)
            if (st->D11[a] && st->D12[mod_sub_t[d][a]]) sigma++;
        int delta_val = 0;
        for (int a = 0; a < M; a++)
            if (st->D12[a] && D22[mod_sub_t[a][d]]) delta_val++;
        int common = sigma + delta_val;

        if (st->D12[d]) {
            int excess = common - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - D1 - D2 + common;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }
    return cost;
}

/* Full independent verification: recompute deltas from scratch */
static int verify_full(const sa_state_t *st) {
    int d12_check[M], d12t_check[M];
    compute_delta_full(st->D12, d12_check);
    compute_delta_full(st->D12T, d12t_check);

    for (int d = 1; d < M; d++) {
        if (d12_check[d] != st->delta_12[d] || d12t_check[d] != st->delta_12T[d]) {
            fprintf(stderr, "BUG: delta mismatch at d=%d (d12: %d vs %d, d12T: %d vs %d)\n",
                    d, st->delta_12[d], d12_check[d], st->delta_12T[d], d12t_check[d]);
            return -1;
        }
    }

    /* Recheck V1V1+V2V2 from scratch */
    for (int d = 1; d < M; d++) {
        int cv11 = st->delta_11[d] + d12_check[d];
        if (st->D11[d]) {
            if (cv11 > RED_THRESH) return -2;
        } else {
            int bc = (N_TOTAL - 2) - 2 * D1 + cv11;
            if (bc > BLUE_THRESH) return -3;
        }
        int d22d = st->delta_11[d] + (M - 2 - 2 * D11_SZ) + 2 * st->D11[d];
        int cv22 = d22d + d12t_check[d];
        if (!st->D11[d]) {
            if (cv22 > RED_THRESH) return -4;
        } else {
            int bc = (N_TOTAL - 2) - 2 * D2 + cv22;
            if (bc > BLUE_THRESH) return -5;
        }
    }

    return verify_v1v2(st) == 0 ? 0 : -6;
}

/* ================================================================
 * SA search for D12 given fixed D11
 * Returns 1 if valid D12 found, 0 otherwise.
 * ================================================================ */
static int sa_search_d12(sa_state_t *st, int trial_seed, int max_iter, int *best_cost_out) {
    rng_seed(&st->rng, (unsigned long long)trial_seed * 999983ULL + 17);

    /* Random D12: always include 0, pick D12_SZ-1 from {1,...,M-1} */
    memset(st->D12, 0, sizeof(int) * M);
    memset(st->D12T, 0, sizeof(int) * M);
    st->D12[0] = 1;
    st->D12T[0] = 1;

    int elems[M];
    int ne = 0;
    for (int i = 1; i < M; i++) elems[ne++] = i;
    for (int i = ne - 1; i > 0; i--) {
        int j = rng_int(&st->rng, i + 1);
        int tmp = elems[i]; elems[i] = elems[j]; elems[j] = tmp;
    }
    for (int i = 0; i < D12_SZ - 1; i++) {
        st->D12[elems[i]] = 1;
        st->D12T[mod_neg_t[elems[i]]] = 1;
    }
    rebuild_d12_lists(st);

    compute_delta_full(st->D12, st->delta_12);
    compute_delta_full(st->D12T, st->delta_12T);

    int current_cost = compute_cost(st);
    int best_cost = current_cost;

    double T = 8.0;
    double T_min = 0.0001;
    double alpha = exp(log(T_min / T) / (double)max_iter);

    for (int it = 0; it < max_iter; it++) {
        if (current_cost == 0) {
            if (verify_v1v2(st) == 0) {
                int vf = verify_full(st);
                if (vf == 0) {
                    *best_cost_out = 0;
                    return 1;
                }
                fprintf(stderr, "verify_full failed (%d), continuing\n", vf);
            }
        }

        int r = rng_int(&st->rng, 100);

        if (r < 80) {
            /* ---- Single swap ---- */
            if (st->d12_cnt <= 1 || st->d12_ncnt == 0) continue;
            int ri = rng_int(&st->rng, st->d12_cnt);
            int rem = st->d12_mem[ri];
            if (rem == 0) continue;
            int ai = rng_int(&st->rng, st->d12_ncnt);
            int add_el = st->d12_nmem[ai];

            st->D12[rem] = 0; st->D12[add_el] = 1;
            int rem_t = mod_neg_t[rem], add_t = mod_neg_t[add_el];
            st->D12T[rem_t] = 0; st->D12T[add_t] = 1;

            update_delta_inc(st->delta_12, st->D12, rem, add_el);
            update_delta_inc(st->delta_12T, st->D12T, rem_t, add_t);

            int new_cost = compute_cost(st);
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double(&st->rng) < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                list_to_nmem(st, rem);
                list_to_mem(st, add_el);
                if (current_cost < best_cost) best_cost = current_cost;
            } else {
                st->D12[rem] = 1; st->D12[add_el] = 0;
                st->D12T[rem_t] = 1; st->D12T[add_t] = 0;
                update_delta_inc(st->delta_12, st->D12, add_el, rem);
                update_delta_inc(st->delta_12T, st->D12T, add_t, rem_t);
            }
        } else {
            /* ---- Double swap (incremental O(M), not O(M^2)) ---- */
            if (st->d12_cnt <= 2 || st->d12_ncnt < 2) continue;

            int ri1 = rng_int(&st->rng, st->d12_cnt);
            int rem1 = st->d12_mem[ri1];
            if (rem1 == 0) continue;
            int ri2 = rng_int(&st->rng, st->d12_cnt - 1);
            int rem2 = st->d12_mem[ri2 >= ri1 ? ri2 + 1 : ri2];
            if (rem2 == 0) continue;

            int ai1 = rng_int(&st->rng, st->d12_ncnt);
            int ai2 = rng_int(&st->rng, st->d12_ncnt - 1);
            int add1 = st->d12_nmem[ai1];
            int add2 = st->d12_nmem[ai2 >= ai1 ? ai2 + 1 : ai2];

            /* Sub-swap 1: rem1 -> add1 */
            st->D12[rem1] = 0; st->D12[add1] = 1;
            int rem1_t = mod_neg_t[rem1], add1_t = mod_neg_t[add1];
            st->D12T[rem1_t] = 0; st->D12T[add1_t] = 1;
            update_delta_inc(st->delta_12, st->D12, rem1, add1);
            update_delta_inc(st->delta_12T, st->D12T, rem1_t, add1_t);

            /* Sub-swap 2: rem2 -> add2 */
            st->D12[rem2] = 0; st->D12[add2] = 1;
            int rem2_t = mod_neg_t[rem2], add2_t = mod_neg_t[add2];
            st->D12T[rem2_t] = 0; st->D12T[add2_t] = 1;
            update_delta_inc(st->delta_12, st->D12, rem2, add2);
            update_delta_inc(st->delta_12T, st->D12T, rem2_t, add2_t);

            int new_cost = compute_cost(st);
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double(&st->rng) < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                list_to_nmem(st, rem1); list_to_nmem(st, rem2);
                list_to_mem(st, add1); list_to_mem(st, add2);
                if (current_cost < best_cost) best_cost = current_cost;
            } else {
                /* Undo sub-swap 2 then sub-swap 1 */
                st->D12[add2] = 0; st->D12[rem2] = 1;
                st->D12T[add2_t] = 0; st->D12T[rem2_t] = 1;
                update_delta_inc(st->delta_12, st->D12, add2, rem2);
                update_delta_inc(st->delta_12T, st->D12T, add2_t, rem2_t);

                st->D12[add1] = 0; st->D12[rem1] = 1;
                st->D12T[add1_t] = 0; st->D12T[rem1_t] = 1;
                update_delta_inc(st->delta_12, st->D12, add1, rem1);
                update_delta_inc(st->delta_12T, st->D12T, add1_t, rem1_t);
            }
        }

        T *= alpha;
    }

    /* Final check */
    if (current_cost == 0 && verify_v1v2(st) == 0 && verify_full(st) == 0) {
        *best_cost_out = 0;
        return 1;
    }

    *best_cost_out = best_cost;
    return 0;
}

/* ================================================================
 * Output: write one orbit result as JSONL
 * ================================================================ */
static void write_orbit_result(int idx) {
    pthread_mutex_lock(&output_mutex);

    fprintf(output_fp, "{\"orbit_id\":%d,\"d11\":[", idx);
    for (int i = 0; i < D11_SZ; i++) {
        if (i) fputc(',', output_fp);
        fprintf(output_fp, "%d", orbit_d11_sorted[idx][i]);
    }
    fprintf(output_fp, "],\"working\":%s,\"best_cost\":%d,\"trials_run\":%d",
            results[idx].working ? "true" : "false",
            results[idx].best_cost, results[idx].trials_run);
    if (results[idx].working) {
        fprintf(output_fp, ",\"d12\":[");
        for (int i = 0; i < D12_SZ; i++) {
            if (i) fputc(',', output_fp);
            fprintf(output_fp, "%d", results[idx].d12_elems[i]);
        }
        fputc(']', output_fp);
    } else {
        fprintf(output_fp, ",\"d12\":null");
    }
    fprintf(output_fp, "}\n");
    fflush(output_fp);

    pthread_mutex_unlock(&output_mutex);
}

/* ================================================================
 * Resume: load existing JSONL results
 * ================================================================ */
static int load_resume_data(void) {
    FILE *f = fopen(output_path, "r");
    if (!f) return 0;

    int count = 0, w_count = 0;
    char line[8192];
    while (fgets(line, sizeof(line), f)) {
        char *p = strstr(line, "\"orbit_id\":");
        if (!p) continue;
        int orbit_id = atoi(p + 11);
        if (orbit_id < 0 || orbit_id >= MAX_ORBITS) continue;
        atomic_store(&results[orbit_id].done, 1);
        char *w = strstr(line, "\"working\":true");
        if (w) {
            results[orbit_id].working = 1;
            w_count++;
        }
        count++;
    }
    fclose(f);
    fprintf(stderr, "Resume: loaded %d completed orbits (%d working)\n", count, w_count);
    atomic_store(&completed_count, count);
    atomic_store(&working_count, w_count);
    return count;
}

/* ================================================================
 * Worker thread
 * ================================================================ */
static double elapsed_sec(void) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return (now.tv_sec - start_tv.tv_sec) + (now.tv_usec - start_tv.tv_usec) * 1e-6;
}

static void *worker_thread(void *arg) {
    (void)arg;
    sa_state_t st;

    while (1) {
        int idx = atomic_fetch_add(&next_orbit_idx, 1);
        if (idx >= num_orbits) break;

        if (atomic_load(&results[idx].done)) continue;

        /* Set up D11 for this orbit */
        memcpy(st.D11, orbit_d11_ind[idx], sizeof(int) * M);
        compute_delta_full(st.D11, st.delta_11);

        int found = 0;
        int overall_best = 999999;
        int trials_run = 0;

        for (int trial = 0; trial < cfg_trials; trial++) {
            trials_run++;
            int trial_best = 0;
            int seed = idx * 10000 + trial;
            if (sa_search_d12(&st, seed, cfg_iter, &trial_best)) {
                found = 1;
                overall_best = 0;
                break;
            }
            if (trial_best < overall_best) overall_best = trial_best;
        }

        results[idx].working = found;
        results[idx].best_cost = overall_best;
        results[idx].trials_run = trials_run;
        if (found) {
            int ei = 0;
            for (int i = 0; i < M && ei < D12_SZ; i++)
                if (st.D12[i]) results[idx].d12_elems[ei++] = i;
        }
        atomic_store(&results[idx].done, 1);

        int comp = atomic_fetch_add(&completed_count, 1) + 1;
        if (found) atomic_fetch_add(&working_count, 1);

        write_orbit_result(idx);

        /* Progress reporting */
        if (comp % 100 == 0 || found) {
            double t = elapsed_sec();
            int wk = atomic_load(&working_count);
            double pw = comp > 0 ? (double)wk / comp : 0;
            double eta = comp > 0 ? t / comp * (num_orbits - comp) : 0;
            fprintf(stderr, "[%d/%d] working=%d p_working=%.4f p*pw=%.3f "
                    "elapsed=%.0fs ETA=%.0fs%s\n",
                    comp, num_orbits, wk, pw, P * pw, t, eta,
                    found ? " <<FOUND>>" : "");
        }
    }

    return NULL;
}

/* ================================================================
 * Write final summary JSON
 * ================================================================ */
static void write_summary(void) {
    char summary_path[512];
    snprintf(summary_path, sizeof(summary_path), "p43_exhaustive_summary.json");

    int wk = atomic_load(&working_count);
    int comp = atomic_load(&completed_count);
    double pw = comp > 0 ? (double)wk / comp : 0;
    double t = elapsed_sec();

    FILE *sf = fopen(summary_path, "w");
    if (!sf) { fprintf(stderr, "Cannot open %s\n", summary_path); return; }

    fprintf(sf, "{\n");
    fprintf(sf, "  \"p\": %d,\n", P);
    fprintf(sf, "  \"total_orbits\": %d,\n", num_orbits);
    fprintf(sf, "  \"orbits_tested\": %d,\n", comp);
    fprintf(sf, "  \"working_orbits\": %d,\n", wk);
    fprintf(sf, "  \"p_working\": %.6f,\n", pw);
    fprintf(sf, "  \"p_times_p_working\": %.4f,\n", P * pw);
    fprintf(sf, "  \"sa_config\": {\n");
    fprintf(sf, "    \"trials_per_orbit\": %d,\n", cfg_trials);
    fprintf(sf, "    \"iterations_per_trial\": %d,\n", cfg_iter);
    fprintf(sf, "    \"start_temp\": 8.0,\n");
    fprintf(sf, "    \"cooling_rate_formula\": \"exp(log(0.0001/8.0)/iter)\"\n");
    fprintf(sf, "  },\n");
    fprintf(sf, "  \"elapsed_seconds\": %.1f,\n", t);
    fprintf(sf, "  \"threads\": %d\n", cfg_threads);
    fprintf(sf, "}\n");
    fclose(sf);

    /* Also print to stdout */
    printf("\n========================================\n");
    printf("RESULTS: p=%d exhaustive orbit scan\n", P);
    printf("  Total orbits:   %d\n", num_orbits);
    printf("  Orbits tested:  %d\n", comp);
    printf("  Working:        %d\n", wk);
    printf("  p_working:      %.6f\n", pw);
    printf("  p * p_working:  %.4f\n", P * pw);
    printf("  Time:           %.1fs\n", t);
    printf("  Threads:        %d\n", cfg_threads);
    printf("========================================\n");
    printf("Summary written to %s\n", summary_path);
    printf("Full results in %s\n", output_path);
}

/* ================================================================
 * Main
 * ================================================================ */
int main(int argc, char **argv) {
    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc)
            cfg_threads = atoi(argv[++i]);
        else if (strcmp(argv[i], "--trials") == 0 && i + 1 < argc)
            cfg_trials = atoi(argv[++i]);
        else if (strcmp(argv[i], "--iter") == 0 && i + 1 < argc)
            cfg_iter = atoi(argv[++i]);
        else if (strcmp(argv[i], "--resume") == 0)
            cfg_resume = 1;
        else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc)
            strncpy(output_path, argv[++i], sizeof(output_path) - 1);
        else if (strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s [--threads N] [--trials N] [--iter N] [--resume] [--output FILE]\n",
                   argv[0]);
            return 0;
        }
    }

    /* Auto-detect threads */
    if (cfg_threads <= 0) {
        cfg_threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
        if (cfg_threads <= 0) cfg_threads = 4;
        /* Leave 1 core free for system responsiveness */
        if (cfg_threads > 2) cfg_threads--;
    }

    fprintf(stderr, "p43 exhaustive orbit scan\n");
    fprintf(stderr, "  Threads:    %d\n", cfg_threads);
    fprintf(stderr, "  SA trials:  %d per orbit\n", cfg_trials);
    fprintf(stderr, "  SA iter:    %d per trial\n", cfg_iter);
    fprintf(stderr, "  Output:     %s\n", output_path);
    fprintf(stderr, "  Resume:     %s\n", cfg_resume ? "yes" : "no");

    /* Initialize */
    gettimeofday(&start_tv, NULL);
    init_mod_tables();

    fprintf(stderr, "\nPhase 1: Enumerating orbits...\n");
    enumerate_orbits();

    /* Initialize results */
    for (int i = 0; i < num_orbits; i++)
        atomic_store(&results[i].done, 0);
    atomic_store(&next_orbit_idx, 0);
    atomic_store(&completed_count, 0);
    atomic_store(&working_count, 0);

    /* Resume */
    if (cfg_resume) load_resume_data();

    /* Open output file (append mode for resume compatibility) */
    output_fp = fopen(output_path, "a");
    if (!output_fp) {
        fprintf(stderr, "ERROR: cannot open %s\n", output_path);
        return 1;
    }

    /* Phase 2: parallel SA */
    fprintf(stderr, "\nPhase 2: SA search (%d orbits, %d threads)...\n",
            num_orbits, cfg_threads);

    pthread_t *threads = malloc(sizeof(pthread_t) * cfg_threads);
    for (int i = 0; i < cfg_threads; i++)
        pthread_create(&threads[i], NULL, worker_thread, (void *)(long)i);
    for (int i = 0; i < cfg_threads; i++)
        pthread_join(threads[i], NULL);
    free(threads);

    fclose(output_fp);

    /* Summary */
    write_summary();

    return 0;
}
