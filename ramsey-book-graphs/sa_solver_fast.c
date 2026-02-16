/*
 * High-performance SA solver for R(B_{n-1}, B_n) = 4n-1.
 * OPTIMIZED VERSION with:
 *   - Fully incremental O(m) delta updates for ALL move types
 *   - Parallel tempering with multiple temperature chains
 *   - Longer iterations, periodic progress reporting
 *   - Incremental cost updates (delta-cost from changed differences only)
 *
 * Usage: ./sa_solver_fast <n> [max_seeds] [max_iter] [num_chains]
 *
 * Key improvement: The original solver used O(m^2) full recomputation for
 * D11 pair swaps and D12 double swaps. This version uses O(m) incremental
 * updates for ALL move types, giving ~m/2 speedup for those moves.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_M 300
#define MAX_CHAINS 16

/* Global parameters (shared across chains) */
static int N_PARAM, M, N_TOTAL;
static int RED_THRESH, BLUE_THRESH;
static int D12_SIZE;

/* Precomputed modular arithmetic tables */
static int mod_add[MAX_M][MAX_M];
static int mod_sub[MAX_M][MAX_M];
static int mod_neg[MAX_M];

/* Symmetric pairs (shared) */
static int pairs[MAX_M][2];
static int num_pairs;

/* ====== Per-chain state ====== */
typedef struct {
    /* Set indicators */
    int D11[MAX_M];
    int D12[MAX_M];
    int D12T[MAX_M];

    /* Autocorrelation arrays */
    int delta_11[MAX_M];
    int delta_12[MAX_M];
    int delta_12T[MAX_M];

    /* D11 pair membership */
    int pair_in_d11[MAX_M];

    /* D12 member/nonmember lists with index tracking */
    int d12_mem[MAX_M], d12_nmem[MAX_M];
    int d12_cnt, d12_ncnt;
    int d12_mem_idx[MAX_M], d12_nmem_idx[MAX_M];

    /* D11 pair in/out lists with index tracking */
    int d11_in[MAX_M], d11_out[MAX_M];
    int d11_in_cnt, d11_out_cnt;
    int d11_in_idx[MAX_M], d11_out_idx[MAX_M];

    /* Per-difference cost contribution (for incremental cost) */
    int cost_v11[MAX_M];  /* V1V1 cost contribution for difference d */
    int cost_v22[MAX_M];  /* V2V2 cost contribution for difference d */

    /* Current and best state */
    int current_cost;
    int best_cost;
    int best_D11[MAX_M];
    int best_D12[MAX_M];

    /* Temperature */
    double T;

    /* RNG state (xoshiro256**) */
    unsigned long long rng_s[4];

    /* Derived sizes */
    int d11_size;
} Chain;

/* ====== RNG ====== */
static inline unsigned long long rotl(const unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline unsigned long long rng_next(Chain *c) {
    const unsigned long long result = rotl(c->rng_s[1] * 5, 7) * 9;
    const unsigned long long t = c->rng_s[1] << 17;
    c->rng_s[2] ^= c->rng_s[0];
    c->rng_s[3] ^= c->rng_s[1];
    c->rng_s[1] ^= c->rng_s[2];
    c->rng_s[0] ^= c->rng_s[3];
    c->rng_s[2] ^= t;
    c->rng_s[3] = rotl(c->rng_s[3], 45);
    return result;
}

static void rng_seed(Chain *c, unsigned long long seed) {
    c->rng_s[0] = seed ^ 0x9E3779B97F4A7C15ULL;
    c->rng_s[1] = seed ^ 0x6C62272E07BB0142ULL;
    c->rng_s[2] = seed ^ 0xBF58476D1CE4E5B9ULL;
    c->rng_s[3] = seed ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) rng_next(c);
}

static inline int rng_int(Chain *c, int n) {
    return (int)(rng_next(c) % (unsigned long long)n);
}

static inline double rng_double(Chain *c) {
    return (rng_next(c) >> 11) * 0x1.0p-53;
}

/* ====== Modular arithmetic tables ====== */
static void init_mod_tables(void) {
    for (int a = 0; a < M; a++) {
        mod_neg[a] = (M - a) % M;
        for (int b = 0; b < M; b++) {
            mod_add[a][b] = (a + b) % M;
            mod_sub[a][b] = (a - b + M) % M;
        }
    }
}

/* ====== Build symmetric pairs ====== */
static void build_pairs(void) {
    num_pairs = 0;
    for (int x = 1; x < M; x++) {
        int neg_x = mod_neg[x];
        if (x <= neg_x) {
            pairs[num_pairs][0] = x;
            pairs[num_pairs][1] = neg_x;
            num_pairs++;
        }
    }
}

/* ====== Full delta computation (used only for initialization) ====== */
static void compute_delta_full(const int *ind, int *delta) {
    for (int d = 0; d < M; d++) {
        int count = 0;
        for (int a = 0; a < M; a++) {
            if (ind[a] && ind[mod_sub[a][d]])
                count++;
        }
        delta[d] = count;
    }
}

/* ====== Incremental delta update for single-element swap ======
 * When swapping rem -> add_el in set S (S already updated):
 * For d != 0:
 *   change = S[add-d] + S[add+d] - S[rem-d] - S[rem+d]
 *          + [rem-d == add] + [rem+d == add]
 */
static inline void update_delta_swap(int *delta, const int *S, int rem, int add_el) {
    for (int d = 1; d < M; d++) {
        int change = S[mod_sub[add_el][d]] + S[mod_add[add_el][d]]
                   - S[mod_sub[rem][d]] - S[mod_add[rem][d]];
        if (mod_sub[rem][d] == add_el) change++;
        if (mod_add[rem][d] == add_el) change++;
        delta[d] += change;
    }
}

/* ====== Cost computation ======
 * Compute per-difference cost contributions and total cost.
 * This is O(m) and only called during initialization and verification.
 */
static int compute_cost_full(Chain *c) {
    int d11_size = 0, d12_size = 0;
    for (int i = 0; i < M; i++) {
        d11_size += c->D11[i];
        d12_size += c->D12[i];
    }
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int cost = 0;
    for (int d = 1; d < M; d++) {
        int cv11 = c->delta_11[d] + c->delta_12[d];

        int c11 = 0;
        if (c->D11[d]) {
            int excess = cv11 - RED_THRESH;
            if (excess > 0) c11 = excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) c11 = excess;
        }
        c->cost_v11[d] = c11;

        int d22d = c->delta_11[d] + (M - 2 - 2 * d11_size) + 2 * c->D11[d];
        int cv22 = d22d + c->delta_12T[d];
        int c22 = 0;
        if (!c->D11[d]) {
            int excess = cv22 - RED_THRESH;
            if (excess > 0) c22 = excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d2 + cv22;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) c22 = excess;
        }
        c->cost_v22[d] = c22;

        cost += c11 + c22;
    }
    return cost;
}

/* ====== Incremental cost recomputation ======
 * After delta arrays change, recompute cost from the per-difference arrays.
 * Still O(m) but avoids redundant work when we know which differences changed.
 * For simplicity in this version, we just recompute the full cost in O(m).
 * The key savings come from O(m) delta updates instead of O(m^2).
 */
static int recompute_cost(Chain *c) {
    int d11_size = c->d11_size;
    int d12_size = c->d12_cnt;
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int cost = 0;
    for (int d = 1; d < M; d++) {
        int cv11 = c->delta_11[d] + c->delta_12[d];

        if (c->D11[d]) {
            int excess = cv11 - RED_THRESH;
            cost += (excess > 0) ? excess : 0;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            cost += (excess > 0) ? excess : 0;
        }

        int d22d = c->delta_11[d] + (M - 2 - 2 * d11_size) + 2 * c->D11[d];
        int cv22 = d22d + c->delta_12T[d];
        if (!c->D11[d]) {
            int excess = cv22 - RED_THRESH;
            cost += (excess > 0) ? excess : 0;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d2 + cv22;
            int excess = bc - BLUE_THRESH;
            cost += (excess > 0) ? excess : 0;
        }
    }
    return cost;
}

/* ====== V1V2 verification (only called when V1V1+V2V2 cost=0) ====== */
static int verify_v1v2(Chain *c) {
    int d11_size = 0, d12_size = 0;
    for (int i = 0; i < M; i++) {
        d11_size += c->D11[i];
        d12_size += c->D12[i];
    }
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int D22[MAX_M];
    for (int i = 0; i < M; i++)
        D22[i] = (i > 0 && !c->D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < M; d++) {
        int sigma = 0;
        for (int a = 0; a < M; a++) {
            if (c->D11[a] && c->D12[mod_sub[d][a]])
                sigma++;
        }
        int delta = 0;
        for (int a = 0; a < M; a++) {
            if (c->D12[a] && D22[mod_sub[a][d]])
                delta++;
        }
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

/* ====== List management ====== */
static void rebuild_d12_lists(Chain *c) {
    c->d12_cnt = 0;
    c->d12_ncnt = 0;
    memset(c->d12_mem_idx, -1, sizeof(int) * M);
    memset(c->d12_nmem_idx, -1, sizeof(int) * M);
    for (int i = 0; i < M; i++) {
        if (c->D12[i]) {
            c->d12_mem_idx[i] = c->d12_cnt;
            c->d12_mem[c->d12_cnt++] = i;
        } else {
            c->d12_nmem_idx[i] = c->d12_ncnt;
            c->d12_nmem[c->d12_ncnt++] = i;
        }
    }
}

static void rebuild_d11_lists(Chain *c) {
    c->d11_in_cnt = 0;
    c->d11_out_cnt = 0;
    for (int i = 0; i < num_pairs; i++) {
        if (c->pair_in_d11[i]) {
            c->d11_in_idx[i] = c->d11_in_cnt;
            c->d11_in[c->d11_in_cnt++] = i;
        } else {
            c->d11_out_idx[i] = c->d11_out_cnt;
            c->d11_out[c->d11_out_cnt++] = i;
        }
    }
}

static void d12_swap_to_mem(Chain *c, int elem) {
    int ni = c->d12_nmem_idx[elem];
    int last = c->d12_nmem[c->d12_ncnt - 1];
    c->d12_nmem[ni] = last;
    c->d12_nmem_idx[last] = ni;
    c->d12_ncnt--;
    c->d12_mem_idx[elem] = c->d12_cnt;
    c->d12_mem[c->d12_cnt++] = elem;
    c->d12_nmem_idx[elem] = -1;
}

static void d12_swap_to_nmem(Chain *c, int elem) {
    int mi = c->d12_mem_idx[elem];
    int last = c->d12_mem[c->d12_cnt - 1];
    c->d12_mem[mi] = last;
    c->d12_mem_idx[last] = mi;
    c->d12_cnt--;
    c->d12_nmem_idx[elem] = c->d12_ncnt;
    c->d12_nmem[c->d12_ncnt++] = elem;
    c->d12_mem_idx[elem] = -1;
}

static void d11_pair_swap_in(Chain *c, int pair_idx) {
    /* Move pair from out to in */
    int oi = c->d11_out_idx[pair_idx];
    int last = c->d11_out[c->d11_out_cnt - 1];
    c->d11_out[oi] = last;
    c->d11_out_idx[last] = oi;
    c->d11_out_cnt--;
    c->d11_in_idx[pair_idx] = c->d11_in_cnt;
    c->d11_in[c->d11_in_cnt++] = pair_idx;
    c->d11_out_idx[pair_idx] = -1;
}

static void d11_pair_swap_out(Chain *c, int pair_idx) {
    /* Move pair from in to out */
    int ii = c->d11_in_idx[pair_idx];
    int last = c->d11_in[c->d11_in_cnt - 1];
    c->d11_in[ii] = last;
    c->d11_in_idx[last] = ii;
    c->d11_in_cnt--;
    c->d11_out_idx[pair_idx] = c->d11_out_cnt;
    c->d11_out[c->d11_out_cnt++] = pair_idx;
    c->d11_in_idx[pair_idx] = -1;
}

/* ====== Incremental D11 pair swap ======
 * Remove pair (r0, r1), add pair (a0, a1).
 * Decompose into 4 single-element updates to delta_11:
 *   Step 1: remove r0 (set D11[r0]=0), add a0 (set D11[a0]=1) => update delta
 *   Step 2: remove r1 (set D11[r1]=0), add a1 (set D11[a1]=1) => update delta
 * But we need to handle the case where r0==a1 or r1==a0 (shouldn't happen for
 * distinct pairs, but let's be safe).
 *
 * Actually for a symmetric pair swap, we remove {r0,r1} and add {a0,a1}.
 * We can do this as two sequential swaps in the delta:
 *   First:  remove r0, add a placeholder that's already out... no, simpler:
 *
 *   We update delta_11 for the combined effect. The standard formula for
 *   a single swap (rem->add) works for one element at a time. We need to
 *   apply two swaps sequentially.
 *
 * Approach: temporarily set D11 to intermediate state after first sub-swap,
 * update delta, then do second sub-swap.
 */
static inline void update_delta_d11_pair(Chain *c, int r0, int r1, int a0, int a1) {
    /* State: D11 already has r0,r1 OUT and a0,a1 IN.
     * We need to get from old delta to new delta.
     *
     * Strategy: restore to old state, do two sequential incremental updates.
     */

    /* Temporarily restore old state */
    c->D11[r0] = 1; c->D11[r1] = 1;
    c->D11[a0] = 0; c->D11[a1] = 0;

    /* Sub-swap 1: remove r0, add a0 */
    c->D11[r0] = 0; c->D11[a0] = 1;
    update_delta_swap(c->delta_11, c->D11, r0, a0);

    /* Sub-swap 2: remove r1, add a1 */
    c->D11[r1] = 0; c->D11[a1] = 1;
    update_delta_swap(c->delta_11, c->D11, r1, a1);

    /* D11 is now in new state (r0,r1 out, a0,a1 in) - correct */
}

/* Revert D11 pair swap (undo: put r0,r1 back in, take a0,a1 back out) */
static inline void revert_delta_d11_pair(Chain *c, int r0, int r1, int a0, int a1) {
    /* State: D11 has r0,r1 OUT and a0,a1 IN (new state).
     * We want to revert to old state: r0,r1 IN, a0,a1 OUT.
     * Do the sub-swaps in reverse order.
     */

    /* Revert sub-swap 2: remove a1, add r1 */
    c->D11[a1] = 0; c->D11[r1] = 1;
    update_delta_swap(c->delta_11, c->D11, a1, r1);

    /* Revert sub-swap 1: remove a0, add r0 */
    c->D11[a0] = 0; c->D11[r0] = 1;
    update_delta_swap(c->delta_11, c->D11, a0, r0);

    /* D11 is now in old state (r0,r1 in, a0,a1 out) */
}

/* ====== Incremental D12 double swap ======
 * Remove (rem1, rem2), add (add1, add2). Same decomposition approach.
 */
static inline void update_delta_d12_double(Chain *c,
    int rem1, int rem2, int add1, int add2,
    int rem1t, int rem2t, int add1t, int add2t) {

    /* State: D12 already has rem1,rem2 OUT and add1,add2 IN.
     * Restore, then do two sequential swaps.
     */

    /* Restore old state */
    c->D12[rem1] = 1; c->D12[rem2] = 1;
    c->D12[add1] = 0; c->D12[add2] = 0;
    c->D12T[rem1t] = 1; c->D12T[rem2t] = 1;
    c->D12T[add1t] = 0; c->D12T[add2t] = 0;

    /* Sub-swap 1: remove rem1, add add1 */
    c->D12[rem1] = 0; c->D12[add1] = 1;
    c->D12T[rem1t] = 0; c->D12T[add1t] = 1;
    update_delta_swap(c->delta_12, c->D12, rem1, add1);
    update_delta_swap(c->delta_12T, c->D12T, rem1t, add1t);

    /* Sub-swap 2: remove rem2, add add2 */
    c->D12[rem2] = 0; c->D12[add2] = 1;
    c->D12T[rem2t] = 0; c->D12T[add2t] = 1;
    update_delta_swap(c->delta_12, c->D12, rem2, add2);
    update_delta_swap(c->delta_12T, c->D12T, rem2t, add2t);
}

static inline void revert_delta_d12_double(Chain *c,
    int rem1, int rem2, int add1, int add2,
    int rem1t, int rem2t, int add1t, int add2t) {

    /* Revert sub-swap 2: remove add2, add rem2 */
    c->D12[add2] = 0; c->D12[rem2] = 1;
    c->D12T[add2t] = 0; c->D12T[rem2t] = 1;
    update_delta_swap(c->delta_12, c->D12, add2, rem2);
    update_delta_swap(c->delta_12T, c->D12T, add2t, rem2t);

    /* Revert sub-swap 1: remove add1, add rem1 */
    c->D12[add1] = 0; c->D12[rem1] = 1;
    c->D12T[add1t] = 0; c->D12T[rem1t] = 1;
    update_delta_swap(c->delta_12, c->D12, add1, rem1);
    update_delta_swap(c->delta_12T, c->D12T, add1t, rem1t);
}

/* ====== Initialize a chain ====== */
static void init_chain(Chain *c, int d11_size, unsigned long long seed) {
    rng_seed(c, seed);
    c->d11_size = d11_size;

    /* Initialize D11 randomly (symmetric) */
    memset(c->D11, 0, sizeof(int) * M);
    int num_selected = d11_size / 2;

    int perm[MAX_M];
    for (int i = 0; i < num_pairs; i++) perm[i] = i;
    for (int i = num_pairs - 1; i > 0; i--) {
        int j = rng_int(c, i + 1);
        int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }
    for (int i = 0; i < num_selected; i++) {
        c->D11[pairs[perm[i]][0]] = 1;
        c->D11[pairs[perm[i]][1]] = 1;
    }
    memset(c->pair_in_d11, 0, sizeof(int) * num_pairs);
    for (int i = 0; i < num_selected; i++)
        c->pair_in_d11[perm[i]] = 1;
    rebuild_d11_lists(c);

    /* Initialize D12 randomly */
    memset(c->D12, 0, sizeof(int) * M);
    memset(c->D12T, 0, sizeof(int) * M);
    {
        int elems[MAX_M];
        for (int i = 0; i < M; i++) elems[i] = i;
        for (int i = M - 1; i > 0; i--) {
            int j = rng_int(c, i + 1);
            int tmp = elems[i]; elems[i] = elems[j]; elems[j] = tmp;
        }
        for (int i = 0; i < D12_SIZE; i++) {
            c->D12[elems[i]] = 1;
            c->D12T[mod_neg[elems[i]]] = 1;
        }
    }
    rebuild_d12_lists(c);

    /* Compute initial deltas from scratch */
    compute_delta_full(c->D11, c->delta_11);
    compute_delta_full(c->D12, c->delta_12);
    compute_delta_full(c->D12T, c->delta_12T);

    /* Compute initial cost */
    c->current_cost = compute_cost_full(c);
    c->best_cost = c->current_cost;
    memcpy(c->best_D11, c->D11, sizeof(int) * M);
    memcpy(c->best_D12, c->D12, sizeof(int) * M);
}

/* ====== Copy chain state (for parallel tempering swap) ====== */
static void copy_chain_state(Chain *dst, const Chain *src) {
    memcpy(dst->D11, src->D11, sizeof(int) * M);
    memcpy(dst->D12, src->D12, sizeof(int) * M);
    memcpy(dst->D12T, src->D12T, sizeof(int) * M);
    memcpy(dst->delta_11, src->delta_11, sizeof(int) * M);
    memcpy(dst->delta_12, src->delta_12, sizeof(int) * M);
    memcpy(dst->delta_12T, src->delta_12T, sizeof(int) * M);
    memcpy(dst->pair_in_d11, src->pair_in_d11, sizeof(int) * num_pairs);
    dst->d11_size = src->d11_size;
    dst->current_cost = src->current_cost;
    /* Rebuild lists since they contain pointers to internal arrays */
    rebuild_d12_lists(dst);
    rebuild_d11_lists(dst);
}

/* ====== Solve with parallel tempering ====== */
static Chain chains[MAX_CHAINS];
static int global_best_cost;
static int global_best_D11[MAX_M];
static int global_best_D12[MAX_M];

static int solve_pt(int d11_size, int seed, long long max_iter, int num_chains) {
    /* Temperature ladder for parallel tempering */
    double T_hot = 20.0;
    double T_cold = 0.00005;

    /* Initialize chains with geometric temperature spacing */
    for (int ch = 0; ch < num_chains; ch++) {
        unsigned long long ch_seed = (unsigned long long)seed * 1000003ULL + 42 + ch * 7919ULL;
        init_chain(&chains[ch], d11_size, ch_seed);

        if (num_chains > 1) {
            double frac = (double)ch / (double)(num_chains - 1);
            chains[ch].T = T_cold * pow(T_hot / T_cold, 1.0 - frac);
        } else {
            chains[ch].T = T_hot;
        }
    }

    global_best_cost = 999;
    for (int ch = 0; ch < num_chains; ch++) {
        if (chains[ch].best_cost < global_best_cost) {
            global_best_cost = chains[ch].best_cost;
            memcpy(global_best_D11, chains[ch].best_D11, sizeof(int) * M);
            memcpy(global_best_D12, chains[ch].best_D12, sizeof(int) * M);
        }
    }

    /* Cooling schedule for the coldest chain */
    /* For single chain: geometric cooling from T_hot to T_cold */
    double alpha;
    if (num_chains == 1) {
        alpha = exp(log(T_cold / T_hot) / (double)max_iter);
    } else {
        alpha = 1.0; /* PT chains keep fixed temperatures */
    }

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    int pt_swap_interval = 1000;  /* Try PT swaps every N iterations */
    int pt_swaps_accepted = 0, pt_swaps_attempted = 0;

    /* For single-chain mode, use adaptive reheating */
    long long last_imp_iter = 0;
    int reheat_count = 0;

    for (long long it = 0; it < max_iter; it++) {
        /* Run one SA step on each chain */
        for (int ch = 0; ch < num_chains; ch++) {
            Chain *c = &chains[ch];
            int r = rng_int(c, 100);

            if (r < 22) {
                /* ====== D11 pair swap (INCREMENTAL) ====== */
                if (c->d11_in_cnt == 0 || c->d11_out_cnt == 0) continue;

                int ri = rng_int(c, c->d11_in_cnt);
                int oi = rng_int(c, c->d11_out_cnt);
                int rp = c->d11_in[ri];
                int ap = c->d11_out[oi];

                int r0 = pairs[rp][0], r1 = pairs[rp][1];
                int a0 = pairs[ap][0], a1 = pairs[ap][1];

                /* Apply move and incremental delta update */
                c->D11[r0] = 0; c->D11[r1] = 0;
                c->D11[a0] = 1; c->D11[a1] = 1;
                update_delta_d11_pair(c, r0, r1, a0, a1);

                int new_cost = recompute_cost(c);
                int dc = new_cost - c->current_cost;

                if (dc <= 0 || rng_double(c) < exp(-(double)dc / fmax(c->T, 1e-10))) {
                    c->current_cost = new_cost;
                    c->pair_in_d11[rp] = 0;
                    c->pair_in_d11[ap] = 1;
                    d11_pair_swap_out(c, rp);
                    d11_pair_swap_in(c, ap);
                } else {
                    /* Revert */
                    revert_delta_d11_pair(c, r0, r1, a0, a1);
                }

            } else if (r < 82) {
                /* ====== D12 single swap (INCREMENTAL) ====== */
                if (c->d12_cnt == 0 || c->d12_ncnt == 0) continue;

                int ri = rng_int(c, c->d12_cnt);
                int ai = rng_int(c, c->d12_ncnt);
                int rem = c->d12_mem[ri];
                int add_el = c->d12_nmem[ai];

                int rem_t = mod_neg[rem];
                int add_t = mod_neg[add_el];

                c->D12[rem] = 0; c->D12[add_el] = 1;
                c->D12T[rem_t] = 0; c->D12T[add_t] = 1;

                update_delta_swap(c->delta_12, c->D12, rem, add_el);
                update_delta_swap(c->delta_12T, c->D12T, rem_t, add_t);

                int new_cost = recompute_cost(c);
                int dc = new_cost - c->current_cost;

                if (dc <= 0 || rng_double(c) < exp(-(double)dc / fmax(c->T, 1e-10))) {
                    c->current_cost = new_cost;
                    d12_swap_to_nmem(c, rem);
                    d12_swap_to_mem(c, add_el);
                } else {
                    /* Revert */
                    c->D12[rem] = 1; c->D12[add_el] = 0;
                    c->D12T[rem_t] = 1; c->D12T[add_t] = 0;
                    update_delta_swap(c->delta_12, c->D12, add_el, rem);
                    update_delta_swap(c->delta_12T, c->D12T, add_t, rem_t);
                }

            } else {
                /* ====== D12 double swap (INCREMENTAL) ====== */
                if (c->d12_cnt < 2 || c->d12_ncnt < 2) continue;

                int ri1 = rng_int(c, c->d12_cnt);
                int ri2 = rng_int(c, c->d12_cnt - 1);
                if (ri2 >= ri1) ri2++;
                int ai1 = rng_int(c, c->d12_ncnt);
                int ai2 = rng_int(c, c->d12_ncnt - 1);
                if (ai2 >= ai1) ai2++;

                int rem1 = c->d12_mem[ri1], rem2 = c->d12_mem[ri2];
                int add1 = c->d12_nmem[ai1], add2 = c->d12_nmem[ai2];
                int rem1t = mod_neg[rem1], rem2t = mod_neg[rem2];
                int add1t = mod_neg[add1], add2t = mod_neg[add2];

                /* Apply move */
                c->D12[rem1] = 0; c->D12[rem2] = 0;
                c->D12[add1] = 1; c->D12[add2] = 1;
                c->D12T[rem1t] = 0; c->D12T[rem2t] = 0;
                c->D12T[add1t] = 1; c->D12T[add2t] = 1;

                update_delta_d12_double(c, rem1, rem2, add1, add2,
                                        rem1t, rem2t, add1t, add2t);

                int new_cost = recompute_cost(c);
                int dc = new_cost - c->current_cost;

                if (dc <= 0 || rng_double(c) < exp(-(double)dc / fmax(c->T, 1e-10))) {
                    c->current_cost = new_cost;
                    d12_swap_to_nmem(c, rem1);
                    d12_swap_to_nmem(c, rem2);
                    d12_swap_to_mem(c, add1);
                    d12_swap_to_mem(c, add2);
                } else {
                    /* Revert */
                    revert_delta_d12_double(c, rem1, rem2, add1, add2,
                                            rem1t, rem2t, add1t, add2t);
                }
            }

            /* Track best */
            if (c->current_cost < c->best_cost) {
                c->best_cost = c->current_cost;
                memcpy(c->best_D11, c->D11, sizeof(int) * M);
                memcpy(c->best_D12, c->D12, sizeof(int) * M);

                if (c->current_cost < global_best_cost) {
                    global_best_cost = c->current_cost;
                    memcpy(global_best_D11, c->D11, sizeof(int) * M);
                    memcpy(global_best_D12, c->D12, sizeof(int) * M);
                    last_imp_iter = it;
                }

                if (c->best_cost == 0) {
                    int v12_cost = verify_v1v2(c);
                    if (v12_cost == 0) {
                        clock_gettime(CLOCK_MONOTONIC, &ts_now);
                        double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                                       (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
                        printf("  SOLUTION at iter %lldM chain %d (%.1fs)\n",
                               it / 1000000, ch, elapsed);
                        return 0;
                    } else {
                        c->best_cost = v12_cost;
                        if (v12_cost < global_best_cost) {
                            global_best_cost = v12_cost;
                        }
                    }
                }
            }

            /* Cool (single chain mode only) */
            if (num_chains == 1) {
                c->T *= alpha;

                /* Adaptive reheating for single chain */
                if (it - last_imp_iter > max_iter / 5 && c->T < 0.01) {
                    c->T = T_hot * 0.5;
                    last_imp_iter = it;
                    reheat_count++;
                }
            }
        }

        /* ====== Parallel tempering replica exchange ====== */
        if (num_chains > 1 && it % pt_swap_interval == 0 && it > 0) {
            /* Try to swap adjacent chains */
            /* Use a random RNG for PT decisions (use chain 0's RNG) */
            int ch1 = rng_int(&chains[0], num_chains - 1);
            int ch2 = ch1 + 1;
            pt_swaps_attempted++;

            double beta1 = 1.0 / fmax(chains[ch1].T, 1e-15);
            double beta2 = 1.0 / fmax(chains[ch2].T, 1e-15);
            double delta_E = (double)(chains[ch1].current_cost - chains[ch2].current_cost);
            double delta_beta = beta1 - beta2;
            double log_accept = delta_E * delta_beta;

            if (log_accept >= 0 || log(rng_double(&chains[0])) < log_accept) {
                /* Swap the full states between chains (but keep temperatures).
                 * copy_chain_state copies D11, D12, D12T, deltas, pair_in_d11,
                 * d11_size, current_cost, and rebuilds lists. */
                Chain tmp;
                copy_chain_state(&tmp, &chains[ch1]);
                copy_chain_state(&chains[ch1], &chains[ch2]);
                copy_chain_state(&chains[ch2], &tmp);

                pt_swaps_accepted++;
            }
        }

        /* ====== Progress reporting ====== */
        if (it % 5000000 == 0 && it > 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts_now);
            double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                           (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;

            if (num_chains > 1) {
                printf("    %lldM: gbest=%d costs=[",
                       it / 1000000, global_best_cost);
                for (int ch = 0; ch < num_chains; ch++) {
                    printf("%d", chains[ch].current_cost);
                    if (ch < num_chains - 1) printf(",");
                }
                printf("] Ts=[");
                for (int ch = 0; ch < num_chains; ch++) {
                    printf("%.4f", chains[ch].T);
                    if (ch < num_chains - 1) printf(",");
                }
                printf("] PT=%d/%d (%.1f it/s)\n",
                       pt_swaps_accepted, pt_swaps_attempted,
                       (double)it * num_chains / elapsed);
            } else {
                printf("    %lldM: cost=%d best=%d T=%.6f rh=%d (%.1f it/s)\n",
                       it / 1000000, chains[0].current_cost, global_best_cost,
                       chains[0].T, reheat_count,
                       (double)it / elapsed);
            }
            fflush(stdout);
        }

        /* Periodic delta verification (every 50M iterations) to catch drift */
        if (it % 50000000 == 0 && it > 0) {
            for (int ch = 0; ch < num_chains; ch++) {
                Chain *c = &chains[ch];
                int check_delta[MAX_M];

                /* Verify delta_11 */
                compute_delta_full(c->D11, check_delta);
                for (int d = 0; d < M; d++) {
                    if (check_delta[d] != c->delta_11[d]) {
                        printf("  DRIFT: chain %d delta_11[%d] inc=%d full=%d at iter %lldM\n",
                               ch, d, c->delta_11[d], check_delta[d], it / 1000000);
                        /* Fix it */
                        memcpy(c->delta_11, check_delta, sizeof(int) * M);
                        break;
                    }
                }

                /* Verify delta_12 */
                compute_delta_full(c->D12, check_delta);
                for (int d = 0; d < M; d++) {
                    if (check_delta[d] != c->delta_12[d]) {
                        printf("  DRIFT: chain %d delta_12[%d] inc=%d full=%d at iter %lldM\n",
                               ch, d, c->delta_12[d], check_delta[d], it / 1000000);
                        memcpy(c->delta_12, check_delta, sizeof(int) * M);
                        break;
                    }
                }

                /* Verify delta_12T */
                compute_delta_full(c->D12T, check_delta);
                for (int d = 0; d < M; d++) {
                    if (check_delta[d] != c->delta_12T[d]) {
                        printf("  DRIFT: chain %d delta_12T[%d] inc=%d full=%d at iter %lldM\n",
                               ch, d, c->delta_12T[d], check_delta[d], it / 1000000);
                        memcpy(c->delta_12T, check_delta, sizeof(int) * M);
                        break;
                    }
                }

                /* Verify cost */
                int full_cost = compute_cost_full(c);
                if (full_cost != c->current_cost) {
                    printf("  DRIFT: chain %d cost inc=%d full=%d at iter %lldM\n",
                           ch, c->current_cost, full_cost, it / 1000000);
                    c->current_cost = full_cost;
                }
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                   (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
    printf("  seed=%d: gbest=%d (%.1fs, %d chains)\n",
           seed, global_best_cost, elapsed, num_chains);

    /* Dump diagnostic info for near-solutions */
    if (global_best_cost > 0 && global_best_cost <= 8) {
        /* Use the global best state */
        Chain *c = &chains[0]; /* temporary */
        memcpy(c->D11, global_best_D11, sizeof(int) * M);
        memcpy(c->D12, global_best_D12, sizeof(int) * M);
        memset(c->D12T, 0, sizeof(int) * M);
        for (int i = 0; i < M; i++)
            if (c->D12[i]) c->D12T[mod_neg[i]] = 1;
        compute_delta_full(c->D11, c->delta_11);
        compute_delta_full(c->D12, c->delta_12);
        compute_delta_full(c->D12T, c->delta_12T);

        int d11_sz = 0, d12_sz = 0;
        for (int i = 0; i < M; i++) { d11_sz += c->D11[i]; d12_sz += c->D12[i]; }
        int d22_sz = M - 1 - d11_sz;
        int d1 = d11_sz + d12_sz, d2 = d22_sz + d12_sz;

        printf("  VIOLATIONS (V1V1+V2V2):\n");
        for (int d = 1; d < M; d++) {
            int cv11 = c->delta_11[d] + c->delta_12[d];
            if (c->D11[d]) {
                int excess = cv11 - RED_THRESH;
                if (excess > 0) printf("    V1V1 red d=%d: cv=%d thresh=%d excess=%d\n",
                                       d, cv11, RED_THRESH, excess);
            } else {
                int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
                int excess = bc - BLUE_THRESH;
                if (excess > 0) printf("    V1V1 blue d=%d: bc=%d thresh=%d excess=%d\n",
                                       d, bc, BLUE_THRESH, excess);
            }
            int d22d = c->delta_11[d] + (M - 2 - 2 * d11_sz) + 2 * c->D11[d];
            int cv22 = d22d + c->delta_12T[d];
            if (!c->D11[d]) {
                int excess = cv22 - RED_THRESH;
                if (excess > 0) printf("    V2V2 red d=%d: cv=%d thresh=%d excess=%d\n",
                                       d, cv22, RED_THRESH, excess);
            } else {
                int bc = (N_TOTAL - 2) - 2 * d2 + cv22;
                int excess = bc - BLUE_THRESH;
                if (excess > 0) printf("    V2V2 blue d=%d: bc=%d thresh=%d excess=%d\n",
                                       d, bc, BLUE_THRESH, excess);
            }
        }

        int v12 = verify_v1v2(c);
        printf("  V1V2 cost on best state: %d\n", v12);

        /* Dump near-solution to file */
        if (global_best_cost <= 4) {
            char fname[256];
            snprintf(fname, sizeof(fname), "/tmp/near_n%d_seed%d_fast.json", N_PARAM, seed);
            FILE *fp = fopen(fname, "w");
            if (fp) {
                fprintf(fp, "{\"n\": %d, \"m\": %d, \"cost\": %d,\n", N_PARAM, M, global_best_cost);
                fprintf(fp, " \"D11\": [");
                for (int i = 0, first = 1; i < M; i++)
                    if (global_best_D11[i]) { if (!first) fprintf(fp, ", "); fprintf(fp, "%d", i); first = 0; }
                fprintf(fp, "],\n \"D12\": [");
                for (int i = 0, first = 1; i < M; i++)
                    if (global_best_D12[i]) { if (!first) fprintf(fp, ", "); fprintf(fp, "%d", i); first = 0; }
                fprintf(fp, "]}\n");
                fclose(fp);
                printf("  Dumped near-solution to %s\n", fname);
            }
        }
    }

    fflush(stdout);
    return global_best_cost;
}


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [max_seeds] [max_iter_millions] [num_chains]\n", argv[0]);
        fprintf(stderr, "  n: book graph parameter (p = 2n-1)\n");
        fprintf(stderr, "  max_seeds: number of random restarts (default 20)\n");
        fprintf(stderr, "  max_iter_millions: iterations per seed in millions (default 200)\n");
        fprintf(stderr, "  num_chains: parallel tempering chains (default 4)\n");
        return 1;
    }

    N_PARAM = atoi(argv[1]);
    int max_seeds = argc > 2 ? atoi(argv[2]) : 20;
    long long max_iter = argc > 3 ? (long long)atoi(argv[3]) * 1000000LL : 200000000LL;
    int num_chains = argc > 4 ? atoi(argv[4]) : 4;

    if (num_chains > MAX_CHAINS) num_chains = MAX_CHAINS;
    if (num_chains < 1) num_chains = 1;

    M = 2 * N_PARAM - 1;
    N_TOTAL = 2 * M;
    RED_THRESH = N_PARAM - 2;
    BLUE_THRESH = N_PARAM - 1;
    D12_SIZE = N_PARAM - 1;

    printf("=== SA Solver (Fast) ===\n");
    printf("n=%d, m=%d (p=%d), N=%d\n", N_PARAM, M, M, N_TOTAL);
    printf("red_thresh=%d, blue_thresh=%d, |D12|=%d\n", RED_THRESH, BLUE_THRESH, D12_SIZE);
    printf("Seeds: %d, Iter: %lldM, Chains: %d\n", max_seeds, max_iter / 1000000, num_chains);
    fflush(stdout);

    init_mod_tables();
    build_pairs();
    printf("Symmetric pairs: %d\n", num_pairs);
    fflush(stdout);

    /* Determine D11 sizes to try */
    int half = (M - 1) / 2;
    int d11_sizes[10];
    int num_sizes = 0;

    for (int offset = 0; offset <= 3; offset++) {
        int candidates[2] = {half - offset, half + offset};
        for (int ci = 0; ci < 2; ci++) {
            int s = candidates[ci];
            if (s > 0 && s < M - 1 && s % 2 == 0 && (s / 2) <= num_pairs) {
                int dup = 0;
                for (int j = 0; j < num_sizes; j++)
                    if (d11_sizes[j] == s) { dup = 1; break; }
                if (!dup && num_sizes < 10) d11_sizes[num_sizes++] = s;
            }
        }
    }

    /* Sort by distance from half */
    for (int i = 0; i < num_sizes; i++) {
        for (int j = i + 1; j < num_sizes; j++) {
            if (abs(d11_sizes[j] - half) < abs(d11_sizes[i] - half)) {
                int tmp = d11_sizes[i];
                d11_sizes[i] = d11_sizes[j];
                d11_sizes[j] = tmp;
            }
        }
    }

    printf("D11 sizes to try:");
    for (int i = 0; i < num_sizes; i++) printf(" %d", d11_sizes[i]);
    printf("\n\n");
    fflush(stdout);

    int overall_best = 999;

    for (int si = 0; si < num_sizes; si++) {
        int d11_size = d11_sizes[si];
        printf("--- |D11| = %d ---\n", d11_size);
        fflush(stdout);

        for (int seed = 0; seed < max_seeds; seed++) {
            int result = solve_pt(d11_size, seed, max_iter, num_chains);
            if (result < overall_best) overall_best = result;

            if (result == 0) {
                printf("\n*** SOLUTION FOUND! ***\n");
                printf("n=%d, m=%d\n", N_PARAM, M);
                printf("D11 = [");
                for (int i = 0, first = 1; i < M; i++)
                    if (global_best_D11[i]) { if (!first) printf(", "); printf("%d", i); first = 0; }
                printf("]\n");
                printf("D12 = [");
                for (int i = 0, first = 1; i < M; i++)
                    if (global_best_D12[i]) { if (!first) printf(", "); printf("%d", i); first = 0; }
                printf("]\n");
                printf("|D11|=%d, |D12|=%d\n", d11_size, D12_SIZE);
                fflush(stdout);
                return 0;
            }
        }
    }

    printf("\nOverall best = %d\n", overall_best);
    return 1;
}
