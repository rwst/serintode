#ifndef SERINTODE_SOLVER_H
#define SERINTODE_SOLVER_H

#include "gmp.h"

/* Computes the right nullspace via IML's kernelMP, then applies the
 * spurious-zero check (sum of |entries| == nulldim => trivial-1 columns).
 * Returns the nullity if non-spurious; otherwise frees the nullspace
 * (setting *out_N to NULL) and returns 0. */
long compute_nullspace(long rows, long cols, mpz_t *M, mpz_t **out_N);

void free_nullspace(mpz_t *N, long count);

/* Returns the rank of the integer matrix M (rows x cols, row-major).
 * Computed as cols - nullity via IML's kernelMP. Destroys M's contents
 * (caller still owns the mpz_t storage and must mpz_clear). */
long compute_rank(long rows, long cols, mpz_t *M);

typedef struct {
    long index;
    long max_poly_order;
    long terms_used;
    long final_orders;
    long max_found_ode_order;
} linear_select_result;

linear_select_result select_best_linear(mpz_t *N,
    long ode_order, long max_poly_order, long nulldim,
    long max_ode_order_bound);

typedef struct {
    long index;
    long max_poly_order;
    long max_found_depth;
    long terms_used;
    long final_orders;
    long max_found_ode_order;
} nonlin_select_result;

nonlin_select_result select_best_nonlin(mpz_t *N,
    long **orderexp,
    long ode_order, long numterms, long max_poly_order, long max_depth,
    long nulldim);

#endif
