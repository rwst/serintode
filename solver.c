#include "solver.h"

#include <stdlib.h>
#include "iml.h"

void free_nullspace(mpz_t *N, long count)
{
    for (long i = 0; i < count; i++) mpz_clear(N[i]);
    free(N);
}

long compute_nullspace(long rows, long cols, mpz_t *M, mpz_t **out_N)
{
    long nulldim = kernelMP(rows, cols, M, out_N, 1L);

    if (nulldim <= 0L) {
        free(*out_N);
        *out_N = NULL;
        return 0L;
    }

    mpz_t sum, t;
    mpz_inits(sum, t, NULL);
    for (long i = 0L; i < cols; i++) {
        for (long j = 0L; j < nulldim; j++) {
            mpz_abs(t, (*out_N)[i*nulldim + j]);
            mpz_add(sum, sum, t);
        }
    }
    int spurious = (mpz_cmp_ui(sum, (unsigned long) nulldim) == 0);
    mpz_clears(sum, t, NULL);

    if (spurious) {
        free_nullspace(*out_N, cols * nulldim);
        *out_N = NULL;
        return 0L;
    }
    return nulldim;
}

linear_select_result select_best_linear(mpz_t *N,
    long ode_order, long max_poly_order, long nulldim,
    long max_ode_order_bound)
{
    linear_select_result r = { -1L, 0L, 0L, 0L, 0L };

    long min_max_poly_order = max_poly_order + 1L;
    long min_terms_used = 0L;
    long best = 0L;
    long final_orders = 0L;
    long last_max_found_ode_order = 0L;

    for (long k = 0L; k < nulldim; k++) {
        long ordersused = 0L;
        long termsused = 0L;
        long nonzeroterms = 0L;
        long max_found_poly_order = 0L;
        long max_found_ode_order = 0L;
        for (long i = 0L; i < ode_order + 1L; i++) {
            for (long j = 0L; j < max_poly_order + 1L; j++) {
                if (mpz_cmp_ui(N[(i + j*(ode_order + 1L))*nulldim + k], 0L) != 0L) {
                    nonzeroterms++;
                    if (max_found_poly_order < j) max_found_poly_order = j;
                }
            }
            if (termsused < nonzeroterms) {
                termsused = nonzeroterms;
                ordersused++;
                max_found_ode_order = i;
            }
        }
        if (max_found_ode_order == ode_order) {
            if (min_max_poly_order > max_found_poly_order) {
                min_max_poly_order = max_found_poly_order;
                min_terms_used = termsused;
                final_orders = ordersused;
                best = k;
            } else if (min_max_poly_order == max_found_poly_order) {
                if (min_terms_used > termsused) {
                    min_terms_used = termsused;
                    final_orders = ordersused;
                    best = k;
                }
            }
        }
        last_max_found_ode_order = max_found_ode_order;
    }

    if (final_orders < 2L) return r;
    if (min_terms_used < final_orders + 1L) return r;
    if (max_ode_order_bound < ode_order) return r;

    r.index = best;
    r.max_poly_order = min_max_poly_order;
    r.terms_used = min_terms_used;
    r.final_orders = final_orders;
    r.max_found_ode_order = last_max_found_ode_order;
    return r;
}

nonlin_select_result select_best_nonlin(mpz_t *N,
    long **orderexp,
    long ode_order, long numterms, long max_poly_order, long max_depth,
    long nulldim)
{
    nonlin_select_result r = { -1L, 0L, 0L, 0L, 0L, 0L };

    long min_terms_used = numterms + 1L;
    long min_max_poly_order = max_poly_order + 1L;
    long min_max_found_depth = max_depth + 1L;
    long best = 0L;
    long final_orders = numterms + 1L;
    long max_found_depth = 0L;
    long last_max_found_ode_order = 0L;

    for (long k = 0L; k < nulldim; k++) {
        long ordersused = 0L;
        long termsused = 0L;
        long nonzeroterms = 0L;
        long max_found_poly_order = 0L;
        long max_found_ode_order = 0L;
        for (long i = 0L; i < numterms; i++) {
            for (long j = 0L; j < max_poly_order + 1L; j++) {
                if (mpz_cmp_ui(N[(i + j*numterms)*nulldim + k], 0L) != 0L) {
                    nonzeroterms++;
                    if (max_found_poly_order < j) max_found_poly_order = j;
                }
            }
            if (termsused < nonzeroterms) {
                termsused = nonzeroterms;
                ordersused++;
                for (long l = 0L; l < ode_order + 1L; l++) {
                    if (orderexp[i][l] != 0L) {
                        if (max_found_ode_order < l) max_found_ode_order = l;
                        if (max_found_depth < orderexp[i][l]) max_found_depth = orderexp[i][l];
                    }
                }
            }
        }
        if (max_found_ode_order == ode_order) {
            if ((final_orders > ordersused) && (final_orders > 1L)) {
                final_orders = ordersused;
                min_max_found_depth = max_found_depth;
                min_max_poly_order = max_found_poly_order;
                min_terms_used = termsused;
                best = k;
            } else if ((final_orders == ordersused) && (final_orders > 1L)) {
                if (min_max_found_depth > max_found_depth) {
                    min_max_found_depth = max_found_depth;
                    min_max_poly_order = max_found_poly_order;
                    min_terms_used = termsused;
                    best = k;
                } else if (min_max_found_depth == max_found_depth) {
                    if (min_max_poly_order > max_found_poly_order) {
                        min_max_poly_order = max_found_poly_order;
                        min_terms_used = termsused;
                        best = k;
                    }
                }
            }
        }
        last_max_found_ode_order = max_found_ode_order;
    }

    if (final_orders < 2L) return r;
    if (last_max_found_ode_order < ode_order) return r;
    if (min_terms_used < final_orders + 1L && min_max_found_depth < 2L) return r;

    r.index = best;
    r.max_poly_order = min_max_poly_order;
    r.max_found_depth = min_max_found_depth;
    r.terms_used = min_terms_used;
    r.final_orders = final_orders;
    r.max_found_ode_order = last_max_found_ode_order;
    return r;
}
