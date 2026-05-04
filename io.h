#ifndef SERINTODE_IO_H
#define SERINTODE_IO_H

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"

#define MAX_LINE_LENGTH 100000L

void *s_malloc(size_t n);
void *s_calloc(size_t n, size_t s);

long read_series(const char *finname, long max_coeffs, mpz_t **out_S);
void free_series(mpz_t *S, long count);

/* Parses "<index> <value>" tuples per line. Returns max_index+1; allocates
 * *out_S of that length with mpz_init on every entry (zero for unknown
 * indices), and *out_known of the same length (1 = known, 0 = unknown).
 * Refuses negative indices, indices >= max_coeffs, and duplicates. Comment
 * lines beginning with '#' are skipped. Caller frees with free_series and
 * free(*out_known). */
long read_sparse_series(const char *finname, long max_coeffs,
                        mpz_t **out_S, char **out_known);

typedef void (*monomial_printer)(FILE *eqs, FILE *out,
                                 long term_index, const void *ctx);

void print_ode(FILE *eqs, FILE *out,
               mpz_t *N, long numterms, long max_poly_order,
               long stride, long bestnulldim,
               const char *header_eqs,
               monomial_printer print_monomial,
               const void *mode_ctx);

#endif
