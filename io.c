#include "io.h"

#include <errno.h>
#include <string.h>

void *s_malloc(size_t n)
{
    void *p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "out of memory (malloc %zu): %s\n", n, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return p;
}

void *s_calloc(size_t n, size_t s)
{
    void *p = calloc(n, s);
    if (p == NULL) {
        fprintf(stderr, "out of memory (calloc %zu*%zu): %s\n", n, s, strerror(errno));
        exit(EXIT_FAILURE);
    }
    return p;
}

long read_series(const char *finname, long max_coeffs, mpz_t **out_S)
{
    FILE *fin = fopen(finname, "r");
    if (fin == NULL) {
        fprintf(stderr, "Error: Could not open input file %s. %s\n",
                finname, strerror(errno));
        exit(EXIT_FAILURE);
    }

    mpz_t *S = (mpz_t *) s_malloc(max_coeffs * sizeof(mpz_t));
    char *input_string = (char *) s_malloc((MAX_LINE_LENGTH + 1L) * sizeof(char));

    mpz_t temp;
    mpz_init(temp);

    long num_coeffs = 0L;
    long nonzeroterms = 0L;
    long line_no = 0L;
    while (num_coeffs < max_coeffs) {
        char *fgcheck = fgets(input_string, MAX_LINE_LENGTH, fin);
        if (fgcheck == NULL) {
            if (feof(fin)) break;
            if (ferror(fin)) {
                fprintf(stderr, "Error while reading the file. %s\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
        line_no++;
        if (nonzeroterms == 0L) {
            if (mpz_set_str(temp, input_string, 10) != 0) {
                fprintf(stderr, "Error: malformed integer at line %ld of %s.\n",
                        line_no, finname);
                exit(EXIT_FAILURE);
            }
            if (mpz_cmp_ui(temp, 0) != 0) {
                nonzeroterms++;
                mpz_init_set(S[num_coeffs], temp);
                num_coeffs++;
            }
        } else {
            if (mpz_init_set_str(S[num_coeffs], input_string, 10) != 0) {
                fprintf(stderr, "Error: malformed integer at line %ld of %s.\n",
                        line_no, finname);
                exit(EXIT_FAILURE);
            }
            num_coeffs++;
        }
    }

    mpz_clear(temp);
    fclose(fin);
    free(input_string);
    *out_S = S;
    return num_coeffs;
}

void free_series(mpz_t *S, long count)
{
    for (long i = 0; i < count; i++) mpz_clear(S[i]);
    free(S);
}

void print_ode(FILE *eqs, FILE *out,
               mpz_t *N, long numterms, long max_poly_order,
               long stride, long bestnulldim,
               const char *header_eqs,
               monomial_printer print_monomial,
               const void *mode_ctx)
{
    mpz_t temp, temp2;
    mpz_inits(temp, temp2, NULL);

    fprintf(eqs, "%s", header_eqs);

    long firstorder = 0L;
    for (long i = 0L; i < numterms; i++) {
        mpz_set_ui(temp, 0L);
        for (long j = 0L; j < max_poly_order + 1L; j++) {
            mpz_abs(temp2, N[(i + j*numterms)*stride + bestnulldim]);
            mpz_add(temp, temp, temp2);
        }
        if (mpz_cmp_ui(temp, 0L) <= 0L) continue;

        firstorder++;
        if (firstorder > 1L) {
            fprintf(eqs, "+");
            fprintf(out, "+");
        }
        fprintf(eqs, "(");
        fprintf(out, "(");

        long firstterm = 0L;
        for (long k = 0L; k < max_poly_order + 1L; k++) {
            long idx = (i + k*numterms)*stride + bestnulldim;
            if (mpz_cmp_ui(N[idx], 0L) == 0L) continue;
            firstterm++;

            int sign_pos = mpz_cmp_ui(N[idx], 0L) > 0L;
            mpz_abs(temp, N[idx]);
            int abs_is_one = (mpz_cmp_ui(temp, 1L) == 0L);

            if (firstterm > 1L && sign_pos) {
                fprintf(eqs, "+");
                fprintf(out, "+");
                if (!abs_is_one) {
                    gmp_fprintf(eqs, "%Zd*", N[idx]);
                    gmp_fprintf(out, "%Zd*", N[idx]);
                }
            } else if (firstterm > 1L && !sign_pos) {
                if (!abs_is_one) {
                    gmp_fprintf(eqs, "%Zd*", N[idx]);
                    gmp_fprintf(out, "%Zd*", N[idx]);
                } else {
                    fprintf(eqs, "-");
                    fprintf(out, "-");
                }
            } else {
                if (k > 0) {
                    if (!abs_is_one) {
                        gmp_fprintf(eqs, "%Zd*", N[idx]);
                        gmp_fprintf(out, "%Zd*", N[idx]);
                    } else if (!sign_pos) {
                        fprintf(eqs, "-");
                        fprintf(out, "-");
                    }
                } else {
                    gmp_fprintf(eqs, "%Zd", N[idx]);
                    gmp_fprintf(out, "%Zd", N[idx]);
                }
            }

            if (k > 0) {
                if (k == 1) {
                    fprintf(eqs, "x");
                    fprintf(out, "x");
                } else {
                    fprintf(eqs, "x^%ld", k);
                    fprintf(out, "x^%ld", k);
                }
            }
        }
        fprintf(eqs, ")");
        fprintf(out, ")");

        print_monomial(eqs, out, i, mode_ctx);
    }

    fprintf(eqs, ":\n");
    fprintf(out, "=0\n");

    mpz_clears(temp, temp2, NULL);
}
