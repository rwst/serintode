#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "gmp.h"

#include "serintode.h"
#include "io.h"
#include "modes_nonlin.h"

static int makelookup_run(int argc, char *argv[])
{
    long const MIN_ODE_ORDER = 1L;
    long const MIN_DEPTH = 1L;

    setvbuf(stdout, NULL, _IONBF, 0);

    if (argc != 3) {
        fprintf(stderr, "Usage: serintode makelookup <max-ode-order> <num-coeffs>\n");
        return EXIT_FAILURE;
    }
    long MAX_ODE_ORDER = atol(argv[1]);
    long NUM_COEFFS = atol(argv[2]);
    if (MAX_ODE_ORDER <= 0L || NUM_COEFFS <= 0L) {
        fprintf(stderr, "Error: both arguments must be positive integers. Got MAX_ODE_ORDER=%ld, NUM_COEFFS=%ld.\n",
                MAX_ODE_ORDER, NUM_COEFFS);
        return EXIT_FAILURE;
    }
    long maxnumterms = NUM_COEFFS;

    const char *dirname = "lookuptables";
    int dcheck = mkdir(dirname, S_IRWXU | S_IRWXG | S_IRWXO);
    if (dcheck == -1 && errno != EEXIST) {
        fprintf(stderr, "ERROR: Could not create output directory. %s\n", strerror(errno));
        return EXIT_FAILURE;
    }

    mpz_t temp, temp2;
    mpz_inits(temp, temp2, NULL);

    for (long n = MIN_ODE_ORDER; n <= MAX_ODE_ORDER; n++) {
        long ODE_ORDER = n;
        long MAX_DEPTH_POSS = 1L;
        for (MAX_DEPTH_POSS = 1L; MAX_DEPTH_POSS <= maxnumterms; MAX_DEPTH_POSS++) {
            mpz_fac_ui(temp, ODE_ORDER + 1L + MAX_DEPTH_POSS);
            mpz_fac_ui(temp2, ODE_ORDER + 1L);
            mpz_divexact(temp, temp, temp2);
            mpz_fac_ui(temp2, MAX_DEPTH_POSS);
            mpz_divexact(temp, temp, temp2);
            if ((long) mpz_get_ui(temp) > maxnumterms) {
                MAX_DEPTH_POSS--;
                break;
            }
        }
        for (long p = MIN_DEPTH; p <= MAX_DEPTH_POSS; p++) {
            mpz_fac_ui(temp, ODE_ORDER + 1L + p);
            mpz_fac_ui(temp2, ODE_ORDER + 1L);
            mpz_divexact(temp, temp, temp2);
            mpz_fac_ui(temp2, p);
            mpz_divexact(temp, temp, temp2);
            long numterms = (long) mpz_get_ui(temp) - 1L;

            long *s = (long *) s_calloc(ODE_ORDER + 1L, sizeof(long));
            long **orderexp = (long **) s_malloc((numterms + 1L) * sizeof(long *));
            for (long i = 0L; i < numterms + 1L; i++) {
                orderexp[i] = (long *) s_calloc(ODE_ORDER + 1L, sizeof(long));
            }

            long arrindex = 0L;
            combs(orderexp, s, ODE_ORDER + 1L, p, 0L, &arrindex);

            char fname[4096];
            snprintf(fname, sizeof(fname), "%s/o%ldd%ld.txt", dirname, n, p);
            FILE *fout = fopen(fname, "w");
            if (fout == NULL) {
                fprintf(stderr, "Couldn't open the output file %s. %s\n", fname, strerror(errno));
                return EXIT_FAILURE;
            }
            for (long i = 0L; i < numterms; i++) {
                for (long j = 0L; j <= ODE_ORDER; j++) {
                    fprintf(fout, "%ld ", orderexp[i][j]);
                }
                fprintf(fout, "\n");
            }
            fclose(fout);

            for (long i = 0L; i < numterms + 1L; i++) free(orderexp[i]);
            free(orderexp);
            free(s);
        }
    }

    mpz_clears(temp, temp2, NULL);
    return EXIT_SUCCESS;
}

const mode_descriptor mode_makelookup = {
    .name = "makelookup",
    .description = "Generate term-exponent lookup tables for nonlin --lookup-dir",
    .run = makelookup_run,
};
