#include "serintode.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const mode_descriptor *all_modes[] = {
    &mode_linear,
    &mode_nonlin,
    &mode_mahler,
    &mode_makelookup,
    NULL,
};

static void usage(const char *progname)
{
    fprintf(stderr, "Usage: %s <mode> [args]\n\nModes:\n", progname);
    for (const mode_descriptor **m = all_modes; *m; m++) {
        fprintf(stderr, "  %-12s %s\n", (*m)->name, (*m)->description);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    for (const mode_descriptor **m = all_modes; *m; m++) {
        if (strcmp(argv[1], (*m)->name) == 0) {
            return (*m)->run(argc - 1, argv + 1);
        }
    }
    fprintf(stderr, "Unknown mode: %s\n\n", argv[1]);
    usage(argv[0]);
    return EXIT_FAILURE;
}
