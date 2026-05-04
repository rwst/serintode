#ifndef SERINTODE_H
#define SERINTODE_H

typedef struct {
    const char *name;
    const char *description;
    int (*run)(int argc, char *argv[]);
} mode_descriptor;

extern const mode_descriptor mode_linear;
extern const mode_descriptor mode_nonlin;
extern const mode_descriptor mode_mahler;
extern const mode_descriptor mode_makelookup;

#endif
