# Pre-refactor Makefile (Stage 1 of plan.md).
# Builds the four current in-scope programs and runs regressions
# against tests/expected/ baselines.

CFLAGS       ?= -Wall -O2
IML_DIR      ?= /home/ralf/math/iml
CPPFLAGS      = -I$(IML_DIR)/include
LDFLAGS       = -L$(IML_DIR)/lib64
IML_LDLIBS    = -liml -lcblas -lgmp -lm
PLAIN_LDLIBS  = -lgmp -lm

PROGS = serintode_iml.o serintode_iml_nonlin.o serintode_iml_nonlin_lookup.o makelookup

.PHONY: all clean test
all: $(PROGS)

io.o: io.c io.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

solver.o: solver.c solver.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

serintode_iml.o: serintode_iml.c io.o solver.o io.h solver.h
	$(CC) $(CFLAGS) $(CPPFLAGS) serintode_iml.c io.o solver.o -o $@ $(LDFLAGS) $(IML_LDLIBS)

serintode_iml_nonlin.o: serintode_iml_nonlin.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(IML_LDLIBS)

serintode_iml_nonlin_lookup.o: serintode_iml_nonlin_lookup.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(IML_LDLIBS)

makelookup: makelookup.c
	$(CC) $(CFLAGS) $< -o $@ $(PLAIN_LDLIBS)

test: serintode_iml.o serintode_iml_nonlin.o
	@./serintode_iml.o tests/central_binomials.txt > /dev/null
	@diff -u tests/expected/central_binomials_linear.txt \
	         tests/central_binomials.txt_solution_6-checks.txt \
	    && echo "PASS: linear"
	@./serintode_iml_nonlin.o tests/central_binomials.txt > /dev/null
	@diff -u tests/expected/central_binomials_nonlin.txt \
	         tests/central_binomials.txt_nonlinsol_0-checks.txt \
	    && echo "PASS: nonlin"

clean:
	rm -f $(PROGS) io.o solver.o
