# Stage 4 of plan.md: dispatcher binary `serintode` subsumes linear,
# nonlin (with optional --lookup-dir), and makelookup. Old standalone
# programs are gone.

CFLAGS       ?= -Wall -O2
IML_DIR      ?= /home/ralf/math/iml
CPPFLAGS      = -I$(IML_DIR)/include
LDFLAGS       = -L$(IML_DIR)/lib64
IML_LDLIBS    = -liml -lcblas -lgmp -lm

OBJS = io.o solver.o modes_linear.o modes_nonlin.o modes_mahler.o modes_makelookup.o

.PHONY: all clean test
all: serintode

io.o: io.c io.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

solver.o: solver.c solver.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

modes_linear.o: modes_linear.c serintode.h io.h solver.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

modes_nonlin.o: modes_nonlin.c serintode.h io.h solver.h modes_nonlin.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

modes_mahler.o: modes_mahler.c serintode.h io.h solver.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

modes_makelookup.o: modes_makelookup.c serintode.h io.h modes_nonlin.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

serintode: serintode.c serintode.h $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) serintode.c $(OBJS) -o $@ $(LDFLAGS) $(IML_LDLIBS)

test: serintode
	@./serintode linear --checks=6 tests/central_binomials.txt > /dev/null
	@diff -u tests/expected/central_binomials_linear.txt \
	         tests/central_binomials.txt_linear_6-checks.txt \
	    && echo "PASS: linear"
	@./serintode nonlin --checks=0 tests/central_binomials.txt > /dev/null
	@diff -u tests/expected/central_binomials_nonlin.txt \
	         tests/central_binomials.txt_nonlin_0-checks.txt \
	    && echo "PASS: nonlin"
	@./serintode mahler --checks=6 tests/thue_morse.txt > /dev/null
	@diff -u tests/expected/thue_morse_mahler.txt \
	         tests/thue_morse.txt_mahler_6-checks.txt \
	    && echo "PASS: mahler (Thue-Morse, k-automatic)"
	@./serintode mahler --checks=6 tests/stern.txt > /dev/null
	@diff -u tests/expected/stern_mahler.txt \
	         tests/stern.txt_mahler_6-checks.txt \
	    && echo "PASS: mahler (Stern, k-regular)"
	@./serintode mahler --checks=6 tests/fibonacci.txt > /dev/null
	@diff -u tests/expected/fibonacci_mahler.txt \
	         tests/fibonacci.txt_mahler_6-checks.txt \
	    && echo "PASS: mahler (Fibonacci, trivial rational)"

clean:
	rm -f serintode $(OBJS)
