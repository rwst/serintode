INTRODUCTION
serintode is a C program for searching for ODEs which annihilate an integer series up to a certain order.

It ships as a single dispatcher binary `serintode` with five subcommands:
  serintode linear     <input> [--checks=N] [--min-order=N] [--max-coeffs=N]
  serintode nonlin     <input> [--checks=N] [--min-order=N] [--max-coeffs=N]
                               [--min-depth=N] [--max-depth=N] [--lookup-dir=PATH]
  serintode mahler     <input> [--checks=N] [--min-order=N] [--max-coeffs=N]
                               [--k-min=N] [--k-max=N]
  serintode kkernel    <input> [--k=N] [--max-depth=N] [--max-coeffs=N]
  serintode makelookup <max-ode-order> <num-coeffs>

`linear` searches for linear ODEs. `nonlin` searches for algebraic (nonlinear) ODEs up to the given depth. `mahler` searches for Mahler functional equations `sum_i p_i(x) f(x^(k^i)) = 0`, the structural detector for k-regular sequences (Allouche-Shallit). `kkernel` is a diagnostic: it computes the rank of the Q-span of the k-kernel `{(a_{k^i n + j})_n}` per depth — the rank stabilises iff the sequence is k-regular. `makelookup` precomputes the term-exponent tables that `nonlin --lookup-dir=lookuptables` reads (a speed optimisation; without the flag, the tables are enumerated in-process).

A separate FLINT-based prototype `serintode_flint.c` exists but is out of scope of the current refactor; it has no internal verification and can give spurious results.

The program automatically determines the number of coefficients in the input file and searches for ODEs of increasing order, polynomial-coefficient degree, and (for nonlin) nonlinearity depth, according to the number of coefficients and the number of checks required. It outputs the result both to stdout and to a file `<input>_<mode>_<NUM_CHECKS>-checks.txt` if a solution is found.

Input files should have one base-10 integer coefficient per line. Leading zero coefficients are stripped automatically.

If there are more than 10,000 coefficients, raise `--max-coeffs` accordingly.

If the coefficients are larger than 100,000 digits, raise `MAX_LINE_LENGTH` (a #define in `io.h`) and recompile.

The program looks for polynomial coefficients all of the same degree, the degree being determined by the formula:
MAX_POLY_ORDER=floor((NUM_COEFFS-NUM_CHECKS-ODE_ORDER)/(ODE_ORDER+1L))-1;
The maximum ODE order is found from the formula:
MAX_ODE_ORDER=floor((NUM_COEFFS-NUM_CHECKS)/2)-1;
since an ODE with only constant coefficients is not considered a solution, the degree of each polynomial coefficient should be at least 1. 

The `nonlin` mode searches for algebraic ODEs up to an order determined by the number of coefficients and up to n-products with the "depth" n determined also by the number of coefficients. An order 2 depth 3 ODE would have terms:
p1*y + p2*Dx + p3*Dx^2 + p4*y^3 + p5*y^2*Dx + p6*y^2*Dx^2 + p7*y*(Dx)^2 + p8*y*Dx*Dx^2 + p9*y*(Dx^2)^2 + p10*(Dx)^3 + p11*(Dx)^2*Dx^2 + p12*Dx*(Dx^2)^2 + p13*(Dx^2)^3
where p_k are polynomial coefficients.


INSTALLATION OF IML VERSION
To install the IML versions in Linux, ATLAS and GMP are needed. To install ATLAS:
apt-get -y install libatlas-base-dev

To install GMP, go to some directory, then:
apt-get -y install m4
curl https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2 -o gmp-6.1.2.tar.bz2
tar -xjvf gmp-6.1.2.tar.bz2
cd gmp-6.1.2
./configure
make 
make install
make check
make clean

To install IML, go to some directory, then:
curl http://www.cs.uwaterloo.ca/~astorjoh/iml-1.0.5.tar.bz2 -o ~/iml-1.0.5.tar.bz2 -L
tar -xjvf iml-1.0.5.tar.bz2
cd iml-1.0.5
./configure
make
make install
make clean


INSTALLATION OF FLINT VERSION
flint depends on GMP and MPFR. GMP can be installed as above. To install MPFR:
apt-get -y install libmpfr-dev

To install flint, go to some directory, then:
curl http://www.flintlib.org/flint-2.5.2.tar.gz -o flint-2.5.2.tar.gz
tar -zxvf flint-2.5.2.tar.gz
cd flint-2.5.2
./configure
make
make install
make clean


DIFFERENCES BETWEEN IML AND FLINT VERSIONS
Both the IML and flint versions rely on the construction of a matrix nullspace, the output solution being one of the null vectors (all are solutions). Due to the way each library constructs the nullspace, the output from each program can be different. Both outputs give a correct ODE. The flint solutions may have smaller polynomial coefficient orders, although I haven't tested the differences extensively enough to say that's always or even usually the case. The IML documentation states that rather than "only compute a basis for the rational nullspace", it can "produce an integer basis for the sublattice of all integer vectors in the right kernel of the input matrix". Perhaps the differences in the null space output between flint and IML are due to IML going beyond the rational nullspace and considering the integer sublattice of the nullspace. 

The IML version is faster than the flint version. IML also checks the result to see if indeed the nullspace it constructed is valid. The flint version doesn't appear to do this.


VERSIONS THAT RELY ON IML
serintode_iml uses IML, the Integer Matrix Library, https://cs.uwaterloo.ca/~astorjoh/iml.html, which relies on GMP, CBLAS, and possibly ATLAS.
IML uses GMP, but not fully. The functions we make use of from IML are kernelMP or nullspaceMP. These functions reduces the input matrix mod p, where p is a random prime. From nullspace.c, we see that these functions both have the following line:
p = RandPrime(15, 19);
The function RandPrime is defined in basisop.c, where the statement above means that the random prime is found between 2^15 and 2^19-1. 

If any input matrix element is 0 mod p, then it assigns it the value p, performing the mod p procedure according to this line in nullspace.c:
DA[i] = (double) ((temp = (A[i] % ((long) p))) >= 0 ? temp : ((long) p) + temp);
We see in the line above that the result is converted to a double, which is why the random prime needs to be in the range 2^15 and 2^19-1, or rather, as stated in the RowEchelonTransform function, it has the precondition:
ceil(n/2)*(p-1)^2+(p-1) <= 2^53-1 = 9007199254740991 (n >= 2)
where n is the number of rows.
This is because IML uses CBLAS in its RowEchelonTransform function, which requires double types.

If IML finds a solution, it checks it, so any solution returned by IML is correct. However, if the check fails, it will attempt a new random prime and repeat the process. In principle, this could go on for a while, but in practice the input matrix would have to be quite pathological. A particular prime p is bad iff the gcd of all minors are divisible by p. Since in this case IML chooses another random prime, the condition that the input matrix will totally fail the algorithm is if the gcd of all minors is divisible by all primes in the range 2^15 and 2^19-1. In practice, if the gcd of all minors is divisible by a sufficient number of these primes, the algorithm could take a long while to find a suitable prime. This would still be rare. The theoretical algorithm is described here:
https://cs.uwaterloo.ca/~astorjoh/jscdense.pdf

The IML version is faster than the flint version.

Compilation:
make
(override IML_DIR=/path/to/iml if IML lives somewhere other than /home/ralf/math/iml)


VERSION THAT RELIES ON FLINT
serintode_flint uses the flint library, http://www.flintlib.org/, which uses fully arbitrary precision integers via GMP for all of its calculations, without using intermediary mod p calculations. In testing, it is slower than IML. For simply determining whether to construct the nullspace, that is, whether COLUMNS-RANK>0, it is only roughly 5% slower. But for the actual construction of the nullspace, it can be 3x slower for a matrix of size 400 or so. IML checks whether the matrix annihilates the nullspace vectors, but flint doesn't appear to do that. I plan to add that in myself. That may be the reason why flint sometimes gives spurious results unless enough checks are given. More testing is needed to determine that question. 

Compilation can be done via:
gcc -Wall serintode_flint.c -o serintode_flint.o -lflint -lgmp -lm -I /usr/local/include/flint