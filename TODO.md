# TODO

Findings from an audit of the C sources. None are exploitable through the
data file — they are robustness, correctness, and quality items.

## Correctness

- [ ] `serintode_flint.c` does not verify that the returned null-space vectors
      actually annihilate the input matrix. README and the file's own header
      comment note that this is why the FLINT version "can give spurious
      results unless enough checks are given." Port IML's annihilation check
      (`kernelMP` / `nullspaceMP` already does it internally) into the FLINT
      path before treating a null vector as a candidate ODE.

## Robustness

- [ ] Check return value of `mpz_init_set_str` / `mpz_set_str` (and
      `fmpz_set_str` in the flint version) when parsing each line. Currently a
      malformed line silently becomes a zero coefficient.
  - `serintode_iml.c:86`, `serintode_iml.c:96`
  - `serintode_flint.c:79`
  - `serintode_iml_nonlin.c:120`, `serintode_iml_nonlin.c:130`
  - `serintode_iml_nonlin_lookup.c` (matching lines)

- [ ] Check `fscanf` return value when reading lookup tables; a truncated file
      leaves `orderexp[i][j]` undefined.
  - `serintode_iml_nonlin_lookup.c:320`

- [ ] Detect lines longer than `MAX_LINE_LENGTH-1` (no newline read by
      `fgets`) and either error out or join the continuation. Currently such a
      line is silently split into two coefficients. Only matters for
      coefficients with ≥99999 digits.

- [ ] `input_string` is a `MAX_LINE_LENGTH+1` (= 100001 byte) stack array in
      every program. Safe today under the default 8 MB Linux stack, but a
      footgun if anyone bumps `MAX_LINE_LENGTH`. Move to `malloc`/`free`.
  - `serintode_iml.c:34`, `serintode_flint.c:36`,
    `serintode_iml_nonlin.c:60`, `serintode_iml_nonlin_lookup.c:36`

## Latent buffer overflows (not reachable from input data)

- [ ] Replace `sprintf` into the 64-byte `fouteqsname` / `fname` stack buffers
      with `snprintf`. Triggers only if `finname` is edited to a long string
      in the source, but the fix is one-line per call site.
  - `serintode_iml.c:354`
  - `serintode_flint.c:318`
  - `serintode_iml_nonlin.c:695`
  - `serintode_iml_nonlin_lookup.c:308`, `serintode_iml_nonlin_lookup.c:681`
  - `makelookup.c:145`

## `makelookup.c`

- [ ] Validate `atol(argv[1])` / `atol(argv[2])` — reject ≤0 and absurdly
      large values before allocating.
- [ ] `(dcheck==-1)&(errno!=EEXIST)` uses bitwise `&`; should be `&&`. Works
      by coincidence today.
  - `makelookup.c:72`, `makelookup.c:81`
- [ ] `mkdir(dirname, S_IRWXU|S_IRWXG|S_IRWXO)` creates the lookup dir as
      0777. Tighten to 0755 (or 0700) unless world-writable is intentional.
  - `makelookup.c:71`

## Ergonomics

- [ ] Take `finname`, `NUM_CHECKS`, `MIN_ODE_ORDER`, `MAX_COEFFS` from `argv`
      / a config rather than recompiling. README already calls this out as a
      known limitation.

## Maintainability

- [ ] The four `serintode_*` programs each carry their own copy of the input
      reader, the output writer, and the main null-space loop (~500–900 LOC
      each). Every fix in this file has to be applied 3–4 times. Factor the
      shared pieces (file parsing, ODE pretty-printing, IML/FLINT wrappers)
      into a small library or a single program with a backend switch.

## Out of scope of this audit

- [ ] `old/serintode_iml.c`, `old/serintode_nonlin2.c`,
      `old/serintode_nonlin_twoterms.c` were not audited. They likely share
      the same parser-return and `sprintf`-into-fixed-buffer issues. Either
      audit them or delete the directory if it's truly archival.
