# Theoretical batch method and measured limitations

This note records a derived experimental alternative to the SIMD Lerch search;
no literature-priority claim is made.
It gives a complete batch formula, not a faster production backend.
The best completed one-million experiment was **6.196 s**, compared with
**1.147 s** for the unchanged SIMD implementation. At nine million, SIMD
completed in **67.156 s**; the batch run was **cancelled by the user after
72 min 20.548 s**, still in its row phase. There is no completed
nine-million batch result.

The compact measurements are in
[`theoretical-batch-results.json`](theoretical-batch-results.json).
This documentation does not require importing the experimental native
dependencies into the SIMD production build.

## Exact reduction to unique small fractions

For a prime $p$, write

$$
q_p(a)=\frac{a^{p-1}-1}{p}\pmod p,\qquad
Q_j(p)=\sum_{a=1}^{p-1}q_p(a)^j\pmod p.
$$

The representatives in this definition are the canonical integers
$1,\ldots,p-1$. Changing their lifts is not harmless:
$q_p(a+kp)=q_p(a)-k/a$, whereas multiplication gives
$q_p(ab)=q_p(a)+q_p(b)$.

For all primes through $N$, fix a denominator bound $B\ge2$ and set
$A_p=\lfloor(p-1)/B\rfloor$. Treat $p\le\max(B,3)$ separately.
The positive reduced fractions

$$
1\le a\le A_p,\qquad 1\le b\le B,\qquad \gcd(a,b)=1
$$

are injective modulo $p$, since the absolute value of the determinant
of two such fractions is less than $p$. Their signs cover every nonzero
residue: pigeonhole the $B+1$ multiples $0,c,\ldots,Bc$ into $B$ intervals
of at most $A_p+1$ integers, subtract two in the same interval, and reduce
the resulting fraction.

Each sign class therefore has one or two positive representatives.
Choose its representative with the **smallest denominator**. For $b>1$,
a smaller-denominator competitor must have

$$
d=p\,a^{-1}\bmod b,\qquad c=\frac{p-ad}{b},\qquad 1\le d<b.
$$

It is in the rectangle exactly when $c\le A_p$. Thus retain $a/b$ precisely
when $ad<p-bA_p$. For fixed $b<B$ and $r=p\bmod b$, this is a single
activation threshold along that residue class:

$$
p>aB,\qquad p\ge ad+b\left\lceil\frac{ad}{B-b}\right\rceil,
\qquad d=r\,a^{-1}\bmod b.
$$

The denominator-one case has only the eligibility threshold; denominator
$B$ contributes no retained fractions. Sorting thresholds once per
$(b,r)$ makes each prime's selected numerator set a prefix of a fixed list.
There is no unresolved collision correction.

Let $x_a=q_p(a)$ and define the integer

$$
h(a,b,r)=b-2\bigl((-ar^{-1})\bmod b\bigr).
$$

For the canonical residue $z=a/b\bmod p$, writing $bz=a+kp$ gives

$$
T(a,b)=x_a-x_b+\frac{h(a,b,r)}{2a}
      =q_p(z)+\frac1{2z}.
$$

This centered expression is invariant under $z\mapsto p-z$.
For $p>3$, the inverse-square sum vanishes, so
$Q_1=2\sum T$ and $Q_2=2\sum T^2$ over the retained fractions.

For each numerator row let $m_a$ be its selected-denominator count and
$H_a$ its integer sum of $h$. For each denominator column let
$n_b$ be the active prefix length and define

$$
P_b=\prod_{\text{active }a}a,\quad F_b=q_p(P_b),\quad
C_b=\sum_{\text{active }a}\frac h a,\quad
D_b=\sum_{\text{active }a}\frac{h^2}{a^2}.
$$

Then the complete moments are

$$
Q_1=2\sum_b(F_b-n_bx_b)+\sum_b C_b,
$$

$$
Q_2=
2\sum_a m_ax_a^2+2\sum_b n_bx_b^2
-4\sum_b x_bF_b+2\sum_a\frac{x_aH_a}{a}
-2\sum_b x_bC_b+\frac12\sum_bD_b
\pmod p.
$$

The cross term is tractable because
$\sum_{\text{active }a}q_p(a)=q_p(P_b)$. This uses an ordinary product
for an unweighted first moment; it does not recover second moments from
a product alone.

## Full cost, including preprocessing

Columns use constant-size product and accumulating remainder trees in
activation order. The product $P_b$ is retained modulo $p^2$ for its
Fermat quotient; the harmonic channels need modulo $p$. A state
$(P,P^2,P\sum h/a,P^2\sum h^2/a^2)$ composes by fixed-size polynomial
identities. Empty, shared and sparse prefix cuts are handled exactly.
Across contexts, factor input and target-modulus sizes total
$\widetilde O(NB)$ bits.

Rows use blocks of length $L$. Freeze masks already active at the block
start. Their periodic generating functions have denominators $1-z^b$
and numerator degrees below $b$. A shared denominator tree and inverse
series evaluate the count and carry channels in
$\widetilde O(B^2+L)$ work per row block. Apply activations inside the
block by explicit arithmetic-progression updates. The row integers have
known bounds, so a suitably sized auxiliary field gives exact signed
reconstruction, not an unresolved CRT computation.

There are $O(N/B)$ rows and $O(N/L)$ blocks per row. Accounting for the
periodic evaluations, event updates, all small-base quotient evaluations,
column trees, sieving, sorting and shared precomputation gives

$$
\widetilde O\left(\frac{N^2B}{L}+NL+\frac{N^2}{B}+NB\right)
$$

bit operations, **assuming fast integer and polynomial arithmetic**.
The output-series length $L$ and all modulus bit lengths are counted.
Choosing $B\asymp N^{1/3}$ and $L\asymp N^{2/3}$ gives
$\widetilde O(N^{5/3})$. This is neither the requested
$\widetilde O(N^{3/2})$ bound nor a single-prime square-root algorithm.

A streaming implementation can have near-linear working space, but the
measured speed-oriented implementation additionally caches prime-base
quotients. Its explicit memory cost is not covered by that streaming
space statement. The asymptotic argument does not establish a practical
crossover against SIMD.

## What was implemented

The research code computes complete $Q_1$, $Q_2$ and the Lerch residue
$(Q_2+Q_1^2-2Q_1)/2\bmod p$. Completed million-bound sweeps matched every
one of the 78,498 full tuples against the existing NEON implementation,
with direct-power checks through 10,000.

Implementations included FLINT polynomial rows, a Rust NTT/Newton-series
row engine, bounded packed periodic sums, three-coefficient word packing,
LCM grouping, sparse output access, explicit NEON register tiling and
validated input-cycle reuse. The fast periodic paths are gated by
**$B\le128$** (three-field packing by $B\le127$); larger bounds use the
general NTT path.

Column variants included GMP, pure-Rust `num-bigint`, mixed precision,
compressed prefix cuts, cached word Montgomery remainders and proved
omission of modulus products larger than every possible exact prefix.
The Rust-only version improved rows but slowed the wide columns.
An optional GMP hybrid retained the faster overall measured combination.

The investigation also tested the MIT-licensed
`nkKolja/Prime-field-arithmetic` division routines, pinned at
`c10836d835028746bee8c9d364a385662515c71f`. Factoring divisor normalization
and the 3/2 preinverse out of each call passed extensive correctness and
sanitizer checks. Warm-kernel wins did not survive preparation and safe
conversion at actual node reuse counts; the dispatcher was **not enabled**.
A separate wide-Barrett experiment likewise lost after cold setup.

## Measurements and the nine-million cancellation

All listed full runs used eight workers on the same Apple M1 Max.
Cold preprocessing is included; compilation, independent validation and
result serialization are excluded. The initial 178.736 s build used
ordinary optimized release; later comparisons used native CPU flags.
The SIMD comparator is the unchanged production NEON16 kernel, not the
separate inverse-grouping optimization.

| Bound | Implementation | Complete wall seconds |
|---:|---|---:|
| 1,000,000 | Initial FLINT/native batch | 178.736 |
| 1,000,000 | Optimized rows with GMP columns | 6.916 |
| 1,000,000 | Rust-only bounded columns | 8.778 |
| 1,000,000 | Latest cached-word/GMP hybrid, median of three cold runs | **6.196** |
| 1,000,000 | Same-comparison old NEON, median of three cold runs | **1.147** |
| 9,000,000 | Old NEON, one completed cold run | **67.156** |
| 9,000,000 | Default-parameter batch | **User-cancelled; no completed result** |

The latest completed million-bound hybrid is about **5.4 times slower
than old NEON**, despite a large improvement over the first batch port.
Pure-Rust wide arithmetic and the rejected division alternatives must not
be described as faster than the measured GMP paths.

At nine million the unchanged defaults are $B=209$, $L=43{,}371$ and
43,062 numerator rows. This crosses the packed-path limit and invokes
the general NTT fallback. Its prime-basis cache contains
1,382,125,257 `u32` entries: **5,528,501,028 bytes (5.15 GiB)**.
That allocation and approximately 17.3 seconds of cold setup were charged.
Three representative NTT blocks took **9.31–10.42 ms each**; these are
component observations, not a completed sweep time.

The user stopped the actual batch on **2026-09-23 at 22:02 CEST**.
The owned process and wrapper exited and were not restarted.
The receipt records **4340.547991583 seconds (72 min 20.548 s)** elapsed,
termination by SIGTERM, and **not** a timeout. Peak process RSS was
5,874,647,040 bytes (5.47 GiB).

The final logged progress, at 4337.8 seconds, was **30,848 of 43,062 rows**
and **4,123,866 of 4,489,413 row blocks (91.86%)**. Those are completed
**row-work blocks**, not completed primes or final moment coverage.
The column phase and final moment assembly had not started. No complete
nine-million batch tuple file or comparison exists.

The nine-million SIMD oracle is complete and remains valid: 602,489
primes, last prime 8,999,993, computed hits `[3,103,839,2237]`.
A hybrid progress-only correction changed the executable hash after the
SIMD measurement; the baseline function and SIMD kernel files were checked
byte-identical and both measured binaries/source snapshots were preserved.
The SIMD sweep was not rerun.

These results do **not** justify an overnight billion-bound search or a
week-long four-billion search. The cache, arithmetic bounds and fallback
regime require separate validation, and the cancelled run cannot supply
an extrapolated success claim. The experimental dependency tree is not
needed for, or promoted into, production SIMD search.
