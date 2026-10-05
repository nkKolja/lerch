# Lerch prime search

An optimized version of [Vranic's Lerch-prime search](https://github.com/veljkovranic/lerch),
with ARM NEON and x86 AVX-512 implementations. It uses the same mathematical
test, reduces the arithmetic needed for each prime, and extends the
completed search to one billion.

## Result

We checked **every prime above 200 million through one billion**:
**39,768,597 primes, with no new Lerch primes found**.

Vranic's earlier search covered the range through 200 million and found
the fifth Lerch prime, 42,447,347. Combining their result with this search,
the only Lerch primes up to one billion are:

$$
3,\quad 103,\quad 839,\quad 2237,\quad 42{,}447{,}347.
$$

The run used **eight cores of an AMD EPYC 9R45** on an AWS c8a.2xlarge
instance. It completed on 14 September 2026 after **53 hours 24 minutes
of wall-clock time**, including controller overhead and gaps between chunks. At $0.43108 per hour,
this corresponds to approximately **$23 in compute**.

## How the test works

For an odd prime $p$ and an integer $a$ between 1 and $p-1$, Fermat's
little theorem says that $a^{p-1}-1$ is divisible by $p$. The resulting
integer is the **Fermat quotient**:

$$
q_p(a)=\frac{a^{p-1}-1}{p}.
$$

The Lerch test uses two sums: the quotients and their squares, both
reduced modulo $p$:

$$
Q_1=\sum_{a=1}^{p-1}q_p(a)\pmod p,\qquad
Q_2=\sum_{a=1}^{p-1}q_p(a)^2\pmod p.
$$

Then $p$ is a Lerch prime exactly when

$$
Q_1^2+Q_2-2Q_1\equiv0\pmod p.
$$

Evaluating the definition separately for each $a$ would require $p-1$
modular exponentiations. The implementation in Vranic's repository uses
a primitive-root recurrence: repeatedly multiply by a suitable number
$g$ and reduce modulo $p$ to visit every integer from 1 to $p-1$ once.
A recurrence updates the current Fermat quotient from the previous one,
using the multiplication's carry and the residue's modular inverse.
This turns the main calculation into a loop of additions,
multiplications and reductions.

Our implementation computes the same sums with fewer loop steps and
cheaper updates. The
[mathematical details](REPRODUCING.md#9-mathematical-details)
derive the test and the identities used below.

## Improvements to the calculation

1. **Use doubling for cheaper state updates.**
   Replace multiplication by the general value $g$ in the main loop with
   multiplication by 2. Reducing $2a$ modulo $p$ then requires at most one
   subtraction, and updating the inverse requires a conditional addition
   followed by a right shift. When doubling visits only part of the
   nonzero residues, the code follows each of the remaining cycles in
   turn, covering the complete set.

2. **Process a residue and its negative together.**
   The quotients for $a$ and $p-a$ are related by
   $q_p(p-a)\equiv q_p(a)+a^{-1}\pmod p$.
   We track their combined value $x=2q_p(a)+a^{-1}$ and recover the
   complete sums from the values of $x$ and $x^2$.
   For $p>3$, these are $Q_1=\sum x$ and $Q_2=\tfrac12\sum x^2$,
   with one representative per pair; $p=3$ has a separate correction.
   This reduces the traversal from $p-1$ to **$(p-1)/2$ states**.

3. **Share state between inverse pairs on ARM.**
   Let $r=2^{64}\bmod p$. Group the two sign pairs
   $\{a,p-a\}$ and $\{r/a,p-r/a\}$, taking representatives modulo $p$.
   One residue doubles while its partner halves, so the two quotient
   streams share their residue and inverse updates. Each quotient still
   contributes its own square. The `neon-inverse` backend uses this
   grouping; the carry-word implementation remains available as `neon`.
   Over all 78,498 primes through one million on the M1 Max, with eight
   workers, the grouping took **0.971 s** versus 1.182 s for carry words,
   **1.22x faster**, with every result matching
   ([timings and result digest](evidence/publication/inverse-grouping-1m.json)).

4. **Update independent states together with SIMD.**
   One recurrence step depends on its previous state. We split cycles
   into segments and calculate each segment's starting state separately,
   giving several independent streams of work.
   SIMD instructions apply the same operation to multiple streams at
   once. The inverse-grouped ARM kernel processes eight groups of four
   residues in two NEON vector groups. The carry-word ARM backend uses
   16 sign-pair streams, and the x86 kernel uses 64 sign-pair streams.
   All streams within a worker use the same prime and modulus.

5. **Accumulate products before reducing them.**
   Quotient states use Montgomery representation, an encoding that
   allows modular products to be reduced with multiplication and shifts.
   The loop forms ordinary 64-bit square products and sums a block of
   them before applying one Montgomery reduction. The block size is
   chosen from $p$ so the sum fits in 64 bits.
   Wide accumulators also hold the first-moment sum and the reduced
   square-block sums until final normalization. This removes reductions
   from most loop iterations.

6. **Generate doubling carries in blocks in the ARM carry-word backend.**
   Each doubling step needs to know whether $2a$ crossed $p$.
   Dividing $2^{32}a$ by $p$ gives a 32-bit quotient whose bits are the
   next 32 doubling carries. The `neon` kernel consumes these bits one at
   a time while updating the quotient state, replacing the repeated
   residue doubling and comparison with one division per 32 steps.

7. **Run different primes on different workers.**
   A sieve finds the primes in each interval. Workers take separate
   primes from that list, and each worker uses SIMD for its own prime.
   This combines the original implementation's prime-level parallelism with the
   faster inner loop. Completed intervals are saved as checkpoints.

## Measured speed

These measurements use all **5,286 primes from 200M through 200.1M**,
with eight workers on each machine.

| Machine | Version | Wall-clock time |
|---|---|---:|
| M1 Max | Vranic's original implementation | 794.923 s |
| M1 Max | Optimized NEON carry-word kernel | 27.839 s |
| M4 Max | Vranic's original implementation | 598.759 s |
| M4 Max | Optimized NEON carry-word kernel | 21.339 s |
| EPYC 9R45 | Optimized AVX2 kernel | 18.273 s |
| EPYC 9R45 | Optimized AVX-512 kernel | **7.314 s** |

The optimized ARM version is **28.55x faster on the same M1 Max** and
**28.06x faster on the same M4 Max**.
On the EPYC, AVX-512 is **2.50x faster than AVX2**.
Comparing Vranic's original implementation on the M1 Max with the optimized EPYC
version gives **108.68x overall**, combining software and hardware gains.

The original baseline is one run on each Mac; the optimized times are
medians of three runs. The table records the September 11 implementations.
[Measurements and source versions](evidence/publication/benchmarks.json)
are retained for reproduction. The M4 Max rows rebuild those exact sources
on 30 September with Rust 1.95.0; every result matched the M1 Max oracle
([M4 Max measurements](evidence/publication/m4-max-benchmarks.json)).

## Theoretical batch result

Complete $Q_1$ and $Q_2$ for **all primes through $N$** can be computed in
**$\widetilde O(N^{5/3})$ bit operations**, assuming fast integer and
polynomial arithmetic. The construction combines unique small-fraction
representatives, periodic row sums, and product/remainder trees.
[The theorem and proof](docs/theoretical-batch.md) give the full cost
accounting and the balanced parameter choice.

In the completed comparison through **$N=1{,}000{,}000$**, covering all
78,498 primes with eight workers on the M1 Max and fresh setup included,
the batch implementation's median was **6.196073291 s**, versus
**1.147294291 s** for the prior carry32 NEON kernel: about **5.4x slower**.
This comparator is distinct from the later inverse-grouped NEON backend.

Setup, polynomial transforms, memory traffic, and multiprecision
product/remainder-tree costs dominate the practical calculation.
For the intended range $p<2^{32}$, these overheads make this implementation
impractical as an acceleration of the SIMD search. The production commands
therefore use SIMD and retain their two-billion input limit.

## Build and run

Use an AArch64 CPU with NEON or an x86-64 CPU with AVX-512F for the fast
path. `auto` selects inverse grouping on ARM and AVX-512 on x86.
Use `--backend neon` for the previous ARM carry-word implementation.
The original reference calculation is also available.

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release --locked
cargo test --release --locked

# Check one prime with the appropriate fast backend.
target/release/lerch-prime-search check --prime 42447347 --backend auto

# Independently verify a small prime using the definitions.
target/release/lerch-prime-search verify --prime 103 --backend auto

# Time the same interval used in the benchmark table.
target/release/lerch-prime-search benchmark \
  --start 200000000 --end 200100000 --threads 8 --backend auto \
  --oracle evidence/publication/original-200000000-200100000.json \
  --output benchmark.json
```

See `--help` and [REPRODUCING.md](REPRODUCING.md) for search, resume and
audit commands. The search supervisor and evidence tools use Python 3.12
or newer and its standard library.

The fast backends support inputs up to two billion. Current checks cover
native ARM execution, including large-prime comparisons, and x86
cross-compilation. The earlier x86 builds produced the completed search
and EPYC benchmarks.
[Version-specific checks](evidence/publication/current-code-validation.json).

## Data and reproducibility

The [private results release](https://github.com/nkKolja/lerch/releases/tag/results-1b-v1)
contains the complete **899 MB results archive**, exact historical source
versions and checksums. The data is stored as release assets.

All 8,001 checkpoints passed a separate audit of their hashes, prime
coverage and recorded identities. Nine widely spaced primes were also
recomputed with the original recurrence. The original completed data
is preserved unchanged.

## Credit and license

Vranic's code, their lower-range results and their discovery of the
fifth prime are credited in
[CITATION.cff](CITATION.cff) and the
[upstream evidence record](evidence/publication/upstream.json).
The original [MIT license](LICENSE) is preserved.

For background, see [Sondow's paper](https://arxiv.org/abs/1110.3113),
[Dobson's paper](https://arxiv.org/abs/1311.2242), and
[OEIS A197632](https://oeis.org/A197632).
