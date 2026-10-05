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

The original loop, kept in `src/reference_recurrence.rs`, visits all
$p-1$ residues in one sequence by repeatedly multiplying by a primitive
root $g$. Every step performs four hardware divisions: to reduce $gc$, to
update the inverse, to update the quotient, and to add the square to
$Q_2$. The optimized code changes that loop in six places.

1. **Multiply by 2 instead of by $g$** (`src/doubling_cycles.rs`).
   The residue now doubles, so reducing $2c$ needs at most one
   subtraction of $p$. Its inverse halves, which is a conditional
   addition of $p$ and a shift. The quotient update adds $q_p(2)$, plus
   the new inverse when the doubling passed $p$. No step divides.
   Doubling alone reaches every residue only when 2 is a primitive root,
   so the code computes the order of 2 and walks every cycle. The
   original recurrence is used once per cycle, to find each cycle's
   starting point, not once per residue.

2. **Visit one residue of each pair $a$, $p-a$** (`src/moments.rs`).
   Because $q_p(p-a)\equiv q_p(a)+a^{-1}\pmod p$, the loop keeps one value
   per pair, $x=2q_p(a)+a^{-1}$, the sum of both quotients. For $p>3$,
   $Q_1=\sum x$ and $Q_2=\tfrac12\sum x^2$ over one residue per pair;
   $p=3$ gets a separate correction. The loop runs **$(p-1)/2$ steps
   instead of $p-1$**, with the same work per step.

3. **Let one inverse serve two residues** (`src/inverse_grouping.rs`).
   Inverses are stored in Montgomery form with $R=2^{64}$, so the stored
   inverse of $a$ is the number $b=r/a \bmod p$, where $r=2^{64}\bmod p$.
   That number is itself a residue the loop must visit, and its stored
   inverse is $a$. As $a$ doubles, $b$ halves. The kernel therefore keeps
   one pair $(a,b)$ and two quotient values, $x$ for $\{a,p-a\}$ and $y$
   for $\{b,p-b\}$: one residue update and one inverse update serve two
   pairs, so each step covers four residues. Over all 78,498 primes
   through one million on the M1 Max with eight workers, this took
   **0.971 s** versus 1.182 s for the carry-word kernel, **1.22x faster**,
   with every result matching
   ([timings and result digest](evidence/publication/inverse-grouping-1m.json)).

4. **Run 16 to 64 sequences at once with SIMD** (`src/neon.rs`,
   `src/inverse_grouping.rs`, `src/avx512.rs`).
   Each step depends on the previous one, so a single sequence cannot use
   SIMD. The code gives every SIMD lane its own sequence: whole cycles
   when a prime has many, otherwise one cycle cut into equal segments
   whose starting values are computed directly. All lanes of a worker
   share the same prime.

5. **Reduce sums once per block, not every step** (`src/moments.rs`,
   `src/reduction.rs`).
   The original reduced $Q_2$ modulo $p$ after every square. The kernels
   add plain 64-bit squares of the Montgomery-encoded values into 64-bit
   accumulators and apply one Montgomery reduction per block. The block
   length, at most $\lfloor(2^{64}-1)/(p-1)^2\rfloor$ and 8192, rules out
   overflow; it is 18 near one billion. $Q_1$ and the reduced block sums
   also stay in 64-bit accumulators until one final reduction, so the
   sums need no reduction inside the loop.

6. **Get 32 doubling carries from one division** (`src/carry_words.rs`).
   Each doubling step only needs to know whether $2c$ reached $p$. Those
   answers are the binary digits of $c/p$, so dividing $2^{32}c$ by $p$
   gives the next 32 at once, and the remainder is the residue 32 steps
   later. The kernel reads one bit per step, which replaces the per-step
   doubling, comparison and subtraction on $c$ with one division every
   32 steps. In the 200M benchmark this was 1.06x faster than doubling
   $c$ every step on the M1 Max and 1.16x on the M4 Max.

Each backend combines these changes as follows:

| Change | `neon` | `neon-inverse` (ARM default) | `avx512` (x86 default) |
|---|:-:|:-:|:-:|
| 1. Multiply by 2 | yes | yes | yes |
| 2. One residue per pair | yes | yes | yes |
| 3. One inverse for two residues | no | yes | no |
| 4. Sequences per worker | 16 | 16 pairs, as 8 × 2 | 64 |
| 5. Block reduction | yes | yes | yes |
| 6. Carry words | yes | no | no |

`neon-inverse` does not use carry words: it needs the value of $a$
itself to update $y$. Parallelism across primes is as in the original:
a sieve lists the primes, each worker takes whole primes, and completed
intervals are saved as checkpoints.

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
