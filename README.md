# Lerch prime search

A faster version of [Veljko Vranic's Lerch-prime search](https://github.com/veljkovranic/lerch).
It keeps the original mathematical test and adds faster ARM NEON and
x86 AVX-512 implementations.

## Result

We checked **every prime above 200 million through one billion**:
**39,768,597 primes, with no new Lerch primes found**.

Veljko's earlier search covered the range through 200 million and found
the fifth Lerch prime, 42,447,347. Combining his result with this search,
the only Lerch primes up to one billion are:

$$
3,\quad 103,\quad 839,\quad 2237,\quad 42{,}447{,}347.
$$

| Search detail | Result |
|---|---|
| Machine | AWS c8a.2xlarge, AMD EPYC 9R45 |
| Resources | 8 physical CPU cores, 16 GiB RAM |
| Implementation | AVX-512, 64 streams per worker, 8 workers |
| Build | Rust 1.98.1, release with native CPU features |
| Elapsed time | **53 hours 24 minutes**, including interruptions |
| Sum of chunk computation times | **51 hours 41 minutes** |
| Completed | 14 September 2026, 03:17 UTC |
| Saved checkpoints | **8,001**, with a result for every prime |

At the published On-Demand price of **$0.43108/hour**, the elapsed run
corresponds to about **$23 in compute**.

## What is faster?

The original code already uses a recurrence instead of computing each
power separately. This version improves that approach:

- **Doubling:** use multiplication by 2 and modular halving in the inner
  loop. Separate cycles cover primes where 2 is not a primitive root.
- **Pairs:** process $a$ and $p-a$ together, cutting the traversal in half.
- **Independent streams:** split cycles into pieces that can run together.
- **SIMD:** process several states per instruction, using NEON on ARM
  and AVX-512 on x86.
- **Batched arithmetic:** keep states in Montgomery form and accumulate
  products in wide registers before reducing them.
- **Deferred sums:** avoid reducing both moment sums after every term.
- **Carry words on ARM:** calculate 32 doubling carries at once, rather
  than calculate a new carry in every step.
- **Parallel search:** give different primes to different CPU workers;
  SIMD runs within each worker, without another layer of threads.

The cleaned code also has resumable searches, deadline limits,
independent verification tools and checked arithmetic bounds up to 2B.

## Measured speed

These measurements use all **5,286 primes from 200M through 200.1M**,
with eight workers on the M1 Max and EPYC.

| Machine | Version | Time |
|---|---|---:|
| M1 Max | Original generic recurrence | 794.923 s |
| M1 Max | Optimized NEON carry-word kernel | 27.839 s |
| EPYC 9R45 | Optimized AVX2 kernel | 18.273 s |
| EPYC 9R45 | Optimized AVX-512 kernel | **7.314 s** |

- **28.55x faster on the same M1 Max.**
- **2.50x faster for AVX-512 versus AVX2 on the same EPYC.**
- **108.68x faster overall**, comparing the original M1 run with the
  optimized EPYC run. That includes the hardware change.

The original baseline is one run; the optimized times are medians of
three runs. These are the recorded September 11 builds, before the
publication cleanup. [Full measurements and source versions](evidence/publication/benchmarks.json)
are retained so the comparisons can be repeated.

## How the test works

For an odd prime $p$, the Fermat quotient is

$$
q_p(a)=\frac{a^{p-1}-1}{p}.
$$

Compute the two sums over $a=1,\ldots,p-1$:

$$
Q_1=\sum_a q_p(a)\pmod p,\qquad
Q_2=\sum_a q_p(a)^2\pmod p.
$$

Then $p$ is a Lerch prime exactly when

$$
Q_1^2+Q_2-2Q_1\equiv0\pmod p.
$$

The fast loop generates the quotients by doubling residues. It uses the
relationship between $a$ and $p-a$ to recover both sums after visiting
only half the residues. See the
[mathematical explanation](REPRODUCING.md#9-mathematical-details)
for the derivation and pairing identity.

## Build and run

Use an AArch64 CPU with NEON or an x86-64 CPU with AVX-512F for the fast
path. The original reference calculation is also available.

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

Choose a new output filename. See `--help` and
[REPRODUCING.md](REPRODUCING.md) for search, resume, deadline and audit
commands. The supervisor uses Python's standard library; the documented
reproduction commands use Python 3.12 or newer.

Both kernels support inputs up to two billion. Current checks cover
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

For cloud runs, set a computation deadline and stop the instance when finished.

## Credit and license

This is an improvement of Veljko Vranic's original software. His code,
lower-range results and discovery of the fifth prime are credited in
[CITATION.cff](CITATION.cff) and the
[upstream evidence record](evidence/publication/upstream.json).
The original [MIT license](LICENSE) is preserved.

For background, see [Sondow's paper](https://arxiv.org/abs/1110.3113),
[Dobson's paper](https://arxiv.org/abs/1311.2242), and
[OEIS A197632](https://oeis.org/A197632).
