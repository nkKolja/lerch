# No additional Lerch primes between 200 million and one billion

## Computational result

An optimized Rust search tested **all 39,768,597 primes in
$200{,}000{,}000 < p \le 1{,}000{,}000{,}000$ and found no additional Lerch
primes**. It completed on 14 September 2026 using an eight-core AMD EPYC
9R45 instance. Every completed checkpoint is retained: 8,001 compressed
archives, 886,227,329 bytes of raw compressed results, and one canonical
result tuple per prime.

This is an extension of **Veljko Vranic's upstream search through
200,000,000**, not a new discovery of the fifth Lerch prime and not a fresh
search of the entire range below one billion. Combining the two computations
gives exactly the five already known Lerch primes through that bound:

$$
3,\quad 103,\quad 839,\quad 2237,\quad 42{,}447{,}347.
$$

| Evidence | Interval | Distinct primes | Lerch primes found |
|---|---|---:|---|
| [Upstream v0.3.0](https://github.com/veljkovranic/lerch/releases/tag/v0.3.0), Veljko Vranic | Through 200,000,000 | 11,078,937 | The five listed above |
| [This extension](evidence/publication/audit.json) | Above 200,000,000 through 1,000,000,000 | 39,768,597 | None |
| Combined computational coverage | Through 1,000,000,000 | 50,847,534 | The same five |

The upstream count includes 2, which is not an odd Lerch-prime candidate.
Its six manifest counts sum to 11,078,938 because the prime 32,452,867
appears at a shared boundary; the distinct count removes that duplication.
[The attribution ledger](evidence/publication/upstream.json) preserves the
original intervals, manifest hashes, archive hashes and source revision.
**This machine did not recompute all 50,847,534 primes.**

The contribution is a **finite computational extension and reproducible
implementation**, not a theorem excluding further Lerch primes, a
sublinear algorithm, or a claim of mathematical priority. The complete
large campaign has not been independently repeated by another researcher.

## Definition and arithmetic method

For an odd prime $p$ and an integer $a$ not divisible by $p$, define

$$
q_p(a)=\frac{a^{p-1}-1}{p},\qquad
W_p=\frac{(p-1)!+1}{p},\qquad
\ell_p=\frac{\sum_{a=1}^{p-1}q_p(a)-W_p}{p}.
$$

Lerch's congruence makes $\ell_p$ an integer. A Lerch prime satisfies
$\ell_p\equiv0\pmod p$, equivalently

$$
\sum_{a=1}^{p-1}a^{p-1}-(p-1)!-p\equiv0\pmod{p^3}.
$$

Let $Q_1=\sum q_p(a)\bmod p$ and $Q_2=\sum q_p(a)^2\bmod p$ over the
canonical representatives $a=1,\ldots,p-1$. The moment test is

$$
L_p=\frac{Q_1^2+Q_2-2Q_1}{2}\pmod p,\qquad
p\text{ is Lerch}\ \Longleftrightarrow\ L_p=0.
$$

Division by 2 denotes its modular inverse. To see the identity, expand
$\prod_a(1+p\,q_p(a))=((p-1)!)^{p-1}$ modulo $p^3$ and use
$(p-1)!=-1+pW_p$ and $Q_1\equiv W_p\pmod p$. This avoids computing
an exponentiation for every residue in the search, while the independent
definition-level verifier retains that slower route.

### Doubling cycles cover every odd prime

Write $2c=c'+kp$, with canonical $1\le c,c'<p$ and $k\in\{0,1\}$.
For $u=q_p(c)\bmod p$, $v=c^{-1}\bmod p$ and $G=q_p(2)\bmod p$,

$$
v'=v/2\pmod p,\qquad u'=u+G+kv'\pmod p.
$$

This recurrence does **not** require 2 to be a primitive root.
If $m=\operatorname{ord}_p(2)$ and $d=(p-1)/m$, choose a primitive root
$g$. The $d$ cosets with representatives $1,g,\ldots,g^{d-1}$
partition all nonzero residues into doubling cycles. Representatives are
streamed using the original general-base recurrence; no table of all
residues is needed.

### Centered sign pairs halve the traversal

Here the partner of the canonical positive integer $a$ is the canonical
positive integer **$p-a$**, not the literal negative integer $-a$.
Binomial expansion gives

$$
q_p(p-a)=q_p(a)+a^{-1}\pmod p.
$$

With $x=2u+v$, both members of a sign pair have the same centered state.
Their first- and second-moment contributions are $x$ and
$(x^2+v^2)/2$. For $p>3$, the full inverse-square sum vanishes; hence,
after covering the entire half-domain,

$$
Q_1=\sum_{\text{pairs}}x,\qquad
Q_2=\frac12\sum_{\text{pairs}}x^2\pmod p.
$$

For $p=3$, add 1 to the centered square sum before halving.
Partial square sums are **not** individually treated as partial $Q_2$:
the global conversion follows the check that exactly $(p-1)/2$ pairs
were consumed.

If $m$ is even, $2^{m/2}=-1$: traverse half of every coset.
If $m$ is odd, $d$ is even and negation pairs cosets $r$ and $r+d/2$:
traverse all of the first $d/2$ cosets. Each sign pair is visited once.
The centered update is

$$
c'=2c-kp,\qquad v'=v/2,\qquad
x'=x+2G+(2k-1)v'\pmod p.
$$

The traversal remains linear in $p$ per prime. The improvements are
constant-factor reductions in work and instruction cost.

### Selected ARM and x86 kernels

The publication implementation keeps two production paths: **AArch64 NEON,
16 streams with 32-step carry words**, and **x86-64 AVX-512, 64 streams**.
Workers process different primes; there is no nested thread pool per prime.
Unsupported SIMD selections fail explicitly rather than reporting scalar
work as SIMD. The slower generic recurrence is an explicit reference path.

On ARM, one quotient/remainder calculation generates 32 doubling carries:

$$
2^{32}c=Kp+c_{32},\qquad
K=\left\lfloor 2^{32}c/p\right\rfloor.
$$

The bits of $K$, most significant first, give those carries. Every centered
state and moment update is still performed; this is not a 32-step jump over
the arithmetic. Carry-word consumption is independent of square-block
flushes. The division implementation was faster than the tested reciprocal
variant on the M1 Max, so the latter is not a production option.

Both architectures retain Montgomery radix **$R=2^{64}$**, even though
canonical vector states occupy 32-bit lanes. Widening products accumulate
squares in 64-bit lanes, with per-lane batch cap

$$
B=\min\left(8192,\left\lfloor
\frac{2^{64}-1}{(p-1)^2}\right\rfloor\right).
$$

Reduction is deferred across this bounded block, not performed for each
square. Balanced segments and masked tails preserve exact coverage.
The historical campaign kernels also had scalar cleanup for incomplete
coset groups; the cleaned production kernels use SIMD for those tails.
That implementation change is not retroactively attributed to the campaign.

The cleaned production kernels support **two billion on NEON and
AVX-512**; the historical NEON measurements used its earlier one-billion
limit. For $p\le2\cdot10^9$, doubled states and
modular-halving numerators are below $2p<2^{32}$. A centered update must
normalize before adding its final correction: a single three-operand
sum could overflow. First-moment totals are bounded by $(p-1)^2$ and
deferred square totals by $p(p-1)$, both below $4\cdot10^{18}<2^{64}$.
Coset-seed products fit below $p^2$; exponentiation modulo $p^2$ uses
128-bit products below $p^4<2^{128}$. The square cap is only 4 near
two billion. These bounds, not an unchecked change of an input constant,
justify the AVX-512 extension.

**Current cleaned-code validation (22 September 2026):** native ARM checks
included primes near the two-billion limit; see the
[validation receipt](evidence/publication/current-code-validation.json).
**Cleaned x86 native execution after cleanup: NOT RUN**, because host
authentication is blocked. Prior native AVX-512 evidence belongs to the
archived code, not this cleaned snapshot.

## Completed campaign: machine, time and cost

The producing build was source
`616b5e3f406aba620e9d6780e93275619490c7da`, not the later two-billion
fix or the publication cleanup. Its binary SHA-256 is
`27aa1036c7a32c0e8519c74d88b8511776242a8782495da8f3cb99b76b80151e`.
The separately pinned controller was
`ae192d18e877d1796503af3b114765ff5e3ed1a8`.

| Item | Recorded configuration |
|---|---|
| Instance / region | AWS EC2 c8a.2xlarge, us-east-1 |
| Processor / allocation | AMD EPYC 9R45; 8 allocated cores, 1 hardware thread per core |
| Memory / OS | 16 GiB; x86-64 Ubuntu Linux (26.04 LTS observed at the post-run audit) |
| Toolchain / build | Rust 1.98.1, LLVM 22.1.8; native CPU target, release, thin LTO, one codegen unit |
| Kernel / parallelism | `avx512-64`; 8 outer workers, 1 inner worker |
| Partition | Inclusive 100,000-integer chunks, plus the empty final singleton 1,000,000,000 |
| Logical work | 23,664,304,852,488,604 residue terms; 11,832,152,426,244,302 primary sign-pair steps |
| Search start / final progress | 2026-09-11 21:53:00 UTC / 2026-09-14 03:17:33 UTC |

[The full audit report](evidence/publication/audit.json) supplies the exact
timestamps and unrounded measurements. A later manifest update is not the
search finish time.

| Timing definition | Seconds | Hours | Rate-based compute estimate |
|---|---:|---:|---:|
| Sum of multithreaded chunk-process elapsed times | 186,036.409 | 51.676780 | USD 22.28 |
| Sum of controller-measured chunk wall times | 186,422.066 | 51.783907 | USD 22.32 |
| Campaign start to final progress, including gaps and overhead | 192,273.376 | 53.409271 | USD 23.02 |

These are **elapsed instance times, not summed worker CPU-hours**.
Chunk wall time includes about 6 minutes 26 seconds beyond process
elapsed time. Another 1 hour 37 minutes 31 seconds lies outside those
chunk wall intervals; this residual includes inter-chunk bookkeeping and
interruptions, not a separately measured breakdown of either.

The estimate uses the official **Linux On-Demand c8a.2xlarge rate of
USD 0.43108 per instance-hour**, effective 1 September 2026, from the
21 September price-list publication. The
[selected price-list fields](evidence/publication/aws-price.json) and
[calculation inputs and results](evidence/publication/cost.json) are
retained. This is **not an invoice**; storage, taxes, data transfer,
discounts and other charges are not established by these calculations.

The user reported **approximately USD 64 overall**, including an instance
left running after the search finished. That self-reported total is not
an audited algorithm cost or an itemized development bill.
**Ending a job does not stop EC2 billing.** Use a persisted compute deadline
and a separate, deliberate instance-stop/budget procedure; neither the
historical campaign nor a process deadline implies automatic VM shutdown.

## Measured improvements

All rows below use the **same inclusive interval
200,000,000..200,100,000**, all 5,286 primes, eight outer workers and no
nested parallelism. These are **historical 11 September builds**.
Timing includes fresh pool/sieve/setup, traversal, decoding, sorted
collection and teardown; it excludes compilation, fixture warmups,
hashing, tuple comparison and evidence I/O.

| Host | Method | Samples | Median seconds (min..max) |
|---|---|---:|---:|
| Apple M1 Max | Original generic recurrence | 1 | 794.922981 (single pass, not a median) |
| Apple M1 Max | Centered NEON paired | 3 | 29.570979 (29.570754..29.789657) |
| Apple M1 Max | NEON carry32, division | 3 | 27.839479 (27.830115..27.852997) |
| Apple M1 Max | NEON carry32, reciprocal | 3 | 28.772852 (28.764703..28.867703) |
| AMD EPYC 9R45 | AVX2, 16 streams | 3 | 18.273070 (18.235986..18.363819) |
| AMD EPYC 9R45 | AVX-512, 16 streams | 3 | 16.092824 (16.092035..16.237224) |
| AMD EPYC 9R45 | AVX-512, 32 streams | 3 | 9.330806 (9.323259..9.378327) |
| AMD EPYC 9R45 | AVX-512, 64 streams | 3 | 7.314034 (7.303450..7.511043) |

Every optimized pass matched every six-field tuple from the saved original
generic oracle. M1 methods used rotating order; EPYC methods used recorded,
seeded shuffled order. Exact samples, source revisions and timing scope are
in [the benchmark summary](evidence/publication/benchmarks.json), with
[M1 records](evidence/publication/m1-max-200m.jsonl),
[EPYC records](evidence/publication/epyc-200m.jsonl) and the
[complete generic oracle](evidence/publication/original-200000000-200100000.json).

The **same-M1** original-to-carry32 comparison is **28.554x**, but its
baseline has only one pass. Carry32 improves on paired by **1.062x**.
The **same-EPYC** AVX2-to-AVX-512 comparison is **2.498x**.
M1 carry32 versus EPYC AVX-512 is **3.806x**, while M1 generic versus
EPYC AVX-512 is **108.685x**: both latter comparisons change hardware
and must not be described as algorithm-only speedups.
The EPYC generic full-interval baseline was interrupted; **no completed
same-host EPYC generic timing exists**. These observations do not establish
a universal fastest backend or a confidence interval for performance.

### Two-billion support is not two-billion search coverage

A separately versioned follow-up,
`dc4079671dc81f013c9a543f24b5e8cf63f6e6f0`, tested the narrow interval
**1,999,900,000..2,000,000,000**: 4,586 primes, three identical result
sets, no hits, median **116.084170 seconds** on the eight-core EPYC.
Only two large primes received full original-generic cross-checks; the
whole upper interval did not receive a generic or definition-level rerun.
An alternating old/new 200M comparison measured a **1.253% median
slowdown** for that overflow-safe change.
[Exact samples and scope](evidence/publication/two-billion-bounded.json)
are preserved separately. **No full one-to-two-billion campaign or new
24-hour continuation was launched.** Supported inputs and tested narrow
windows are not exhaustive coverage of that larger range.

## Evidence, reproducibility and limitations

The [standalone Python auditor](scripts/audit_results.py) read every
completed archive and checked all compressed SHA-256 hashes, exact
contiguous chunk boundaries, canonical hashes, field ranges, moment/Lerch
identities, summaries, candidate lists and progress totals. A separate
segmented sieve regenerated the **entire ordered prime list**, detecting
missing, duplicated or reordered primes. The empty composite endpoint is
included. All 8,001 checkpoints passed.

The audit is independent of the Rust sieve and archive writer, but it
**does not recompute all Fermat quotients or certify all primitive roots**.
Hashes establish byte identity, not arithmetic correctness. The original
campaign's chunk flags correctly record that no full generic oracle was
run for those chunks. As an additional bounded cross-check,
[nine fixed stratified primes](evidence/publication/generic-samples/summary.json)
from the start through 999,999,937 matched the original generic recurrence
in all six fields. That reference shares some arithmetic with the optimized
implementation; it is not another researcher's full replication.

The fifth prime and the lower interval are credited to the upstream
project, whose retained verification includes definition-level Rust,
CPython bigint and exact FLINT Bernoulli-number computations. Those are
not new computations performed as part of this extension.

See **[REPRODUCING.md](REPRODUCING.md)** for the data format, complete asset
checksums, safe extraction, independent audit, pinned source builds,
benchmark recipes and full-run commands. The raw campaign is a
**898,508,800-byte release asset**, not a large Git blob. The exact
producing source, lockfile, controller and full MIT notice are separately
archived. Absolute operational paths are removed only from explicitly
labelled derived metadata; original manifests and canonical archives are
unchanged. The package is prepared for private distribution; access to
release assets requires repository permission and publisher upload.

## Build and use the publication implementation

The arithmetic core is Rust; the resumable search supervisor and evidence
tools use Python's standard library. The recorded toolchain is Rust
1.98.1. The reproduction commands use Python 3.12 or newer. Native fast
search requires AArch64 NEON or x86-64 AVX-512F; reference computations
remain available without either. Build on the target machine:

```sh
RUSTFLAGS="-C target-cpu=native" cargo build --release --locked
cargo test --release --locked
python3 -m unittest discover -s scripts -p 'test_*.py'

target/release/lerch-prime-search check --prime 42447347 --backend auto
target/release/lerch-prime-search verify --prime 103 --backend auto
target/release/lerch-prime-search validate --limit 1000 --backend auto
target/release/lerch-prime-search benchmark \
  --start 200000000 --end 200100000 --threads 8 --backend auto \
  --oracle evidence/publication/original-200000000-200100000.json \
  --output benchmark.json
```

Output paths must be new. `benchmark` is the explicit timed range command;
ordinary checks and searches do not trigger hidden benchmark campaigns.
`reference --prime P` selects the slower original recurrence;
`verify --prime P` uses independent definitions, while its explicit
`--generic` option is only a generic-recurrence cross-check.
Do not use a CPU-native binary on a different, unsupported processor.
Current commands produce version-2 evidence; historical version-1 data is
read-only and is never silently migrated or resumed.

## Attribution, license and references

Novak Kaluderovic is responsible for this computational extension and its
reproducible-engineering package. The original software, fifth-prime result
and completed lower-range data are by **Veljko Vranic**. This attribution
does not assert joint authorship of a paper. The original
**MIT license, Copyright (c) 2026 Veljko Vranic**, is preserved in full in
[LICENSE](LICENSE) and in the source capsules. See [CITATION.cff](CITATION.cff)
for software citation metadata.

1. M. Lerch (1905), *Zur Theorie des Fermatschen Quotienten*,
   Mathematische Annalen 60, 471-490.
   [Original reference](https://eudml.org/doc/158206).
2. J. Sondow (2014), *Lerch Quotients, Lerch Primes, Fermat-Wilson
   Quotients, and the Wieferich-non-Wilson Primes 2, 3, 14771*,
   in Combinatorial and Additive Number Theory, 243-255.
   [arXiv:1110.3113v5](https://arxiv.org/abs/1110.3113v5).
   This also records the earlier searches by Marek Wolf.
3. J. B. Dobson (2016), *A characterization of Wilson-Lerch primes*,
   Integers 16, A51.
   [arXiv:1311.2242v7](https://arxiv.org/abs/1311.2242v7);
   earlier versions circulated under the title *A note on Lerch primes*.
4. OEIS Foundation, [A197632](https://oeis.org/A197632), accessed
   22 September 2026: the five terms above and the upstream 200M bound.
   This report does not claim that OEIS has adopted the new 1B bound.
5. V. Vranic (2026), *Lerch prime search and verification*, upstream
   [v0.3.0](https://github.com/veljkovranic/lerch/releases/tag/v0.3.0),
   pinned commit `a42ef4064e0f2f443b74de36d5535c7c40e39736`.
