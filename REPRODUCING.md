# Reproducing and auditing the Lerch evidence

The [README](README.md) states the scientific result and method.
This guide distinguishes four operations: checking downloaded bytes,
auditing coverage and recorded arithmetic identities, performing bounded
independent arithmetic checks, and actually repeating a large search.
**Only the last would independently repeat the entire computation.**

The archived campaign is complete. Audit its files read-only and use a
new output directory for any reproduction or continuation.

## 1. Obtain and authenticate the evidence

Compact evidence is in `evidence/publication/`. From the repository root:

```sh
shasum -a 256 -c SHA256SUMS
python3 -m unittest discover -s scripts -p 'test_*.py'
```

The publication assets are listed in
[`assets.json`](evidence/publication/assets.json).
All six are attached to the **private**
[`results-1b-v1` release](https://github.com/nkKolja/lerch/releases/tag/results-1b-v1)
in `nkKolja/lerch`, with their uploaded sizes and SHA-256 digests checked
against the retained manifest. An account with repository access can
download them:

```sh
gh release download results-1b-v1 --repo nkKolja/lerch --dir release-data
(cd release-data && shasum -a 256 -c ../evidence/publication/RELEASE-SHA256SUMS)
```

`sha256sum -c` can replace `shasum -a 256 -c` on Linux. Do not proceed if a
checksum differs. The files are comfortably below GitHub's per-asset
2 GiB limit; the large dataset is stored as a release asset.

| Asset | Bytes | Purpose |
|---|---:|---|
| `lerch-200m-1b-campaign-v1.tar` | 898,508,800 | Every completed original chunk, manifest and progress file, plus audit provenance |
| `lerch-producing-source-616b5e3.tar.gz` | 77,518 | Exact Rust crate, Cargo.lock, examples, tests, original license/citation and campaign controller |
| `lerch-producing-source-dc40796.tar.gz` | 75,261 | Exact later two-billion-safe benchmark source |
| `lerch-producing-source-41236ec.tar.gz` | 66,681 | Exact M1 carry-word benchmark source |
| `lerch-producing-source-8eaaab0.tar.gz` | 62,523 | Exact original generic range-oracle harness source |
| `lerch-bounded-two-billion-dc40796.tar.gz` | 1,034,881 | Nine follow-up result files for the 200M and near-2B benchmark windows |

All full asset SHA-256 values are in
[`RELEASE-SHA256SUMS`](evidence/publication/RELEASE-SHA256SUMS).
The principal dataset digest is:

```text
8af31ca5ef02a33102d1199fb34944c13bb4f7cb5d97c570c2025b6b2b108fd1
```

It contains 8,010 **regular files** under `lerch-campaign-v1/`:
8,001 manifest-listed gzip chunks; `manifest.json`, `progress.jsonl`,
`progress.log`; five provenance files; and `SHA256SUMS`. No lock,
orphan/partial checkpoint, morning partial snapshot, deployment
configuration, service log or SSH log is included. The 8,001 gzip files
are unmodified, totalling **886,227,329 bytes**; their uncompressed JSON
totals **4,357,914,763 bytes**.

### Safe extraction

Inspect and validate paths before extraction. The following Python 3.12+
example refuses links, special files, duplicate names and escaping paths,
and requires a new destination. It does not execute any archived script.
Allow roughly 2 GiB of free space for the tar and its still-compressed
chunk files together.

```sh
python3 - <<'PY'
from pathlib import Path, PurePosixPath
import tarfile

archive_path = Path("release-data/lerch-200m-1b-campaign-v1.tar")
destination = Path("extracted-campaign")
if destination.exists():
    raise SystemExit("destination must be new")
with tarfile.open(archive_path) as archive:
    members = archive.getmembers()
    names = set()
    for member in members:
        path = PurePosixPath(member.name)
        if (not member.isfile() or path.is_absolute() or ".." in path.parts
                or not path.parts or path.parts[0] != "lerch-campaign-v1"
                or member.name in names):
            raise SystemExit(f"unsafe archive member: {member.name}")
        names.add(member.name)
    if len(members) != 8010:
        raise SystemExit("unexpected member count")
    destination.mkdir()
    archive.extractall(destination, members=members, filter="data")
PY

(cd extracted-campaign/lerch-campaign-v1 && shasum -a 256 -c SHA256SUMS)
```

Apply the same regular-file/path checks to source capsules before
extraction, using the capsule's own top-level directory from its filename.
The bounded follow-up archive has exactly nine regular basename-only
members: `below1b-old-{1,2,3}.json`, `below1b-new-{1,2,3}.json` and
`near2b-{1,2,3}.json`.

## 2. Audit every completed checkpoint without rerunning the search

Use the reviewed script in this repository, not code executed automatically
from a downloaded archive:

```sh
python3 scripts/audit_results.py extracted-campaign/lerch-campaign-v1 \
  --workers 8 \
  --manifest-sha256 f01b8a20fba9565025a21419423a81ebc675877c4ebfc42f5edbfad254b015ef \
  --output independent-audit.json
```

The output must be outside the input directory and must not already
exist. The auditor does not change any input or launch a Rust computation.
Use fewer workers on a smaller machine. The report's timestamp, Python
version and audit duration will change; its coverage, hashes, totals and
candidate list must match the
[retained report](evidence/publication/audit.json).

The audited result is:

| Field | Expected |
|---|---|
| Campaign format / status | `lerch-campaign-v1` / `complete` |
| Requested inclusive interval | 200,000,000..1,000,000,000 |
| Chunks / progress blocks | 8,001 / 81 |
| Ordered primes | 39,768,597 |
| Residue terms / primary pair steps | 23,664,304,852,488,604 / 11,832,152,426,244,302 |
| Candidate hits | Empty list |
| Next unprocessed integer | 1,000,000,001 |
| Compressed archive-index digest | `9098368d3592e4cf944b22c6e2f6067830526c8d53d51a1a96d8f9f0bb409a2c` |
| Canonical chunk-index digest | `4def4d82a6c50b2bb42c05c145a87cadbf3b3cc6206cfdbd28e37476622c2b7e` |

The separate Python segmented sieve verifies **each actual prime in order**,
not merely a total count. The reader rejects missing/overlapping/duplicated
intervals, missing/reordered primes, wrong sizes or digests, invalid field
types/ranges, inconsistent Lerch remainders or flags, inconsistent summaries,
and incomplete or contradictory progress. The last chunk contains only the
composite integer 1,000,000,000 and has zero rows; dropping it is not accepted.
Tests deliberately damage fixtures and compare serial/parallel results.

The original audit used Python 3.14.4 and eight workers on the EPYC,
taking **15.275 seconds internally** (15.31 seconds under GNU time).
[Resource measurements](evidence/publication/audit-resources.json) are
for the audit, **not the 53-hour search**; the reported maximum RSS is
not the sum of concurrent worker memory.

### Audit coverage and arithmetic cross-checks

The reader checks the stored $Q_1,Q_2,L_p$ identity using Python integers
and independently regenerates the complete prime sequence. The underlying
arithmetic is cross-checked with the generic recurrence at nine selected
primes and with definition-level tests on small inputs. Section 4 gives
those checks; section 6 explains how to repeat the full computation.

## 3. Version-1 data and hash specification

`manifest.json` is the authoritative original `lerch-campaign-v1` record.
Its configuration pins source, binary and controller fingerprints. Chunks
are named `chunks/chunk-LLLLLLLLLL-HHHHHHHHHH.json.gz`; both endpoints
are inclusive. Each archive is a `lerch-range-v1` JSON object:

```text
format, interval {start, end}, requested_batch,
sample {method, threads, timing, summary, full_oracle_tuples_matched},
canonical [{p, primitive_root, q1, q2, lerch_remainder, is_lerch}, ...]
```

Rows are ordered strictly by increasing prime. The four integer residue
fields are canonical; the root lies in `1..p-1`. `lerch_remainder` is
$L_p$ modulo $p$, **not the full integer Lerch quotient**. For odd $p$,
`is_lerch` is exactly `lerch_remainder == 0`.
The special row for 2, supported by the format but absent from this
campaign, has zero moments, null remainder and a false flag.

The canonical digest starts with ASCII `lerch-canonical-v1` followed by
one zero byte. For every row append, in order:

```text
p, primitive_root, q1, q2        four unsigned 64-bit little-endian integers
remainder-present              one byte, 0 or 1
lerch_remainder or zero         one unsigned 64-bit little-endian integer
is_lerch                       one byte, 0 or 1
```

Apply SHA-256 to that byte stream. It is not a hash of JSON formatting.
By contrast, `archive_sha256` is the hash of the **original compressed
gzip bytes**, and the anchored manifest digest is the hash of the exact
manifest file bytes.

The audit adds two aggregate commitments, not present in the historical
manifest. In manifest chunk order, the archive-index digest starts with
`lerch-archive-index-v1` plus a zero byte, then appends
`ARCHIVE_SHA256`, two ASCII spaces, the relative archive path and a newline
for each chunk. The canonical-index digest starts with
`lerch-canonical-index-v1` plus a zero byte, then appends the start and end
as two little-endian `u64` values followed by the **32 decoded digest
bytes** of that chunk's canonical hash. Neither is a flat canonical hash
of all rows concatenated.

The current CLI writes `lerch-range-v2` / `lerch-campaign-v2`, with explicit
backend/build provenance. The six canonical fields and
`lerch-canonical-v1` hash stay compatible. Its oracle reader accepts old
range data, but it does not resume or overwrite an old campaign. The
standalone historical auditor here deliberately accepts **v1 only**.

## 4. Build the exact producing version

The capsule `lerch-producing-source-616b5e3.tar.gz` is distinct from the
current publication source. Its `SOURCE_PROVENANCE.json` lists every
selected file, original commit and byte-level SHA-256. It includes
`Cargo.toml`, `Cargo.lock`, the complete Rust `src/`, `examples/`, Rust
tests, verification helpers, original `CITATION.cff` and unchanged MIT
license. The exact controller is at `controller/epyc_campaign.py`, taken
from its later controller commit. No Git history or deployment configuration
is necessary to build the frozen arithmetic.

| Provenance item | SHA / SHA-256 |
|---|---|
| Producing source commit | `616b5e3f406aba620e9d6780e93275619490c7da` |
| Controller source commit | `ae192d18e877d1796503af3b114765ff5e3ed1a8` |
| Historical binary SHA-256 | `27aa1036c7a32c0e8519c74d88b8511776242a8782495da8f3cb99b76b80151e` |
| Controller file SHA-256 | `a93f7794d27ff23996f5a0ee52e71d5bd865c18dc8c92c57ee9dfad3f8d44eb1` |
| Cargo.lock SHA-256 | `7dc5e0bd9f233be5da367ebc23cfc323b55edfe8d0aee9490960a8066b45c398` |
| Original LICENSE SHA-256 | `5649081d81a1b58b3d313231c5292e564022ab848244569cc0415432657d7f65` |

After verifying and safely extracting the capsule, review the source and
build on a suitable AVX-512F x86-64 host:

```sh
cd lerch-producing-source-616b5e3
shasum -a 256 -c SHA256SUMS
RUSTFLAGS="-C target-cpu=native" cargo +1.98.1 build --release --locked \
  --example cycles_comparison
cargo +1.98.1 test --release --locked
target/release/examples/cycles_comparison \
  range 200000000 200100000 avx512-64 reproduced-200m.json 8
```

Rust 1.98.1 / LLVM 22.1.8, native `znver5`, thin LTO and one codegen unit
were recorded for the EPYC. The source/lockfile capsule preserves
dependencies; it is not a promise that another OS, compiler installation
or linker environment yields a bit-identical executable. A different
build should still reproduce the canonical tuples. The old example
commands are supported **inside these historical capsules**, not by
the cleaned current checkout.

### Bounded cross-checks

The original generic oracle for all 5,286 primes in the 200M benchmark
interval is retained as
`evidence/publication/original-200000000-200100000.json`.
Its file SHA-256 is
`62e2f827b8e5f092c8a8d99fb44745dd8e3d5c6d21ea59d3b1e4d6d80fa81ee5`,
and its canonical digest is
`ed4e5e6f577278e72a61a207041679e42d84d3396d0d224bc83b8bc8d078242b`.
The original generic calculation was performed once on the M1 Max,
not rerun automatically for each optimized benchmark.

The later fixed sample protocol selects the first prime in chunks
0, 1000, ..., 7000 and the last prime in chunk 7999. It checks
200000033, 300000007, 400000009, 500000003, 600000001, 700000001,
800000011, 900000011 and 999999937 with the **full original generic
recurrence**. All matched. Exact raw rows and hashes are under
[`generic-samples/`](evidence/publication/generic-samples).
These concurrent single-worker durations are not benchmark samples.

`scripts/check_archived_samples.py` records that protocol and verifies the
historical binary fingerprint before executing it. A rebuilt binary need
not have that fingerprint; in that case run the same nine explicit
`cycles_comparison range P P generic NEW_FILE.json 1` commands and compare
each entire `canonical` row to its archived row. Do not falsely label a
rebuilt binary as the original. Definition-level checks use the current
`verify` command or the retained independent Python implementation, not
the generic timing oracle.

## 5. Repeat the historical performance experiments

Build each capsule's range harness with the recorded toolchain and native
flags. Use a new output path for every invocation. The source identifiers
and every recorded timing sample are in
[`benchmarks.json`](evidence/publication/benchmarks.json).
The saved oracle avoids silently repeating the expensive generic baseline.

In the M1 carry-word capsule (`41236ec`), on an AArch64 NEON host:

```sh
target/release/examples/cycles_comparison range-compare \
  200000000 200100000 paired,carry32,carry32-reciprocal \
  m1-comparison.jsonl PATH_TO_SAVED_ORACLE.json 3 8
```

In the EPYC capsule (`616b5e3`), on an AVX-512F/AVX2 x86-64 host:

```sh
target/release/examples/cycles_comparison range-compare \
  200000000 200100000 avx2-16,avx512-16,avx512-32,avx512-64 \
  epyc-comparison.jsonl PATH_TO_SAVED_ORACLE.json 3 8
```

To deliberately repeat the one-pass generic baseline, build `8eaaab0`
and select `range ... generic NEW_FILE.json 8`. Do not substitute its
M1 elapsed time for an unmeasured same-EPYC baseline.

For the later `dc40796` capsule, run three separate invocations of:

```sh
target/release/examples/cycles_comparison range \
  1999900000 2000000000 avx512-64 NEW_FILE.json 8
```

Each pass has 4,586 primes. Compare all canonical rows across passes, and
do not call their agreement an independent full generic oracle. The
recorded large-prime generic checks were at 1,499,999,957 and 1,999,999,973.
The raw follow-up archive also retains three old and three new 200M passes
for the overflow-fix comparison. It has **no wider 1B-2B coverage**.

For a measurement of the **current** publication code instead, use
`lerch-prime-search benchmark` as shown in the README. Record its version,
backend, binary hash, machine and full output; do not attach old timings
to that new binary.

## 6. Repeat a full search only with an explicit compute budget

The historical 200M-1B search occupied about **53.4 instance-hours** on
the recorded machine. A faithful full arithmetic repetition is costly;
the much shorter audit above is not a substitute if independent arithmetic
replication is the goal.

For an explicitly approved historical reproduction on a **new directory**,
the frozen controller can run the producing binary:

```sh
# Historical, unbudgeted controller. Review before use; this is not a command
# to run against the original completed campaign or a production service.
python3 controller/epyc_campaign.py \
  --binary target/release/examples/cycles_comparison \
  --source-sha 616b5e3f406aba620e9d6780e93275619490c7da \
  --output NEW_HISTORICAL_REPRODUCTION_DIRECTORY \
  --start 200000000 --end 1000000000 \
  --chunk-width 100000 --progress-width 10000000 \
  --threads 8 --method avx512-64
```

There is **no persisted deadline or instance shutdown in that frozen
controller**. It is supplied for provenance and exact historical
reproduction, not as a recommended unattended cloud launcher. The original
scheduled start and autotuning are not required to reproduce the arithmetic;
new timestamps and compressed/file hashes may differ even when canonical
results agree.

Prefer a separately approved, deadline-limited current-version search:

New searches through 1B persist a default 24-hour deadline if none is
supplied. Above 1B an explicit absolute deadline is required; a supplied
deadline cannot exceed 24 hours after the campaign's original creation.

```sh
# Set DEADLINE_UTC to an explicit, timezone-qualified absolute cutoff
# permitted by the current CLI; never compute a fresh deadline on restart.
target/release/lerch-prime-search search \
  --start 200000000 --end 1000000000 \
  --threads 8 --backend avx512 --chunk-size 100000 \
  --output-dir NEW_PUBLICATION_SEARCH_DIRECTORY \
  --source-sha CURRENT_PUBLICATION_COMMIT \
  --deadline-utc "$DEADLINE_UTC"
```

Use `--resume` only for the **same current-version** output directory and
identity. A persisted deadline is not renewed by restarting. A campaign
paused at its budget limit certifies only its completed prefix; it must
not be reported as a complete 1B search. To authorize later work, use a
new explicit window and directory starting at the recorded next unprocessed
integer. A 24-hour limit is not a promise that a 53-hour search fits in
that window.

The supervisor stops processes, not the EC2 instance. Arrange an explicit
instance-stop action and billing alerts separately, confirm no other job
needs the machine, and verify the instance state after stopping it.

## 7. Preserve upstream credit and evidence

The lower-range result belongs to
[Vranic's upstream v0.3.0](https://github.com/veljkovranic/lerch/releases/tag/v0.3.0),
pinned at `a42ef4064e0f2f443b74de36d5535c7c40e39736`.
[The upstream ledger](evidence/publication/upstream.json) records the six
intervals, byte-level and aggregate manifest hashes, raw archive URLs and
checksums. The two upstream archives are:

```text
1d11f507d43250e04cd75f33960c321d0becbcd17f1cfc706217e5e2801793fe  completed-search-results-through-50000000.tar.zst
e602f3279996bd7129c5fe9077abb248f25709764e4ea3902a65666723a8bd49  completed-search-results-50000000-200000000.tar.zst
```

The first is a tracked asset at the pinned upstream commit; the second
is also a v0.3.0 release attachment. After verifying and safely extracting
them, `scripts/audit_search_results.py` audits their original segment
format, not the newer campaign formats:

```sh
for interval in 2-4496112 4496113-18816869 18816870-18977772 \
                18977773-32452867 32452867-50000000 50000000-200000000
do
  python3 scripts/audit_search_results.py "results/$interval"
done
```

The shared prime endpoint 32,452,867 must be counted only once when
combining these intervals. The shared endpoint 50,000,000 is composite.
This extension relies on the attributed upstream evidence below 200M;
the commands are available for readers who wish to audit it independently.

The fifth-prime verification transcripts and exact Bernoulli verification
are preserved in the upstream
[`evidence/42447347/`](https://github.com/veljkovranic/lerch/tree/a42ef4064e0f2f443b74de36d5535c7c40e39736/evidence/42447347)
and [Bernoulli release](https://github.com/veljkovranic/lerch/releases/tag/bernoulli-42447346).
The retained local `scripts/independent_verify.py` and FLINT wrapper
provide separate arithmetic routes, but large definition-level checks
are deliberately not hidden inside routine commands.

## 8. Packaging and provenance policy

`scripts/package_results.py` requires a matching, anchored full audit
report, validates progress again, and checks every gzip file while
packing. It selects **only** manifest-listed chunks and the documented
metadata. The tar is deterministic for those input bytes: USTAR format,
fixed member order, uid/gid/mtime zero, mode 0644, no links. Gzip chunk
bytes are copied, not recompressed.

```sh
python3 scripts/package_results.py extracted-campaign/lerch-campaign-v1 \
  --audit-report evidence/publication/audit.json \
  --license LICENSE --output NEW_DATASET_COPY.tar
```

With the matching auditor/packager version and retained audit report, this
reproduces the release tar's byte-level hash. A newly timed audit report
would intentionally produce a different tar hash. It must be labelled
as a new package, not passed off as the unchanged original.

`scripts/package_source.py` reads the exact historical Git objects, not a
modified worktree, and packages only the selected complete Rust source
tree plus verification helpers and license. It records a hash for every
file. Its default is the 1B-producing version; `--revision` accepts the
three explicitly retained benchmark revisions. It requires those Git
objects to generate the capsules, but **using the already verified
capsules does not require preserving experimental Git history**.

```sh
python3 scripts/package_source.py --output NEW_PRODUCING_SOURCE_COPY.tar.gz
```

[`provenance.json`](evidence/publication/provenance.json) labels every
transformation of selected metadata. In particular, the EPYC benchmark's
absolute `oracle_path` was replaced by a basename; the original and
derived SHA-256 values are different and both are recorded. Numeric
samples and all other bytes in that JSONL file are unchanged.
Original campaign manifest/progress and canonical gzip archives are
not redacted or rewritten. Private operational logs and orphan/partial
files are excluded rather than republished.

## 9. Mathematical details

### Why two moments suffice

For an odd prime $p$, write

$$
q_p(a)=\frac{a^{p-1}-1}{p},\qquad
W_p=\frac{(p-1)!+1}{p}.
$$

Use the full integer sums $S=\sum_{a=1}^{p-1}q_p(a)$ and
$T=\sum_{a=1}^{p-1}q_p(a)^2$ below. Their residues modulo $p$ are
$Q_1$ and $Q_2$.

The definition gives $1+p q_p(a)=a^{p-1}$. Multiply for every $a$:

$$
\prod_{a=1}^{p-1}(1+p q_p(a))=((p-1)!)^{p-1}.
$$

Expand the left side modulo $p^3$. Terms using three or more
nonconstant factors vanish. The terms using two factors sum to
$\sum_{a<b}q_p(a)q_p(b)=(S^2-T)/2$, giving

$$
1+pS+\frac{p^2}{2}(S^2-T)\pmod{p^3}.
$$

For the right side, $(p-1)!=-1+pW_p$ and $p-1$ is even:

$$
((p-1)!)^{p-1}
=(1-pW_p)^{p-1}
\equiv1-(p-1)pW_p+\binom{p-1}{2}p^2W_p^2
\equiv1+pW_p+p^2(W_p^2-W_p)\pmod{p^3}.
$$

Lerch's congruence makes $\ell_p=(S-W_p)/p$ an integer. Equating the
expansions, cancelling and dividing by $p^2$ gives

$$
2\ell_p+S^2-T\equiv2W_p^2-2W_p\pmod p.
$$

Substitute $S\equiv W_p\pmod p$:

$$
\boxed{2\ell_p\equiv Q_1^2+Q_2-2Q_1\pmod p.}
$$

This moment test is part of the upstream method. The implementation
improvements speed up the calculation of those same two sums.

### Doubling cycles and sign pairs

Write $2c=c'+kp$ with canonical $1\le c,c'<p$ and $k\in\{0,1\}$.
For $u=q_p(c)$, $v=c^{-1}$ and $G=q_p(2)$, all taken modulo $p$,

$$
v'=v/2,\qquad u'=u+G+kv'\pmod p.
$$

Let $m$ be the smallest positive integer with $2^m\equiv1\pmod p$,
and let $d=(p-1)/m$. If $g$ is a primitive root, the cosets with
representatives $1,g,\ldots,g^{d-1}$ give all $d$ doubling cycles.
This covers every nonzero residue even when 2 is not primitive.

The partner of the canonical integer $a$ is the canonical integer
$p-a$, not the literal negative integer $-a$. Binomial expansion gives

$$
q_p(p-a)\equiv q_p(a)+a^{-1}\pmod p.
$$

The centered state $x=2u+v$ is the same for both members of a pair.
Their contributions to the moments are $x$ and $(x^2+v^2)/2$.
For $p>3$, the sum of $v^2$ over one representative of every sign pair
is zero modulo $p$, so

$$
Q_1=\sum_{\text{pairs}}x,\qquad
Q_2=\frac12\sum_{\text{pairs}}x^2\pmod p.
$$

For $p=3$, add 1 to the centered square sum before halving.
The conversion is made only after checking that all $(p-1)/2$ pairs
have been consumed, not separately for arbitrary partial blocks.

If $m$ is even, $2^{m/2}\equiv-1$, so traverse half of every cycle.
If $m$ is odd, $d$ is even and negation pairs cosets $r$ and $r+d/2$;
traverse all of the first $d/2$ cycles. In either case the update is

$$
c'=2c-kp,\qquad v'=v/2,\qquad
x'=x+2G+(2k-1)v'\pmod p.
$$

There is no exponentiation inside this loop.

### Carry words and safe accumulation

For the ARM kernel, compute

$$
2^{32}c=Kp+c_{32},\qquad K=\left\lfloor 2^{32}c/p\right\rfloor.
$$

The bits of $K$, most significant first, are the next 32 doubling
carries. This saves the per-step residue update; it does not skip
the quotient or moment updates.

Both kernels use Montgomery radix $R=2^{64}$. Encoded squares have
scale $R^2$, so each block needs one Montgomery reduction to restore
scale $R$. The safe per-lane block size is

$$
B=\min\left(8192,\left\lfloor
\frac{2^{64}-1}{(p-1)^2}\right\rfloor\right).
$$

For $p\le2\cdot10^9$, a canonical two-operand sum is below
$2p<2^{32}$. The centered update must reduce between additions:
forming a raw three-operand sum first can overflow 32 bits.
Wide first-moment and reduced-square totals are bounded by $(p-1)^2$
and $p(p-1)$, both below $4\cdot10^{18}<2^{64}$.
Seed products are below $p^2$; modular powers modulo $p^2$ use 128-bit
products below $p^4<2^{128}$.

The square block size is only 4 near two billion. These checked bounds
explain the supported limit; raising the limit alone is not safe.
The historical search used scalar cleanup for incomplete coset groups;
the cleaned kernels use SIMD tails. Old benchmark timings are therefore
not silently relabelled as measurements of the new source.
