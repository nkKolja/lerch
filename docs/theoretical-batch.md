# A theoretical batch result for Fermat-quotient moments

Complete $Q_1$ and $Q_2$ for every prime through $N$ can be computed in
$\widetilde O(N^{5/3})$ bit operations, assuming fast integer and
polynomial arithmetic. The construction batches unique small-fraction
representatives, periodic row sums, and product/remainder-tree columns.

## Theorem

Let $N\ge3$, and choose integers $2\le B<N$ and $1\le L\le N$. For every
prime $p\le N$, define the moments over the canonical integers by

$$
q_p(a)=\frac{a^{p-1}-1}{p}\pmod p,\qquad
Q_j(p)=\sum_{a=1}^{p-1}q_p(a)^j\pmod p,\quad j\in\{1,2\}.
$$

**Theorem.** Assuming fast integer and polynomial arithmetic, the complete
values $Q_1(p)$ and $Q_2(p)$ for all primes $p\le N$ can be computed in

$$
\widetilde O\left(\frac{N^2B}{L}+NL+\frac{N^2}{B}+NB\right)
$$

bit operations. Choosing $B\asymp N^{1/3}$ and $L\asymp N^{2/3}$ gives
the bound **$\widetilde O(N^{5/3})$**.

The tilde suppresses logarithmic factors in $N$, including operand bit
lengths, sorting and modular exponentiation costs. Integer multiplication
and division, and polynomial multiplication and power-series inversion,
use fast algorithms. The accounting below includes preprocessing and the
lengths of the polynomial outputs and multiprecision operands.

## Unique small-fraction representatives

Fix the denominator bound $B$ and set $A_p=\lfloor(p-1)/B\rfloor$.
Compute the primes $p\le\max(B,3)$ directly, including
$Q_1(2)=Q_2(2)=0$. For the remaining primes, consider the positive reduced
fractions

$$
1\le a\le A_p,\qquad 1\le b\le B,\qquad \gcd(a,b)=1.
$$

These fractions are injective modulo $p$: the absolute value of a
nonzero cross determinant is less than $p$. Their signs cover every
nonzero residue. To see this, pigeonhole the $B+1$ distinct residues
$0,c,\ldots,Bc$ into $B$ intervals of at most $A_p+1$ integers.
Subtracting two in the same interval gives a numerator of absolute value
at most $A_p$ and a denominator at most $B$; reducing the fraction
preserves these bounds.

Each sign class therefore has one or two positive representatives.
Choose the one with the **smallest denominator**. For $b>1$, an opposite
representative $c/d$ with $d<b$ satisfies $ad+bc=p$, so it must have

$$
d=p\,a^{-1}\bmod b,\qquad
c=\frac{p-ad}{b},\qquad 1\le d<b.
$$

Here $a^{-1}$ is the inverse modulo $b$. The competitor is a positive
reduced fraction, and it lies in the rectangle exactly when $c\le A_p$.
Thus retain $a/b$ precisely when $ad<p-bA_p$. For fixed $b<B$ and
$r=p\bmod b$, eligibility and retention are equivalently

$$
p>aB,\qquad
p\ge ad+b\left\lceil\frac{ad}{B-b}\right\rceil,
\qquad d=r\,a^{-1}\bmod b.
$$

The denominator-one case has only the eligibility threshold $p>aB$.
Denominator $B$ contributes no retained fractions. Sorting these
thresholds once for each $(b,r)$ makes every prime's selected numerator
set a prefix of a fixed list. The least-denominator rule selects exactly
one representative of every sign pair.

## Exact formulas for both moments

Canonical representatives matter because the quotient identities are

$$
q_p(a+kp)=q_p(a)-\frac{k}{a}\pmod p,\qquad
q_p(ab)=q_p(a)+q_p(b)\pmod p.
$$

Let $x_a=q_p(a)$. For $b>1$ and $r=p\bmod b$, define the integer

$$
h(a,b,r)=b-2\bigl((-ar^{-1})\bmod b\bigr),
$$

where $r^{-1}$ is taken modulo $b$; for $b=1$, set $h=1$.
If $z$ is the canonical residue of $a/b\bmod p$, writing $bz=a+kp$
gives $k=(-ar^{-1})\bmod b$ and hence

$$
T(a,b)=x_a-x_b+\frac{h(a,b,r)}{2a}
      =q_p(z)+\frac1{2z}\pmod p.
$$

This centered expression is invariant under $z\mapsto p-z$.
The two members contribute $2T$ to $Q_1$ and
$2T^2+1/(2z^2)$ to $Q_2$. For $p>3$, the inverse-square sum over
the complete half-domain vanishes, so

$$
Q_1=2\sum_{\text{retained }a/b}T(a,b),\qquad
Q_2=2\sum_{\text{retained }a/b}T(a,b)^2\pmod p.
$$

For each numerator row, let $m_a$ be its selected-denominator count and
$H_a$ its integer sum of $h$. For each denominator column, let $n_b$ be
the active prefix length and define

$$
P_b=\prod_{\text{active }a}a,\quad F_b=q_p(P_b),\quad
C_b=\sum_{\text{active }a}\frac h a,\quad
D_b=\sum_{\text{active }a}\frac{h^2}{a^2}.
$$

Expanding the centered sums gives the complete moments:

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

The cross term uses
$\sum_{\text{active }a}q_p(a)=q_p(P_b)$.
The row terms supply the second moments, while the column product
supplies this unweighted first-moment sum.

## Proof of the bit-complexity bound

### Column products and harmonic sums

For each $(b,r)$, product and accumulating remainder trees evaluate the
active prefixes in threshold order. Retain $P_b$ modulo $p^2$ to obtain
its Fermat quotient; the harmonic channels need only modulo $p$.
Writing $U=P\sum h/a$ and $V=P^2\sum h^2/a^2$, the integral state
$(P,P^2,U,V)$ combines two disjoint blocks by

$$
P=P_1P_2,\qquad
U=U_1P_2+U_2P_1,\qquad
V=V_1P_2^2+V_2P_1^2.
$$

This is a constant number of product-tree channels, with empty, shared
and sparse prefix cuts handled by the same construction. Across the
contexts, factor inputs and target moduli total $\widetilde O(NB)$ bits.
Fast product/remainder trees therefore cost $\widetilde O(NB)$ bit
operations, including column preprocessing.

### Periodic row sums and activations

There are $O(N/B)$ numerator rows. Split the range of possible prime
values into blocks of length $L$. In each block, first freeze the masks active at
its start. For each denominator $b$, their count and carry channels
are periodic modulo $b$, with generating-function denominator $1-z^b$
and numerator degree below $b$. A shared denominator tree and inverse
series evaluate their sum in $\widetilde O(B^2+L)$ bit operations per
row block.

The integer row values satisfy $0\le m_a\le B$ and $|H_a|\le B^2$.
A suitably sized auxiliary field therefore gives exact signed
reconstruction with logarithmic-size coefficients.

There are $O(N/L)$ blocks per row, so the frozen periodic evaluations
cost

$$
\widetilde O\left(\frac NB\frac NL(B^2+L)\right)
=\widetilde O\left(\frac{N^2B}{L}+\frac{N^2}{B}\right).
$$

Each $(a,b,r)$ mask activates once. Apply activations inside a block as
explicit arithmetic-progression updates, each costing $O(1+L/b)$.
There are at most $b$ residue classes for denominator $b$, giving

$$
\widetilde O\left(\frac NB\sum_{b\le B}b(1+L/b)\right)
=\widetilde O(NB+NL).
$$

### Quotients, setup and balanced parameters

Evaluating the small-base quotients and assembling the row and column
terms costs $\widetilde O(N^2/B+NB)$ bit operations. The operands have
$O(\log N)$ bits: Fermat quotients are evaluated modulo $p^2$, and the
final moment arithmetic is modulo $p$. Sieving, sorting activation
thresholds and shared precomputation fit within $\widetilde O(NB)$.
The direct treatment of primes at most $\max(B,3)$ also fits this bound.

Adding these costs proves

$$
\widetilde O\left(\frac{N^2B}{L}+NL+\frac{N^2}{B}+NB\right).
$$

With $B\asymp N^{1/3}$ and $L\asymp N^{2/3}$, the first three terms
are $O(N^{5/3})$ and the fourth is $O(N^{4/3})$. This proves the
$\widetilde O(N^{5/3})$ bit-complexity result for the complete collection
of prime moments.

## Implementation and measured performance

The implementation computes complete $Q_1$, $Q_2$ and the Lerch residue
$(Q_2+Q_1^2-2Q_1)/2\bmod p$. It uses bounded packed periodic rows for
**$B\le128$**, with three-field packing for $B\le127$, and a general
NTT/Newton-series path for larger $B$. Cached word remainders serve
small column nodes; GMP handles the wide product/remainder-tree columns.
The implementation also caches small-base quotient data.

The completed comparison covers **all 78,498 primes through
$N=1{,}000{,}000$**, using eight workers on an Apple M1 Max. Each timing
includes fresh setup and preprocessing; compilation, independent
validation and result serialization are excluded. The full moment tuples
matched the carry32 NEON results, with direct-power checks through 10,000.

| Implementation | Samples | Median wall-clock seconds |
|---|---:|---:|
| Cached-word/GMP batch | 3 | 6.196073291 |
| Prior carry32 NEON kernel | 3 | 1.147294291 |

The batch implementation is about **5.4x slower** in this comparison.
The comparator is the prior carry32 kernel, **not the later
inverse-grouped NEON backend**. Exact completed samples, machine and
timing scope are retained in
[`theoretical-batch-results.json`](theoretical-batch-results.json).

The asymptotic saving does not offset setup, polynomial transforms,
memory traffic, and multiprecision product/remainder-tree costs at these
sizes. **For the intended range $p<2^{32}$, the preprocessing, polynomial
and multiprecision overheads make this implementation impractical as an
acceleration of the SIMD search.** This is an implementation-level
assessment; the measured bound is $N=1{,}000{,}000$, and the production
SIMD CLI supports inputs through 2,000,000,000.
