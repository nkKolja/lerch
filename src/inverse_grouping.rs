//! NEON sign-pair grouping under a -> (2^64 mod p)/a.
//! Selected by `neon-inverse` and by automatic backend detection on AArch64.

use crate::arith::{fermat_quotient_mod_p, inverse_mod, pow_mod};
use crate::moments::{PairedSums, Segment, add, finish};
use crate::{Backend, Canonical, DEFAULT_BATCH_SIZE, DoublingCycleContext};
use serde::Serialize;

pub const KERNEL: &str = "neon8-scaled-inverse-r64";
const LANES: usize = 8;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Work {
    pub sign_pairs: u64,
    pub quartet_steps: u64,
    pub short_orbits: u64,
}

#[derive(Debug, Serialize)]
pub struct Evaluation {
    pub canonical: Canonical,
    pub work: Work,
}

pub fn check_prime(p: u64) -> Result<Evaluation, String> {
    crate::neon::check_platform()?;
    if p == 2 {
        return Ok(Evaluation {
            canonical: Canonical::two(),
            work: Work::default(),
        });
    }
    check(&DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE)?)
}

pub fn check(ctx: &DoublingCycleContext) -> Result<Evaluation, String> {
    crate::neon::check_platform()?;
    if ctx.p() == 3 {
        return Ok(Evaluation {
            canonical: ctx.check(Backend::Neon)?,
            work: Work {
                sign_pairs: 1,
                quartet_steps: 0,
                short_orbits: 1,
            },
        });
    }
    let plan = Plan::new(ctx);
    let mut sums = PairedSums::default();
    let mut work = Work::default();
    let mut pending: Vec<Segment> = Vec::with_capacity(LANES);
    plan.visit(|seed, short| {
        if short {
            let state = plan.state(seed);
            if state.a != state.b && state.a + state.b != ctx.p() {
                return Err("inverse grouping scheduled a non-short orbit as short".into());
            }
            merge(
                &mut sums,
                PairedSums {
                    first: state.x,
                    square: ctx.reducer.multiply(state.x, state.x),
                    pairs: 1,
                },
                ctx.p(),
            );
            work.short_orbits += 1;
        } else {
            if pending.first().is_some_and(|first| first.len != seed.len) {
                flush(&plan, &mut pending, &mut sums, &mut work)?;
            }
            pending.push(seed);
            if pending.len() == LANES {
                flush(&plan, &mut pending, &mut sums, &mut work)?;
            }
        }
        Ok(())
    })?;
    flush(&plan, &mut pending, &mut sums, &mut work)?;
    work.sign_pairs = sums.pairs;
    if 2 * work.quartet_steps + work.short_orbits != work.sign_pairs {
        return Err("inverse grouping work counts do not match sign-pair coverage".into());
    }
    Ok(Evaluation {
        canonical: finish(ctx.p(), ctx.primitive_root(), ctx.reducer, sums)?,
        work,
    })
}

fn merge(total: &mut PairedSums, part: PairedSums, p: u64) {
    total.first = add(total.first, part.first, p);
    total.square = add(total.square, part.square, p);
    total.pairs += part.pairs;
}

fn flush(
    plan: &Plan<'_>,
    pending: &mut Vec<Segment>,
    sums: &mut PairedSums,
    work: &mut Work,
) -> Result<(), String> {
    let mut run = |seeds: &[Segment; LANES], quartets| -> Result<(), String> {
        merge(sums, paired_sums(plan, seeds, quartets)?, plan.ctx.p());
        work.quartet_steps += quartets;
        Ok(())
    };
    if pending.len() == LANES {
        let seeds = std::array::from_fn(|i| pending[i]);
        run(&seeds, LANES as u64 * seeds[0].len)?;
    } else {
        for &seed in pending.iter() {
            run(&plan.ctx.segments(seed), seed.len)?;
        }
    }
    pending.clear();
    Ok(())
}

struct Plan<'a> {
    ctx: &'a DoublingCycleContext,
    n: u64,
    d: u64,
    t: u64,
    h: u64,
}

#[derive(Clone, Copy, Debug)]
struct State {
    a: u64,
    b: u64,
    x: u64,
    #[cfg_attr(not(target_arch = "aarch64"), allow(dead_code))]
    y: u64,
}

impl<'a> Plan<'a> {
    fn new(ctx: &'a DoublingCycleContext) -> Self {
        let p = ctx.p();
        let n = ctx.order() / if ctx.order().is_multiple_of(2) { 2 } else { 1 };
        let g = add(ctx.quotient_two, ctx.quotient_two, p);
        // H = 2r*q(r) = 64G + 2*floor(2^64/p); p is odd.
        let h = (64 * g + 2 * (u64::MAX / p % p)) % p;
        Self {
            ctx,
            n,
            d: (p - 1) / (2 * n),
            t: pow_mod(2, 32, p),
            h,
        }
    }

    fn seed(&self, a: u64, len: u64) -> Segment {
        Segment {
            c: a,
            v: self.ctx.reducer.encode(inverse_mod(a, self.ctx.p())),
            u: self
                .ctx
                .reducer
                .encode(fermat_quotient_mod_p(a, self.ctx.p())),
            len,
        }
    }

    fn state(&self, seed: Segment) -> State {
        let p = self.ctx.p();
        let a = seed.c;
        let b = seed.v;
        let x = add(add(seed.u, seed.u, p), b, p);
        let k = a * b / p;
        // ab = r + kp, so X+Y = H+a+b-2k. All setup arithmetic fits u64.
        let y = (self.h + a + b + 3 * p - 2 * k - x) % p;
        State { a, b, x, y }
    }

    fn visit(
        &self,
        mut consume: impl FnMut(Segment, bool) -> Result<(), String>,
    ) -> Result<(), String> {
        let p = self.ctx.p();
        let mut stream = |a, len, short| {
            if len == 0 {
                Ok(())
            } else {
                consume(self.seed(a, len), short)
            }
        };
        stream(self.t, 1, true)?;
        stream(2 * self.t % p, (self.n - 1) / 2, false)?;
        if self.n.is_multiple_of(2) {
            stream(self.t * pow_mod(2, self.n / 2, p) % p, 1, true)?;
        }
        for seed in self.ctx.representatives_from(
            self.t * self.ctx.primitive_root() % p,
            self.n,
            (self.d - 1) / 2,
        ) {
            consume(seed, false)?;
        }
        if self.d.is_multiple_of(2) {
            if self.n.is_multiple_of(2) {
                let a = self.t * sqrt_two(self.ctx)? % p;
                consume(self.seed(a, self.n / 2), false)?;
            } else {
                let a = self.t * pow_mod(self.ctx.primitive_root(), (p - 1) / 4, p) % p;
                consume(self.seed(a, 1), true)?;
                if self.n > 1 {
                    consume(self.seed(2 * a % p, (self.n - 1) / 2), false)?;
                }
            }
        }
        Ok(())
    }
}

/// Tonelli-Shanks, using the already checked primitive root as a nonresidue.
fn sqrt_two(ctx: &DoublingCycleContext) -> Result<u64, String> {
    let p = ctx.p();
    if pow_mod(2, (p - 1) / 2, p) != 1 {
        return Err("inverse grouping requires a square root of 2 in this coset".into());
    }
    let mut m = (p - 1).trailing_zeros();
    let q = (p - 1) >> m;
    let mut c = pow_mod(ctx.primitive_root(), q, p);
    let mut t = pow_mod(2, q, p);
    let mut x = pow_mod(2, q.div_ceil(2), p);
    while t != 1 {
        let mut i = 0;
        let mut power = t;
        while power != 1 && i < m {
            power = power * power % p;
            i += 1;
        }
        if i == m {
            return Err("inverse grouping square-root setup failed".into());
        }
        let b = pow_mod(c, 1u64 << (m - i - 1), p);
        x = x * b % p;
        c = b * b % p;
        t = t * c % p;
        m = i;
    }
    if x * x % p != 2 {
        return Err("inverse grouping square root is inconsistent".into());
    }
    Ok(x)
}

#[cfg(any(target_arch = "aarch64", test))]
fn quartet_batch_size(p: u64, requested: u64) -> Result<u64, String> {
    use crate::moments::square_batch_size;
    Ok(square_batch_size(p, requested)?.min(square_batch_size(p, u64::MAX)? / 2))
}

fn paired_sums(
    plan: &Plan<'_>,
    seeds: &[Segment; LANES],
    quartets: u64,
) -> Result<PairedSums, String> {
    crate::neon::check_platform()?;
    #[cfg(target_arch = "aarch64")]
    {
        let p = plan.ctx.p();
        let batch = quartet_batch_size(p, plan.ctx.batch_size())?;
        debug_assert_eq!(batch, plan.ctx.kernel_batch_size(Backend::NeonInverse));
        debug_assert!(seeds.iter().enumerate().all(|(lane, seed)| {
            seed.len == quartets / LANES as u64 + u64::from((lane as u64) < quartets % LANES as u64)
        }));
        let states = seeds.map(|seed| plan.state(seed));
        // SAFETY: NEON is detected; seeds are canonical, balanced, and square-capped.
        let (first, square) = unsafe {
            native::moments(
                p,
                batch,
                plan.ctx.reducer,
                add(plan.ctx.quotient_two, plan.ctx.quotient_two, p),
                &states,
                quartets,
            )
        };
        Ok(PairedSums {
            first,
            square,
            pairs: 2 * quartets,
        })
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        let _ = (plan, seeds, quartets);
        unreachable!("unsupported platforms were rejected")
    }
}

#[cfg(target_arch = "aarch64")]
mod native {
    use super::*;
    use crate::reduction::Montgomery32;
    use std::arch::aarch64::*;

    #[derive(Clone, Copy)]
    struct Vectors {
        a: uint32x4_t,
        b: uint32x4_t,
        x: uint32x4_t,
        y: uint32x4_t,
        first: uint64x2_t,
        lo: uint64x2_t,
        hi: uint64x2_t,
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn add_mod(a: uint32x4_t, b: uint32x4_t, p: uint32x4_t) -> uint32x4_t {
        let sum = vaddq_u32(a, b);
        vminq_u32(sum, vsubq_u32(sum, p))
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn advance(s: &mut Vectors, g: uint32x4_t, negative_g: uint32x4_t, p: uint32x4_t) {
        let old_a = s.a;
        let doubled = vaddq_u32(old_a, old_a);
        let carry = vcgeq_u32(doubled, p);
        let one = vdupq_n_u32(1);
        let odd = vceqq_u32(vandq_u32(s.b, one), one);
        s.a = vminq_u32(doubled, vsubq_u32(doubled, p));
        s.b = vshrq_n_u32::<1>(vaddq_u32(s.b, vandq_u32(odd, p)));
        let correction_x = vbslq_u32(carry, s.b, vsubq_u32(p, s.b));
        let correction_y = vbslq_u32(odd, vsubq_u32(p, old_a), old_a);
        s.x = add_mod(add_mod(s.x, g, p), correction_x, p);
        s.y = add_mod(add_mod(s.y, negative_g, p), correction_y, p);
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn accumulate(s: &mut Vectors, x: uint32x4_t, y: uint32x4_t) {
        s.first = vpadalq_u32(s.first, vaddq_u32(x, y));
        s.lo = vmlal_u32(s.lo, vget_low_u32(x), vget_low_u32(x));
        s.lo = vmlal_u32(s.lo, vget_low_u32(y), vget_low_u32(y));
        s.hi = vmlal_u32(s.hi, vget_high_u32(x), vget_high_u32(x));
        s.hi = vmlal_u32(s.hi, vget_high_u32(y), vget_high_u32(y));
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn reduce(s: &mut Vectors, reducer: Montgomery32) -> u64 {
        let mut squares = [0; 4];
        // SAFETY: both stores are within the four-element array.
        unsafe {
            vst1q_u64(squares.as_mut_ptr(), s.lo);
            vst1q_u64(squares.as_mut_ptr().add(2), s.hi);
        }
        s.lo = vdupq_n_u64(0);
        s.hi = vdupq_n_u64(0);
        squares
            .into_iter()
            .map(|x| reducer.reduce_bounded_lazy(x))
            .sum()
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn moments(
        p: u64,
        batch: u64,
        reducer: Montgomery32,
        g: u64,
        seeds: &[State; LANES],
        quartets: u64,
    ) -> (u64, u64) {
        debug_assert!(
            seeds
                .iter()
                .all(|s| { s.a > 0 && s.b > 0 && [s.a, s.b, s.x, s.y].into_iter().all(|x| x < p) })
        );
        debug_assert!(g < p && 4 * quartets < p);
        debug_assert!(
            batch > 0 && u128::from(batch) * 2 * u128::from(p - 1).pow(2) <= u128::from(u64::MAX)
        );
        let modulus = vdupq_n_u32(p as u32);
        let negative_g = vdupq_n_u32(if g == 0 { 0 } else { (p - g) as u32 });
        let g = vdupq_n_u32(g as u32);
        let zero = vdupq_n_u64(0);
        let mut states = std::array::from_fn::<_, 2, _>(|group| {
            let a: [u32; 4] = std::array::from_fn(|i| seeds[4 * group + i].a as u32);
            let b: [u32; 4] = std::array::from_fn(|i| seeds[4 * group + i].b as u32);
            let x: [u32; 4] = std::array::from_fn(|i| seeds[4 * group + i].x as u32);
            let y: [u32; 4] = std::array::from_fn(|i| seeds[4 * group + i].y as u32);
            // SAFETY: each array contains four initialized lanes.
            unsafe {
                Vectors {
                    a: vld1q_u32(a.as_ptr()),
                    b: vld1q_u32(b.as_ptr()),
                    x: vld1q_u32(x.as_ptr()),
                    y: vld1q_u32(y.as_ptr()),
                    first: zero,
                    lo: zero,
                    hi: zero,
                }
            }
        });
        let mut remaining = quartets / LANES as u64;
        let mut square = 0;
        while remaining != 0 {
            let count = remaining.min(batch);
            for _ in 0..count {
                for state in &mut states {
                    accumulate(state, state.x, state.y);
                    advance(state, g, negative_g, modulus);
                }
            }
            for state in &mut states {
                square += reduce(state, reducer);
            }
            remaining -= count;
        }
        let extra = (quartets % LANES as u64) as usize;
        if extra != 0 {
            for (group, state) in states.iter_mut().enumerate() {
                let mask: [u32; 4] = std::array::from_fn(|lane| {
                    if 4 * group + lane < extra {
                        u32::MAX
                    } else {
                        0
                    }
                });
                // SAFETY: mask has four initialized elements.
                let mask = unsafe { vld1q_u32(mask.as_ptr()) };
                accumulate(state, vandq_u32(state.x, mask), vandq_u32(state.y, mask));
                square += reduce(state, reducer);
            }
        }
        let first: u64 = states.into_iter().map(|s| vaddvq_u64(s.first)).sum();
        (first % p, square % p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::moments::half;
    use crate::reference::{direct_invariants, generic};
    use crate::sieve::simple_primes;

    fn advance(state: &mut State, plan: &Plan<'_>) {
        let p = plan.ctx.p();
        let g = add(plan.ctx.quotient_two, plan.ctx.quotient_two, p);
        let old_a = state.a;
        let odd = state.b & 1 != 0;
        let carry = 2 * old_a >= p;
        state.a = 2 * old_a % p;
        state.b = half(state.b, p);
        state.x = (state.x + g + if carry { state.b } else { p - state.b }) % p;
        state.y = (state.y + p - g + if odd { p - old_a } else { old_a }) % p;
    }

    fn assert_state(state: State, plan: &Plan<'_>) {
        let p = plan.ctx.p();
        let reducer = plan.ctx.reducer;
        assert_eq!(state.b, reducer.encode(inverse_mod(state.a, p)));
        for (a, x) in [(state.a, state.x), (state.b, state.y)] {
            assert_eq!(
                x,
                reducer.encode((2 * fermat_quotient_mod_p(a, p) + inverse_mod(a, p)) % p),
                "p={p}, a={a}"
            );
        }
    }

    #[test]
    fn inverse_grouping_exhaustive_coverage_states_and_moments_through_2000() {
        let mut parities = [[0; 2]; 2];
        for p in simple_primes(2000).into_iter().filter(|&p| p > 2) {
            let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).unwrap();
            let plan = Plan::new(&ctx);
            parities[(plan.n % 2) as usize][(plan.d % 2) as usize] += 1;
            let mut seen = vec![false; p as usize];
            let mut shorts = 0;
            let mut pairs = 0;
            plan.visit(|seed, short| {
                let mut state = plan.state(seed);
                for _ in 0..seed.len {
                    assert_state(state, &plan);
                    let is_short = state.a == state.b || state.a + state.b == p;
                    assert_eq!(is_short, short, "p={p}, a={}", state.a);
                    let residues = [state.a, p - state.a, state.b, p - state.b];
                    for &a in &residues[..if short { 2 } else { 4 }] {
                        assert!(!seen[a as usize], "duplicate p={p}, a={a}");
                        seen[a as usize] = true;
                    }
                    pairs += if short { 1 } else { 2 };
                    shorts += u64::from(short);
                    advance(&mut state, &plan);
                }
                Ok(())
            })
            .unwrap();
            assert_eq!(shorts, if p % 4 == 1 { 2 } else { 1 });
            assert_eq!(pairs, (p - 1) / 2);
            assert!(seen[1..].iter().all(|&x| x));
            assert!(!seen[0]);
            let reference = generic(p).unwrap();
            let direct = direct_invariants(p, false, false);
            assert_eq!(
                (reference.q1, reference.q2, reference.lerch_remainder),
                (direct.q1, direct.q2, Some(direct.lerch_remainder))
            );
            if crate::neon::check_platform().is_ok() {
                for batch in [1, 2, 7, 32, DEFAULT_BATCH_SIZE] {
                    let ctx = DoublingCycleContext::new(p, batch).unwrap();
                    let result = check(&ctx).unwrap();
                    assert_eq!(result.canonical, reference, "p={p}, batch={batch}");
                    assert_eq!(result.work.short_orbits, shorts);
                    assert_eq!(result.work.sign_pairs, pairs);
                }
            }
            if pow_mod(2, (p - 1) / 2, p) == 1 {
                let root = sqrt_two(&ctx).unwrap();
                assert_eq!(root * root % p, 2);
            } else {
                assert!(sqrt_two(&ctx).is_err());
            }
        }
        assert!(parities.into_iter().flatten().all(|count| count > 0));
    }

    #[test]
    fn inverse_grouping_high_index_known_hit_and_input_guards() {
        if crate::neon::check_platform().is_err() {
            assert!(check_prime(103).is_err());
            return;
        }
        assert_eq!(check_prime(2).unwrap().canonical, Canonical::two());
        for p in [
            0,
            1,
            4,
            9,
            561,
            1_000_000_001,
            crate::MAX_PRIME,
            crate::MAX_PRIME + 1,
        ] {
            assert!(check_prime(p).is_err());
        }
        for p in [
            3, 5, 7, 13, 17, 31, 127, 257, 2237, 8191, 524_287, 42_447_343, 42_447_347,
        ] {
            let got = check_prime(p).unwrap();
            assert_eq!(got.canonical, crate::check_prime(p, Backend::Neon).unwrap());
            if p <= 524_287 {
                assert_eq!(got.canonical, generic(p).unwrap());
            }
            if p == 42_447_347 {
                assert_eq!(
                    (
                        got.canonical.q1,
                        got.canonical.q2,
                        got.canonical.lerch_remainder,
                        got.canonical.is_lerch
                    ),
                    (34_227_565, 10_415_263, Some(0), true)
                );
            }
        }
    }

    #[test]
    fn inverse_grouping_shifted_tails_and_square_caps_through_two_billion() {
        for p in [
            103,
            2237,
            42_447_347,
            999_999_893,
            1_000_000_007,
            1_073_741_827,
            1_431_655_777,
            1_499_999_957,
            1_999_999_811,
            1_999_999_973,
        ] {
            let cap = quartet_batch_size(p, u64::MAX).unwrap();
            let square = 2 * u128::from(p - 1).pow(2);
            assert!(u128::from(cap) * square <= u128::from(u64::MAX));
            assert!(u128::from(cap + 1) * square > u128::from(u64::MAX));
            let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).unwrap();
            let plan = Plan::new(&ctx);
            for a in [1, 2, 119 % p, p / 2, p - 1, pow_mod(2, 33, p)] {
                for len in [
                    0, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 287, 288, 289, 511, 512,
                    513, 1041,
                ] {
                    let len = len.min((p - 1) / 4);
                    let seed = plan.seed(a, len);
                    let mut state = plan.state(seed);
                    let mut expected = PairedSums::default();
                    for _ in 0..len {
                        assert_state(state, &plan);
                        expected.first = (expected.first + state.x + state.y) % p;
                        expected.square = (expected.square
                            + ctx.reducer.multiply(state.x, state.x)
                            + ctx.reducer.multiply(state.y, state.y))
                            % p;
                        advance(&mut state, &plan);
                    }
                    if crate::neon::check_platform().is_ok() {
                        for batch in [1, 2, 3, 7, 31, 32, 33, DEFAULT_BATCH_SIZE] {
                            let ctx = DoublingCycleContext::new(p, batch).unwrap();
                            let plan = Plan::new(&ctx);
                            let got = paired_sums(&plan, &ctx.segments(seed), len).unwrap();
                            assert_eq!(
                                (got.first, got.square, got.pairs),
                                (expected.first, expected.square, 2 * len),
                                "p={p}, a={a}, len={len}, batch={batch}"
                            );
                        }
                    }
                }
            }
        }
        assert_eq!(
            quartet_batch_size(crate::MAX_PRIME, DEFAULT_BATCH_SIZE),
            Ok(2)
        );
        assert!(quartet_batch_size(7, 0).is_err());
    }

    #[test]
    fn inverse_grouping_genuine_overflow_state() {
        let p = 1_999_999_811;
        let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).unwrap();
        let plan = Plan::new(&ctx);
        let mut state = plan.state(plan.seed(119, 1));
        let g = add(ctx.quotient_two, ctx.quotient_two, p);
        assert_eq!(
            (state.x, state.b, g),
            (1_945_355_118, 276_697_716, 543_417_792)
        );
        assert_eq!(state.x + g + p - state.b / 2, 4_350_423_863);
        advance(&mut state, &plan);
        assert_eq!(state.x, 350_424_241);
        assert_state(state, &plan);
    }
}
