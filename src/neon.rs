//! Sixteen centered R64 streams, with 32 ordinary carries per u64 division.

use crate::carry_words::CarryWord32;
use crate::moments::{PairedSums, Segment};
use crate::reduction::Montgomery32;

pub(crate) fn check_platform() -> Result<(), String> {
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("neon") {
        return Ok(());
    }
    Err("NEON requires an AArch64 CPU with NEON; use reference on unsupported CPUs".into())
}

pub(crate) fn paired_sums(
    words: CarryWord32,
    batch: u64,
    reducer: Montgomery32,
    quotient_two: u64,
    segments: &[Segment; 16],
    pairs: u64,
) -> Result<PairedSums, String> {
    check_platform()?;
    #[cfg(target_arch = "aarch64")]
    {
        use crate::moments::add;
        let p = words.modulus();
        debug_assert!(segments.iter().enumerate().all(|(lane, seed)| {
            seed.c < p
                && seed.v < p
                && seed.u < p
                && seed.len == pairs / 16 + u64::from((lane as u64) < pairs % 16)
        }));
        let seeds = std::array::from_fn(|group| native::Seeds {
            c: std::array::from_fn(|lane| segments[group * 4 + lane].c as u32),
            v: std::array::from_fn(|lane| segments[group * 4 + lane].v as u32),
            x: std::array::from_fn(|lane| {
                let seed = segments[group * 4 + lane];
                add(add(seed.u, seed.u, p), seed.v, p) as u32
            }),
        });
        // SAFETY: NEON is detected; the context supplies canonical R64 seeds,
        // balanced lengths, at most (p-1)/2 pairs and a square-safe batch.
        let (first, square) = unsafe {
            native::moments(
                words,
                batch,
                add(quotient_two, quotient_two, p) as u32,
                &seeds,
                pairs / 16,
                (pairs % 16) as usize,
                reducer,
            )
        };
        Ok(PairedSums {
            first,
            square,
            pairs,
        })
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        let _ = (
            words.modulus(),
            batch,
            reducer,
            quotient_two,
            segments,
            pairs,
        );
        unreachable!("unsupported platforms were rejected")
    }
}

#[cfg(target_arch = "aarch64")]
mod native {
    use super::*;
    use std::arch::aarch64::*;

    #[derive(Clone, Copy)]
    pub(super) struct Seeds {
        pub c: [u32; 4],
        pub v: [u32; 4],
        pub x: [u32; 4],
    }

    #[derive(Clone, Copy)]
    struct State {
        bits: uint32x4_t,
        v: uint32x4_t,
        x: uint32x4_t,
        first: uint64x2_t,
        lo: uint64x2_t,
        hi: uint64x2_t,
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn add(a: uint32x4_t, b: uint32x4_t, p: uint32x4_t) -> uint32x4_t {
        let sum = vaddq_u32(a, b);
        // Each two-operand sum is below 2p < 2^32. Wrapped sum-p is larger if sum<p.
        vminq_u32(sum, vsubq_u32(sum, p))
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn advance(state: &mut State, twice_qg: uint32x4_t, p: uint32x4_t) {
        // Signed comparison extracts the carry WORD's top bit; it does not compare residues.
        let carry = vcltq_s32(vreinterpretq_s32_u32(state.bits), vdupq_n_s32(0));
        state.bits = vshlq_n_u32::<1>(state.bits);
        let one = vdupq_n_u32(1);
        let odd = vceqq_u32(vandq_u32(state.v, one), one);
        state.v = vshrq_n_u32::<1>(vaddq_u32(state.v, vandq_u32(odd, p)));
        let correction = vbslq_u32(carry, state.v, vsubq_u32(p, state.v));
        // Normalize before adding the third operand: the unreduced sum can exceed u32.
        state.x = add(add(state.x, twice_qg, p), correction, p);
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn accumulate(state: &mut State, value: uint32x4_t) {
        state.first = vpadalq_u32(state.first, value);
        let low = vget_low_u32(value);
        let high = vget_high_u32(value);
        state.lo = vaddq_u64(state.lo, vmull_u32(low, low));
        state.hi = vaddq_u64(state.hi, vmull_u32(high, high));
    }

    #[inline]
    #[target_feature(enable = "neon")]
    fn reduce_block(state: &mut State, reducer: Montgomery32) -> u64 {
        let mut totals = [0u64; 4];
        // SAFETY: each store writes two elements within totals.
        unsafe {
            vst1q_u64(totals.as_mut_ptr(), state.lo);
            vst1q_u64(totals.as_mut_ptr().add(2), state.hi);
        }
        state.lo = vdupq_n_u64(0);
        state.hi = vdupq_n_u64(0);
        totals
            .into_iter()
            .map(|total| reducer.reduce_bounded_lazy(total))
            .sum()
    }

    /// Square blocks and carry words have independent positions: a batch flush
    /// must preserve unconsumed bits, especially with the four-term cap near 2B.
    #[target_feature(enable = "neon")]
    pub(super) unsafe fn moments(
        words: CarryWord32,
        batch: u64,
        twice_qg: u32,
        seeds: &[Seeds; 4],
        steps: u64,
        extra: usize,
        reducer: Montgomery32,
    ) -> (u64, u64) {
        let p = words.modulus();
        debug_assert!(u64::from(twice_qg) < p);
        debug_assert!(seeds.iter().all(|seed| {
            seed.c
                .iter()
                .chain(&seed.v)
                .chain(&seed.x)
                .all(|&value| u64::from(value) < p)
        }));
        debug_assert!(extra < 16 && 2 * (16 * steps + extra as u64) < p);
        debug_assert!(
            batch != 0 && u128::from(batch) * u128::from(p - 1).pow(2) <= u128::from(u64::MAX)
        );
        let modulus = vdupq_n_u32(p as u32);
        let twice_qg = vdupq_n_u32(twice_qg);
        let zero = vdupq_n_u64(0);
        let mut boundary_seeds = seeds.map(|seed| seed.c);
        let mut states = seeds.map(|seed| State {
            bits: vdupq_n_u32(0),
            // SAFETY: each array contains four initialized u32 values.
            v: unsafe { vld1q_u32(seed.v.as_ptr()) },
            x: unsafe { vld1q_u32(seed.x.as_ptr()) },
            first: zero,
            lo: zero,
            hi: zero,
        });
        let mut carry_remaining = 0;
        let mut remaining = steps;
        let mut square = 0;
        while remaining != 0 {
            let count = remaining.min(batch);
            let mut block_remaining = count;
            while block_remaining != 0 {
                if carry_remaining == 0 {
                    for (state, seed) in states.iter_mut().zip(&mut boundary_seeds) {
                        let mut bits = [0u32; 4];
                        for lane in 0..4 {
                            (bits[lane], seed[lane]) = words.next(seed[lane]);
                        }
                        // SAFETY: bits contains four initialized carry words.
                        state.bits = unsafe { vld1q_u32(bits.as_ptr()) };
                    }
                    carry_remaining = 32;
                }
                let chunk = block_remaining.min(carry_remaining);
                for _ in 0..chunk {
                    for state in &mut states {
                        accumulate(state, state.x);
                        advance(state, twice_qg, modulus);
                    }
                }
                carry_remaining -= chunk;
                block_remaining -= chunk;
            }
            for state in &mut states {
                square += reduce_block(state, reducer);
            }
            remaining -= count;
        }
        if extra != 0 {
            for (group, state) in states.iter_mut().enumerate() {
                let active = std::array::from_fn::<_, 4, _>(|lane| {
                    if group * 4 + lane < extra {
                        u32::MAX
                    } else {
                        0
                    }
                });
                // SAFETY: active contains four initialized masks.
                let mask = unsafe { vld1q_u32(active.as_ptr()) };
                accumulate(state, vandq_u32(state.x, mask));
                square += reduce_block(state, reducer);
            }
        }
        // Q1 <= (p-1)^2 and the sum of raw REDCs <= p*(p-1), both < 4e18.
        let first: u64 = states
            .into_iter()
            .map(|state| vaddvq_u64(state.first))
            .sum();
        (first % p, square % p)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::arith::{fermat_quotient_mod_p, inverse_mod};

        #[test]
        fn genuine_seed_normalizes_before_third_addend() {
            check_platform().unwrap();
            let p = 1_999_999_811;
            let reducer = Montgomery32::new(p);
            let x = reducer.encode((2 * fermat_quotient_mod_p(119, p) + inverse_mod(119, p)) % p);
            let v = reducer.encode(inverse_mod(119, p));
            let twice_qg = reducer.encode(2 * fermat_quotient_mod_p(2, p) % p);
            assert_eq!((x, v, twice_qg), (1_945_355_118, 276_697_716, 543_417_792));
            let raw = x + twice_qg + p - v / 2;
            assert_eq!(raw, 4_350_423_863);
            assert_eq!(u64::from(raw as u32) % p, 55_456_567);
            // SAFETY: runtime NEON check, canonical definition-derived state.
            unsafe {
                let zero = vdupq_n_u64(0);
                let mut state = State {
                    bits: vdupq_n_u32(CarryWord32::new(p).next(119).0),
                    v: vdupq_n_u32(v as u32),
                    x: vdupq_n_u32(x as u32),
                    first: zero,
                    lo: zero,
                    hi: zero,
                };
                advance(
                    &mut state,
                    vdupq_n_u32(twice_qg as u32),
                    vdupq_n_u32(p as u32),
                );
                assert_eq!(vgetq_lane_u32::<0>(state.v), 138_348_858);
                assert_eq!(vgetq_lane_u32::<0>(state.x), 350_424_241);
            }
        }

        #[test]
        fn boundary_states_batches_and_carry_tails_match_u64() {
            check_platform().unwrap();
            for p in [
                3,
                9,
                42_447_347,
                999_999_999,
                1_000_000_007,
                1_073_741_827,
                1_431_655_777,
                1_499_999_957,
                1_999_999_973,
                1_999_999_999,
            ] {
                let reducer = Montgomery32::new(p);
                let seeds = [Seeds {
                    c: [0, 1, (p / 2) as u32, (p - 1) as u32],
                    v: [0, 1, (p - 2) as u32, (p - 1) as u32],
                    x: [(p - 1) as u32, (p - 2) as u32, 1, 0],
                }; 4];
                for requested in [1, 4, 18, 31, 32, 33, 460, 461, 8192] {
                    let batch = crate::moments::square_batch_size(p, requested).unwrap();
                    for len in [
                        0,
                        1,
                        15,
                        16,
                        17,
                        16 * 31 - 1,
                        16 * 31,
                        16 * 31 + 1,
                        16 * 32,
                        16 * 32 + 1,
                        16 * 33,
                        16 * 64 + 17,
                        16 * batch - 1,
                        16 * batch,
                        16 * batch + 17,
                    ] {
                        let len = len.min((p - 1) / 2);
                        let mut expected = (0, 0);
                        for (group, seed) in seeds.iter().enumerate() {
                            for lane in 0..4 {
                                let (mut c, mut v, mut x) = (
                                    u64::from(seed.c[lane]),
                                    u64::from(seed.v[lane]),
                                    u64::from(seed.x[lane]),
                                );
                                let count =
                                    len / 16 + u64::from((group * 4 + lane) < (len % 16) as usize);
                                for _ in 0..count {
                                    expected.0 = (expected.0 + x) % p;
                                    expected.1 = (expected.1 + reducer.multiply(x, x)) % p;
                                    let carry = 2 * c >= p;
                                    c = 2 * c % p;
                                    v = (v + (v & 1) * p) / 2;
                                    x = (x + p - 2 + if carry { v } else { p - v }) % p;
                                }
                            }
                        }
                        // SAFETY: canonical seeds, bounded pair count and checked square cap.
                        let got = unsafe {
                            moments(
                                CarryWord32::new(p),
                                batch,
                                (p - 2) as u32,
                                &seeds,
                                len / 16,
                                (len % 16) as usize,
                                reducer,
                            )
                        };
                        assert_eq!(got, expected, "p={p}, len={len}, batch={batch}");
                    }
                }
            }
        }
    }
}
