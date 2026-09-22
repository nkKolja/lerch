//! Sixty-four centered streams in four AVX-512 groups, retaining R=2^64 encoding.

use crate::moments::{PairedSums, Segment};
use crate::reduction::Montgomery32;

pub(crate) fn check_platform() -> Result<(), String> {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx512f") {
        return Ok(());
    }
    Err("AVX-512 requires an x86_64 CPU with AVX-512F; use reference on unsupported CPUs".into())
}

pub(crate) fn paired_sums(
    p: u64,
    batch: u64,
    reducer: Montgomery32,
    quotient_two: u64,
    segments: &[Segment; 64],
    pairs: u64,
) -> Result<PairedSums, String> {
    check_platform()?;
    #[cfg(target_arch = "x86_64")]
    {
        use crate::moments::add;
        debug_assert!(segments.iter().enumerate().all(|(lane, seed)| {
            seed.c < p
                && seed.v < p
                && seed.u < p
                && seed.len == pairs / 64 + u64::from((lane as u64) < pairs % 64)
        }));
        let seeds = std::array::from_fn(|group| native::Seeds {
            c: std::array::from_fn(|lane| segments[group * 16 + lane].c as u32),
            v: std::array::from_fn(|lane| segments[group * 16 + lane].v as u32),
            x: std::array::from_fn(|lane| {
                let seed = segments[group * 16 + lane];
                add(add(seed.u, seed.u, p), seed.v, p) as u32
            }),
        });
        // SAFETY: AVX-512F is detected; the context supplies canonical R64 seeds,
        // balanced lengths, at most (p-1)/2 pairs and a square-safe batch.
        let (first, square) = unsafe {
            native::moments(
                p,
                batch,
                add(quotient_two, quotient_two, p) as u32,
                &seeds,
                pairs / 64,
                (pairs % 64) as usize,
                reducer,
            )
        };
        Ok(PairedSums {
            first,
            square,
            pairs,
        })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (p, batch, reducer, quotient_two, segments, pairs);
        unreachable!("unsupported platforms were rejected")
    }
}

#[cfg(target_arch = "x86_64")]
mod native {
    use super::*;
    use std::arch::x86_64::*;

    #[derive(Clone, Copy)]
    pub(super) struct Seeds {
        pub c: [u32; 16],
        pub v: [u32; 16],
        pub x: [u32; 16],
    }

    #[derive(Clone, Copy)]
    struct State {
        c: __m512i,
        v: __m512i,
        x: __m512i,
        first: __m512i,
        even: __m512i,
        odd: __m512i,
    }

    #[inline]
    #[target_feature(enable = "avx512f")]
    fn normalize(value: __m512i, p: __m512i) -> __m512i {
        let reduce = _mm512_cmp_epu32_mask::<_MM_CMPINT_NLT>(value, p);
        _mm512_mask_sub_epi32(value, reduce, value, p)
    }

    #[inline]
    #[target_feature(enable = "avx512f")]
    fn accumulate(state: &mut State, value: __m512i) {
        let odd = _mm512_srli_epi64::<32>(value);
        let even = _mm512_and_si512(value, _mm512_set1_epi64(u32::MAX as i64));
        state.first = _mm512_add_epi64(state.first, _mm512_add_epi64(even, odd));
        state.even = _mm512_add_epi64(state.even, _mm512_mul_epu32(value, value));
        state.odd = _mm512_add_epi64(state.odd, _mm512_mul_epu32(odd, odd));
    }

    #[inline]
    #[target_feature(enable = "avx512f")]
    fn advance(state: &mut State, twice_qg: __m512i, p: __m512i) {
        let doubled = _mm512_add_epi32(state.c, state.c);
        let carry = _mm512_cmp_epu32_mask::<_MM_CMPINT_NLT>(doubled, p);
        state.c = _mm512_mask_sub_epi32(doubled, carry, doubled, p);
        let odd = _mm512_test_epi32_mask(state.v, _mm512_set1_epi32(1));
        state.v = _mm512_srli_epi32::<1>(_mm512_mask_add_epi32(state.v, odd, state.v, p));
        let correction = _mm512_mask_sub_epi32(state.v, !carry, p, state.v);
        // Normalize BEFORE the third operand. Every individual sum stays below 2p < 2^32.
        let shifted = normalize(_mm512_add_epi32(state.x, twice_qg), p);
        state.x = normalize(_mm512_add_epi32(shifted, correction), p);
    }

    #[inline]
    #[target_feature(enable = "avx512f")]
    fn reduce_block(state: &mut State, reducer: Montgomery32) -> u64 {
        let mut totals = [0u64; 16];
        // SAFETY: each unaligned store writes eight u64 values within totals.
        unsafe {
            _mm512_storeu_si512(totals.as_mut_ptr().cast(), state.even);
            _mm512_storeu_si512(totals.as_mut_ptr().add(8).cast(), state.odd);
        }
        state.even = _mm512_setzero_si512();
        state.odd = _mm512_setzero_si512();
        totals
            .into_iter()
            .map(|total| reducer.reduce_bounded_lazy(total))
            .sum()
    }

    /// Each lane consumes steps terms and the first extra lanes one more.
    /// Q1 and raw bounded REDC totals stay below p*(p-1) < 4e18.
    #[target_feature(enable = "avx512f")]
    pub(super) unsafe fn moments(
        p: u64,
        batch: u64,
        twice_qg: u32,
        seeds: &[Seeds; 4],
        steps: u64,
        extra: usize,
        reducer: Montgomery32,
    ) -> (u64, u64) {
        debug_assert!((3..=crate::moments::MAX_PRIME).contains(&p) && p & 1 == 1);
        debug_assert!(u64::from(twice_qg) < p);
        debug_assert!(seeds.iter().all(|seed| {
            seed.c
                .iter()
                .chain(&seed.v)
                .chain(&seed.x)
                .all(|&value| u64::from(value) < p)
        }));
        debug_assert!(extra < 64 && 2 * (steps * 64 + extra as u64) < p);
        debug_assert!(
            batch != 0 && u128::from(batch) * u128::from(p - 1).pow(2) <= u128::from(u64::MAX)
        );
        let modulus = _mm512_set1_epi32(p as i32);
        let twice_qg = _mm512_set1_epi32(twice_qg as i32);
        let zero = _mm512_setzero_si512();
        let mut states = seeds.map(|seed| State {
            // SAFETY: each array contains sixteen initialized u32 values.
            c: unsafe { _mm512_loadu_si512(seed.c.as_ptr().cast()) },
            v: unsafe { _mm512_loadu_si512(seed.v.as_ptr().cast()) },
            x: unsafe { _mm512_loadu_si512(seed.x.as_ptr().cast()) },
            first: zero,
            even: zero,
            odd: zero,
        });
        let mut square = 0;
        let mut remaining = steps;
        while remaining != 0 {
            let count = remaining.min(batch);
            for _ in 0..count {
                for state in &mut states {
                    accumulate(state, state.x);
                    advance(state, twice_qg, modulus);
                }
            }
            for state in &mut states {
                square += reduce_block(state, reducer);
            }
            remaining -= count;
        }
        if extra != 0 {
            for (group, state) in states.iter_mut().enumerate() {
                let active = extra.saturating_sub(group * 16).min(16);
                let mask = ((1u32 << active) - 1) as __mmask16;
                accumulate(state, _mm512_maskz_mov_epi32(mask, state.x));
                square += reduce_block(state, reducer);
            }
        }
        let first: u64 = states
            .into_iter()
            .map(|state| _mm512_reduce_add_epi64(state.first) as u64)
            .sum();
        (first % p, square % p)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::arith::{fermat_quotient_mod_p, inverse_mod};

        #[test]
        fn genuine_seed_normalizes_before_third_addend() {
            if check_platform().is_err() {
                return;
            }
            let p = 1_999_999_811;
            let reducer = Montgomery32::new(p);
            let x = reducer.encode((2 * fermat_quotient_mod_p(119, p) + inverse_mod(119, p)) % p);
            let v = reducer.encode(inverse_mod(119, p));
            let twice_qg = reducer.encode(2 * fermat_quotient_mod_p(2, p) % p);
            assert_eq!((x, v, twice_qg), (1_945_355_118, 276_697_716, 543_417_792));
            let raw = x + twice_qg + p - v / 2;
            assert_eq!(raw, 4_350_423_863);
            assert_eq!(u64::from(raw as u32) % p, 55_456_567);
            // SAFETY: runtime AVX-512F check, canonical definition-derived state.
            unsafe {
                let zero = _mm512_setzero_si512();
                let mut state = State {
                    c: _mm512_set1_epi32(119),
                    v: _mm512_set1_epi32(v as i32),
                    x: _mm512_set1_epi32(x as i32),
                    first: zero,
                    even: zero,
                    odd: zero,
                };
                advance(
                    &mut state,
                    _mm512_set1_epi32(twice_qg as i32),
                    _mm512_set1_epi32(p as i32),
                );
                let mut cs = [0u32; 16];
                let mut vs = [0u32; 16];
                let mut xs = [0u32; 16];
                _mm512_storeu_si512(cs.as_mut_ptr().cast(), state.c);
                _mm512_storeu_si512(vs.as_mut_ptr().cast(), state.v);
                _mm512_storeu_si512(xs.as_mut_ptr().cast(), state.x);
                assert_eq!(cs, [238; 16]);
                assert_eq!(vs, [138_348_858; 16]);
                assert_eq!(xs, [350_424_241; 16]);
            }
        }

        #[test]
        fn boundary_states_batches_and_tails_match_u64() {
            if check_platform().is_err() {
                return;
            }
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
                    c: std::array::from_fn(|lane| [0, 1, p / 2, p - 1][lane % 4] as u32),
                    v: std::array::from_fn(|lane| [0, 1, p - 2, p - 1][lane % 4] as u32),
                    x: std::array::from_fn(|lane| [p - 1, p - 2, 1, 0][lane % 4] as u32),
                }; 4];
                for requested in [1, 4, 18, 31, 32, 33, 460, 461, 8192] {
                    let batch = crate::moments::square_batch_size(p, requested).unwrap();
                    for len in [
                        0,
                        1,
                        15,
                        16,
                        17,
                        31,
                        32,
                        33,
                        63,
                        64,
                        65,
                        64 * batch - 1,
                        64 * batch,
                        64 * batch + 1,
                        64 * (batch + 1) + 17,
                    ] {
                        let len = len.min((p - 1) / 2);
                        let mut expected = (0, 0);
                        for (group, seed) in seeds.iter().enumerate() {
                            for lane in 0..16 {
                                let (mut c, mut v, mut x) = (
                                    u64::from(seed.c[lane]),
                                    u64::from(seed.v[lane]),
                                    u64::from(seed.x[lane]),
                                );
                                let count =
                                    len / 64 + u64::from(((group * 16 + lane) as u64) < len % 64);
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
                                p,
                                batch,
                                (p - 2) as u32,
                                &seeds,
                                len / 64,
                                (len % 64) as usize,
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
