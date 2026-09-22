/// Ordinary-residue carry generation, independent of the R64 Montgomery encoding.
#[derive(Clone, Copy)]
pub(crate) struct CarryWord32 {
    p: u64,
}

impl CarryWord32 {
    pub(crate) fn new(p: u64) -> Self {
        assert!((3..=crate::moments::MAX_PRIME).contains(&p) && p & 1 == 1);
        Self { p }
    }

    pub(crate) fn modulus(self) -> u64 {
        self.p
    }

    /// The quotient's bits, MSB first, are the next 32 ordinary doubling carries.
    /// The returned residue is AFTER all 32 advances, not after a consumed prefix.
    #[inline(always)]
    #[cfg(any(target_arch = "aarch64", test))]
    pub(crate) fn next(self, c: u32) -> (u32, u32) {
        debug_assert!(u64::from(c) < self.p);
        let t = u64::from(c) << 32;
        debug_assert!(t < 1u64 << 63);
        let word = t / self.p;
        let remainder = t - word * self.p;
        debug_assert!(word <= u64::from(u32::MAX) && remainder < self.p);
        (word as u32, remainder as u32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn bitwise(p: u64, mut c: u64) -> (u32, u32) {
        let mut word = 0u32;
        for _ in 0..32 {
            c *= 2;
            let carry = c >= p;
            if carry {
                c -= p;
            }
            word = (word << 1) | u32::from(carry);
        }
        (word, c as u32)
    }

    #[test]
    fn words_match_every_bit_and_boundary_seed_through_two_billion() {
        let mut rng = 0x5eed_cafe_4244_7347;
        let random_moduli =
            std::array::from_fn::<_, 128, _>(|_| ((random(&mut rng) % 2_000_000_000) | 1).max(3));
        for p in (3..=127)
            .step_by(2)
            .chain([
                103,
                839,
                2237,
                524_287,
                42_447_347,
                200_000_001,
                999_999_999,
                1_000_000_007,
                1_073_741_827,
                1_431_655_777,
                1_499_999_957,
                1_999_999_811,
                1_999_999_973,
                1_999_999_999,
            ])
            .chain(random_moduli)
        {
            let words = CarryWord32::new(p);
            let boundaries = [0, 1, p / 2, p / 2 + 1, p - 2, p - 1];
            for i in 0..256 {
                let c = boundaries
                    .get(i)
                    .copied()
                    .unwrap_or_else(|| random(&mut rng) % p);
                let got = words.next(c as u32);
                assert_eq!(got, bitwise(p, c), "p={p}, c={c}");
                assert_eq!(u64::from(got.0), (c << 32) / p);
                assert_eq!(u64::from(got.1), (c << 32) % p);
                assert_eq!(u64::from(got.0) * p + u64::from(got.1), c << 32);
            }
            assert_eq!(words.next(0), (0, 0));
            assert!(words.next(1).0.leading_zeros() > 0);
        }
    }

    #[test]
    fn consecutive_words_preserve_shifted_streams() {
        for p in [
            3,
            9,
            103,
            2237,
            524_287,
            999_999_999,
            1_999_999_973,
            1_999_999_999,
        ] {
            let words = CarryWord32::new(p);
            for offset in [0, 1, 7, 18, 31, 32, 33, 460, 461] {
                let mut c = p - 1;
                for _ in 0..offset {
                    c = 2 * c % p;
                }
                let mut seed = c as u32;
                for _ in 0..65 {
                    let (word, next) = words.next(seed);
                    for bit in (0..32).rev() {
                        let doubled = 2 * c;
                        assert_eq!((word >> bit) & 1, u32::from(doubled >= p));
                        c = doubled % p;
                    }
                    assert_eq!(u64::from(next), c);
                    seed = next;
                }
            }
        }
    }
}
