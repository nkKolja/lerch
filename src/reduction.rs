use crate::arith::mul_mod;

/// Montgomery reducer with R=2^64, specialized to odd moduli below 2^32.
#[derive(Clone, Copy)]
pub struct Montgomery32 {
    modulus: u64,
    negative_inverse: u64,
    r2: u64,
}

impl Montgomery32 {
    pub fn new(modulus: u64) -> Self {
        assert!(modulus > 1 && modulus & 1 == 1 && modulus < (1u64 << 32));
        let mut inverse = 1u64;
        for _ in 0..6 {
            inverse = inverse.wrapping_mul(2u64.wrapping_sub(modulus.wrapping_mul(inverse)));
        }
        let r = ((1u128 << 64) % modulus as u128) as u64;
        Self {
            modulus,
            negative_inverse: inverse.wrapping_neg(),
            r2: mul_mod(r, r, modulus),
        }
    }

    #[inline(always)]
    fn reduce(self, product: u128) -> u64 {
        let correction = (product as u64).wrapping_mul(self.negative_inverse);
        let sum = product + correction as u128 * self.modulus as u128;
        let mut result = (sum >> 64) as u64;
        if result >= self.modulus {
            result -= self.modulus;
        }
        result
    }

    /// Reduce a bounded sum of encoded products, preserving Montgomery encoding.
    ///
    /// Every u64 total T satisfies T < p*R for p >= 3 and R = 2^64.
    /// The REDC numerator T + correction*p is below 2^64 + 2^96, so
    /// u128 cannot overflow and one conditional subtraction is sufficient.
    /// Callers must bound accumulation before forming T; this does not use T % p.
    #[inline(always)]
    #[cfg(test)]
    fn reduce_accumulated_products(self, total: u64) -> u64 {
        self.reduce(u128::from(total))
    }

    /// R64 REDC for T < R, including composite odd moduli.
    /// The low-word carry is exactly T != 0; the quotient lies in 0..=p.
    #[inline(always)]
    #[cfg(test)]
    fn reduce_bounded(self, total: u64) -> u64 {
        let result = self.reduce_bounded_lazy(total);
        if result == self.modulus { 0 } else { result }
    }

    /// Congruent to REDC(T), but p is a permitted representation of zero.
    #[inline(always)]
    pub(crate) fn reduce_bounded_lazy(self, total: u64) -> u64 {
        let correction = total.wrapping_mul(self.negative_inverse);
        let high = ((u128::from(correction) * u128::from(self.modulus)) >> 64) as u64;
        high + u64::from(total != 0)
    }

    #[inline(always)]
    pub fn encode(self, value: u64) -> u64 {
        self.reduce(value as u128 * self.r2 as u128)
    }

    #[inline(always)]
    pub fn multiply(self, a: u64, b: u64) -> u64 {
        self.reduce(a as u128 * b as u128)
    }

    #[inline(always)]
    pub fn decode(self, value: u64) -> u64 {
        self.reduce(value as u128)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn bounded_redc_random_and_boundary_bigint_oracle() {
        let mut random = 0xd12a_7347_5eed_1234_u64;
        let random_moduli = std::array::from_fn::<_, 64, _>(|_| {
            random ^= random << 13;
            random ^= random >> 7;
            random ^= random << 17;
            (u64::from(random as u32) | 1).max(3)
        });
        for p in [
            3,
            5,
            9,
            15,
            103,
            42_447_347,
            999_999_893,
            999_999_929,
            999_999_999,
            1_073_741_823,
            1_431_655_765,
            1_499_999_957,
            1_999_999_973,
            1_999_999_999,
            3_999_999_959,
            u64::from(u32::MAX),
        ]
        .into_iter()
        .chain(random_moduli)
        {
            let reducer = Montgomery32::new(p);
            let modulus = BigUint::from(p);
            let radix = BigUint::from(1u8) << 64usize;
            let boundaries = [
                0,
                1,
                p - 1,
                p,
                2 * p,
                p * p,
                p * p - 1,
                (u64::MAX / p) * p,
                u64::MAX - 1,
                u64::MAX,
            ];
            for i in 0..1010 {
                random ^= random << 13;
                random ^= random >> 7;
                random ^= random << 17;
                let total = boundaries.get(i).copied().unwrap_or(random);
                let got = reducer.reduce_bounded(total);
                assert!(got < p);
                assert_eq!(got, reducer.reduce_accumulated_products(total));
                assert_eq!(
                    (BigUint::from(got) * &radix) % &modulus,
                    BigUint::from(total) % &modulus,
                    "p={p}, total={total}"
                );
                let correction = total.wrapping_mul(reducer.negative_inverse);
                let raw =
                    ((u128::from(total) + u128::from(correction) * u128::from(p)) >> 64) as u64;
                assert!(raw <= p);
                assert_eq!(reducer.reduce_bounded_lazy(total), raw);
                assert_eq!(raw == p, total != 0 && total % p == 0);
                let encoded = reducer.encode(total);
                assert!(encoded < p);
                assert_eq!(
                    BigUint::from(encoded),
                    (BigUint::from(total) * &radix) % &modulus
                );
                assert_eq!(reducer.decode(encoded), total % p);
                let a = total % p;
                let b = total.rotate_left(17) % p;
                let product =
                    reducer.decode(reducer.multiply(reducer.encode(a), reducer.encode(b)));
                assert_eq!(
                    BigUint::from(product),
                    (BigUint::from(a) * BigUint::from(b)) % &modulus
                );
            }
            // Even canonical products need normalization for composite moduli.
            assert_eq!(reducer.decode(p), 0);
        }
        let composite = Montgomery32::new(15);
        assert_eq!(composite.multiply(3, 5), 0);
    }

    #[test]
    fn radix32_requires_its_own_product_bound_and_encoding() {
        let r = 1u128 << 32;
        for p in [3u64, 9, 42_447_347, 999_999_893, 999_999_999] {
            let modulus = u128::from(p);
            let inverse = Montgomery32::new(p).negative_inverse as u32;
            let redc = |total: u64| {
                let m = (total as u32).wrapping_mul(inverse);
                (u128::from(total) + u128::from(m) * modulus) >> 32
            };
            let cap = (modulus * r - 1) / (modulus - 1).pow(2);
            assert!(cap * (modulus - 1).pow(2) < modulus * r);
            assert!((cap + 1) * (modulus - 1).pow(2) >= modulus * r);
            for t in [0, p, (modulus * r - 1) as u64] {
                let raw = redc(t);
                assert!(raw < 2 * modulus);
                let canonical = if raw >= modulus { raw - modulus } else { raw };
                assert_eq!(canonical * r % modulus, u128::from(t) % modulus);
            }
            assert_eq!(redc(p), modulus);
            assert!(redc(u64::MAX) >= 2 * modulus);
            let r64 = (r * r) % modulus;
            let wrong_scale = redc((r64 * r64) as u64) % modulus;
            assert_eq!(wrong_scale * r % modulus, r64 * r64 % modulus);
            if p == 999_999_893 {
                assert_eq!(cap, 4);
                assert_ne!(r % modulus, r64);
                assert_ne!(
                    wrong_scale, r64,
                    "R32 cannot multiply R64-encoded ones unchanged"
                );
            }
        }
    }

    #[test]
    fn accumulated_redc_matches_independent_modulo() {
        for p in [3, 5, 103, 42_447_347, 3_999_999_959, 4_294_967_295] {
            let reducer = Montgomery32::new(p);
            let radix = (1u128 << 64) % u128::from(p);
            for total in [0, 1, p - 1, p * p + 1, u64::MAX - 1, u64::MAX] {
                let reduced = reducer.reduce_accumulated_products(total);
                assert!(reduced < p);
                assert_eq!(
                    u128::from(reduced) * radix % u128::from(p),
                    u128::from(total) % u128::from(p),
                    "p={p}, total={total}"
                );
            }
        }
    }

    #[test]
    fn accumulated_encoded_squares_have_one_radix_factor() {
        for p in [3, 5, 103, 42_447_347, 3_999_999_959] {
            let reducer = Montgomery32::new(p);
            let count = (u64::MAX / ((p - 1) * (p - 1))).min(8192);
            let mut total = 0u64;
            let mut expected = 0u64;
            for _ in 0..count {
                let encoded = p - 1;
                total = total.checked_add(encoded * encoded).unwrap();
                let value = reducer.decode(encoded);
                expected = (expected + mul_mod(value, value, p)) % p;
            }
            if count > 2 {
                assert!(total > p * p);
            }
            let encoded_sum = reducer.reduce_accumulated_products(total);
            assert_eq!(encoded_sum, reducer.encode(expected));
            assert_eq!(reducer.decode(encoded_sum), expected);
        }
    }

    #[test]
    fn montgomery_matches_u128_modulo() {
        for modulus in [3, 103, 65_537, 1_000_003, 3_999_999_959] {
            let montgomery = Montgomery32::new(modulus);
            let values = [0, 1, 2, 17, modulus / 2, modulus - 1];
            for a in values {
                for b in values {
                    let expected = mul_mod(a, b, modulus);
                    let got = montgomery
                        .decode(montgomery.multiply(montgomery.encode(a), montgomery.encode(b)));
                    assert_eq!(got, expected);
                }
            }
        }
    }
}
