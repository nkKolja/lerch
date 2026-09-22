use crate::reduction::Montgomery32;
use serde::{Deserialize, Serialize};

pub const MAX_PRIME: u64 = 2_000_000_000;
pub const DEFAULT_BATCH_SIZE: u64 = 8192;

/// The on-disk tuple and hash encoding are shared with the archived v1 results.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Canonical {
    pub p: u64,
    pub primitive_root: u64,
    pub q1: u64,
    pub q2: u64,
    pub lerch_remainder: Option<u64>,
    pub is_lerch: bool,
}

impl Canonical {
    pub(crate) fn two() -> Self {
        Self {
            p: 2,
            primitive_root: 1,
            q1: 0,
            q2: 0,
            lerch_remainder: None,
            is_lerch: false,
        }
    }
}

/// Ordinary c, and canonical Montgomery R=2^64 encodings of its inverse and quotient.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Segment {
    pub c: u64,
    pub v: u64,
    pub u: u64,
    pub len: u64,
}

/// Encoded sums of X=2U+V over sign pairs; square is not yet a partial Q2.
#[derive(Clone, Copy, Default)]
pub(crate) struct PairedSums {
    pub first: u64,
    pub square: u64,
    pub pairs: u64,
}

#[inline(always)]
pub(crate) fn add(a: u64, b: u64, p: u64) -> u64 {
    let sum = a + b;
    if sum >= p { sum - p } else { sum }
}

#[inline(always)]
pub(crate) fn half(a: u64, p: u64) -> u64 {
    (a + (a & 1) * p) >> 1
}

pub fn square_batch_size(p: u64, requested: u64) -> Result<u64, String> {
    if !(3..=MAX_PRIME).contains(&p) {
        return Err(format!("odd prime must be in 3..={MAX_PRIME}"));
    }
    if requested == 0 {
        return Err("batch size must be positive".into());
    }
    // At the upper limit this cap is four, including when requested is u64::MAX.
    Ok(requested.min(u64::MAX / (p - 1).pow(2)))
}

pub(crate) fn finish(
    p: u64,
    primitive_root: u64,
    reducer: Montgomery32,
    sums: PairedSums,
) -> Result<Canonical, String> {
    if 2 * sums.pairs != p - 1 {
        return Err("doubling cycles did not cover exactly (p-1)/2 sign pairs".into());
    }
    // The half-domain inverse-square sum vanishes for p>3, but equals one for p=3.
    let square = if p == 3 {
        add(sums.square, reducer.encode(1), p)
    } else {
        sums.square
    };
    let q2 = half(square, p);
    let twice_q1 = add(sums.first, sums.first, p);
    let negative_twice_q1 = if twice_q1 == 0 { 0 } else { p - twice_q1 };
    let expression = add(
        add(q2, reducer.multiply(sums.first, sums.first), p),
        negative_twice_q1,
        p,
    );
    Ok(Canonical {
        p,
        primitive_root,
        q1: reducer.decode(sums.first),
        q2: reducer.decode(q2),
        lerch_remainder: Some(reducer.decode(half(expression, p))),
        is_lerch: expression == 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_caps_fit_and_are_maximal() {
        for p in [
            3,
            42_447_347,
            1_000_000_007,
            1_499_999_957,
            1_999_999_973,
            MAX_PRIME,
        ] {
            let cap = square_batch_size(p, u64::MAX).unwrap();
            let square = u128::from(p - 1).pow(2);
            assert!(u128::from(cap) * square <= u128::from(u64::MAX));
            assert!(u128::from(cap + 1) * square > u128::from(u64::MAX));
        }
        assert_eq!(square_batch_size(MAX_PRIME, DEFAULT_BATCH_SIZE), Ok(4));
        assert!(square_batch_size(3, 0).is_err());
        assert!(square_batch_size(MAX_PRIME + 1, 1).is_err());
    }
}
