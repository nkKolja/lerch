//! Independent, intentionally slow references; never used by production search.

use crate::arith::{fermat_quotient_mod_p, mul_mod, pow_mod, prime_base};
use crate::moments::Canonical;
use crate::reference_recurrence::{MomentOptions, recurrence_invariants};

/// Original generic-g O(p) implementation, including its original extra invariants.
pub fn generic(p: u64) -> Result<Canonical, String> {
    let base = prime_base(p)?;
    let result = std::hint::black_box(recurrence_invariants(p, &base, MomentOptions::default()));
    Ok(Canonical {
        p,
        primitive_root: result.primitive_root,
        q1: result.q1,
        q2: result.q2,
        lerch_remainder: result.lerch_remainder,
        is_lerch: result.is_lerch,
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DirectInvariants {
    pub q1: u64,
    pub q2: u64,
    pub q3: Option<u64>,
    pub q4: Option<u64>,
    pub wilson: u64,
    pub lerch_remainder: u64,
}

pub fn direct_values(p: u64) -> Vec<u64> {
    (1..p).map(|a| fermat_quotient_mod_p(a, p)).collect()
}

/// Independent definition-based path: one exponentiation modulo p^2 per a.
pub fn direct_invariants(p: u64, q3_enabled: bool, q4_enabled: bool) -> DirectInvariants {
    assert!((3..=4_000_000_000).contains(&p));
    let p2 = p * p;
    let mut q1 = 0;
    let mut q2 = 0;
    let mut q3 = 0;
    let mut q4 = 0;
    for a in 1..p {
        let q = (pow_mod(a, p - 1, p2) - 1) / p;
        q1 = (q1 + q) % p;
        let square = mul_mod(q, q, p);
        q2 = (q2 + square) % p;
        if q3_enabled {
            q3 = (q3 + mul_mod(square, q, p)) % p;
        }
        if q4_enabled {
            q4 = (q4 + mul_mod(square, square, p)) % p;
        }
    }
    let mut factorial = 1u64;
    for a in 2..p {
        factorial = mul_mod(factorial, a, p2);
    }
    let wilson = ((factorial + 1) / p) % p;
    let rhs = (q2 + mul_mod(wilson, wilson, p) + 2 * (p - wilson)) % p;
    let lerch_remainder = mul_mod(rhs, p.div_ceil(2), p);
    DirectInvariants {
        q1,
        q2,
        q3: q3_enabled.then_some(q3),
        q4: q4_enabled.then_some(q4),
        wilson,
        lerch_remainder,
    }
}
