use crate::arith::prime_base;
use crate::moments::Canonical;
use crate::reference::{direct_invariants, generic};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VerificationTranscript {
    pub p: u64,
    pub method: String,
    pub definition_checked: bool,
    pub candidate_status: String,
    pub recurrence_q1: u64,
    pub recurrence_q2: u64,
    pub direct_q1_mod_p: u64,
    pub direct_q2_mod_p: u64,
    pub wilson_mod_p: Option<u64>,
    pub recurrence_matches_direct: bool,
    pub q1_matches_wilson_mod_p: Option<bool>,
    pub q1_minus_wilson_mod_p2: Option<String>,
    pub power_sum_minus_factorial_minus_p_mod_p3: Option<String>,
    pub lerch_verified: Option<bool>,
    pub verified: bool,
}

/// The default definition path additionally checks both bigint congruences for
/// a Lerch candidate. A generic-only cross-check never promotes candidate status.
pub fn verify(fast: &Canonical, generic_only: bool) -> Result<VerificationTranscript, String> {
    let p = fast.p;
    prime_base(p)?;
    if p == 2 {
        return Err("the Lerch quotient is defined here only for odd primes".into());
    }
    let (q1, q2, remainder, wilson, root_matches) = if generic_only {
        let reference = generic(p)?;
        (
            reference.q1,
            reference.q2,
            reference.lerch_remainder,
            None,
            reference.primitive_root == fast.primitive_root,
        )
    } else {
        let reference = direct_invariants(p, false, false);
        (
            reference.q1,
            reference.q2,
            Some(reference.lerch_remainder),
            Some(reference.wilson),
            true,
        )
    };
    let recurrence_matches_direct = root_matches
        && q1 == fast.q1
        && q2 == fast.q2
        && remainder == fast.lerch_remainder
        && fast.is_lerch == (remainder == Some(0));
    let q1_matches_wilson_mod_p = wilson.map(|value| q1 == value);
    let (q1_minus_wilson_mod_p2, power_residue, lerch_verified) = if fast.is_lerch && !generic_only
    {
        let full = verify_lerch_bigint(p);
        (
            Some(full.q1_minus_wilson_mod_p2.to_string()),
            Some(full.power_sum_residue.to_string()),
            Some(full.q1_minus_wilson_mod_p2.is_zero() && full.power_sum_residue.is_zero()),
        )
    } else {
        (None, None, None)
    };
    let verified = recurrence_matches_direct
        && q1_matches_wilson_mod_p.unwrap_or(true)
        && lerch_verified.unwrap_or(true);
    Ok(VerificationTranscript {
        p,
        method: if generic_only {
            "original generic-g recurrence; not definition-level verification".into()
        } else if fast.is_lerch {
            "direct q_p definitions plus independent bigint p^3 power sum".into()
        } else {
            "direct q_p definitions and factorial modulo p^2".into()
        },
        definition_checked: !generic_only,
        candidate_status: if !verified {
            "verification-failed"
        } else if !fast.is_lerch {
            "not-a-candidate"
        } else if generic_only {
            "pending-independent-verification"
        } else {
            "verified-lerch"
        }
        .into(),
        recurrence_q1: fast.q1,
        recurrence_q2: fast.q2,
        direct_q1_mod_p: q1,
        direct_q2_mod_p: q2,
        wilson_mod_p: wilson,
        recurrence_matches_direct,
        q1_matches_wilson_mod_p,
        q1_minus_wilson_mod_p2,
        power_sum_minus_factorial_minus_p_mod_p3: power_residue,
        lerch_verified,
        verified,
    })
}

pub struct BigintLerchVerification {
    pub q1_minus_wilson_mod_p2: BigUint,
    pub power_sum_residue: BigUint,
}

pub fn verify_lerch_bigint(p: u64) -> BigintLerchVerification {
    assert!(p >= 3);
    let pb = BigUint::from(p);
    let p2 = &pb * &pb;
    let p3 = &p2 * &pb;
    let exponent = BigUint::from(p - 1);
    let mut q1_mod_p2 = BigUint::zero();
    let mut power_sum = BigUint::zero();
    for a in 1..p {
        let power = BigUint::from(a).modpow(&exponent, &p3);
        let q = (&power - BigUint::one()) / &pb;
        q1_mod_p2 = (q1_mod_p2 + &q) % &p2;
        power_sum = (power_sum + power) % &p3;
    }
    let mut factorial = BigUint::one();
    for a in 2..p {
        factorial = (factorial * a) % &p3;
    }
    let wilson_mod_p2 = ((&factorial + BigUint::one()) / &pb) % &p2;
    let q1_minus_wilson_mod_p2 = (&q1_mod_p2 + &p2 - wilson_mod_p2) % &p2;
    let power_sum_residue = (power_sum + &p3 - factorial + &p3 - &pb) % &p3;
    BigintLerchVerification {
        q1_minus_wilson_mod_p2,
        power_sum_residue,
    }
}

/// Direct Lerch quotient modulo p for validation, independently using p^3.
pub fn direct_lerch_remainder_bigint(p: u64) -> u64 {
    let proof = verify_lerch_bigint(p);
    (proof.q1_minus_wilson_mod_p2 / p)
        .to_u64()
        .expect("residue fits u64")
}
