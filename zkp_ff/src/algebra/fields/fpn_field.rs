
use crate::algebra::field_utils::bigint::BigInt;
use crate::algebra::fields::field::Field;

/// Common API for prime characteristic finite fields.
///
/// Covers Fp and its extension fields Fp^n.
pub trait FpnField: Field {

    /// Number of limbs in the base modulus.
    const MODULUS_LIMBS: usize;

    /// Base field characteristic.
    fn field_char() -> BigInt<1>;

    /// Euler exponent:
    ///
    /// (q - 1) / 2
    ///
    /// q = p^extension_degree
    fn euler() -> Vec<u64>;

    /// Two-adicity:
    ///
    /// q - 1 = 2^s * t
    fn two_adicity() -> usize;

    /// Odd component of q - 1.
    fn t() -> Vec<u64>;

    /// (t - 1) / 2
    fn t_minus_1_over_2() -> Vec<u64>;

    /// Quadratic non-residue.
    fn nqr() -> Self;

    /// nqr^t
    fn nqr_to_t() -> Self;

    /// Frobenius map.
    ///
    /// Applies the p^power Frobenius automorphism.
    fn frobenius_map(&self, power: usize) -> Self;



    /// Initialize Tonelli-Shanks constants.
    fn init_tonelli_shanks_constants();
}



// pub trait Field {
//     // 通用有限域接口

//     fn square(&mut self);
//     fn invert(&mut self);

//     fn squared(&self) -> Self;
//     fn inverse(&self) -> Self;

//     fn sqrt(&self) -> Self;

//     fn is_zero(&self) -> bool;

//     fn zero() -> Self;
//     fn one() -> Self;

//     fn extension_degree() -> usize;

//     fn to_words(&self) -> Vec<u64>;
//     fn from_words(&mut self, words: &[u64]) -> bool;
// }

// pub trait FpnField: Field {
//     // 素数特征有限域特有接口

//     fn field_char() -> Vec<u64>;

//     fn euler() -> Vec<u64>;
//     fn two_adicity() -> usize;

//     fn t() -> Vec<u64>;
//     fn t_minus_1_over_2() -> Vec<u64>;

//     fn nqr() -> Self;
//     fn nqr_to_t() -> Self;

//     fn frobenius_map(&self, power: usize) -> Self;
// }
