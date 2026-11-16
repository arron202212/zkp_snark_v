
use std::ops::{
    Add, AddAssign,
    Sub, SubAssign,
    Mul, MulAssign,
    Neg,
};

/// Common API for all finite fields.
///
/// Supports both:
/// - Prime characteristic extension fields F_{p^n}
/// - Binary extension fields F_{2^n}
///
/// This trait is an interface contract.
/// It does not prescribe the internal representation.
pub trait Field:
    Clone
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + Sized
{
    /// Addition in place.
    fn add_assign(&mut self, other: &Self);

    /// Subtraction in place.
    fn sub_assign(&mut self, other: &Self);

    /// Multiplication in place.
    fn mul_assign(&mut self, other: &Self);

    /// Exponentiation by u64.
    fn pow_assign(&mut self, exp: u64);

    /// Exponentiation by arbitrary precision integer.
    fn pow_bigint(&mut self, exp: &[u64]);

    /// Square in place.
    fn square(&mut self);

    /// Inverse in place.
    fn invert(&mut self);

    /// Return squared element.
    fn squared(&self) -> Self {
        let mut result = self.clone();
        result.square();
        result
    }

    /// Return inverse element.
    fn inverse(&self) -> Self {
        let mut result = self.clone();
        result.invert();
        result
    }

    /// Square root.
    ///
    /// The input must be a quadratic residue.
    fn sqrt(&self) -> Self;

    /// Exponentiation by u64.
    fn pow(&self, exp: u64) -> Self {
        let mut result = self.clone();
        result.pow_assign(exp);
        result
    }

    /// Exponentiation by arbitrary precision integer.
    fn pow_bigint_ref(&self, exp: &[u64]) -> Self {
        let mut result = self.clone();
        result.pow_bigint(exp);
        result
    }

    /// Check whether this element is zero.
    fn is_zero(&self) -> bool;

    /// Print element.
    fn print(&self);

    /// Convert to little-endian u64 words.
    fn to_words(&self) -> Vec<u64>;

    /// Restore element from little-endian u64 words.
    ///
    /// Returns false if the represented value is invalid.
    fn from_words(&mut self, words: &[u64]) -> bool;

    /// Randomize current element.
    fn randomize(&mut self);

    /// Clear element.
    fn clear(&mut self);

    /// Zero element.
    fn zero() -> Self;

    /// One element.
    fn one() -> Self;

    /// Random element.
    fn random_element() -> Self;

    /// Extension degree.
    ///
    /// Equals 1 for a prime field.
    fn extension_degree() -> usize;

    /// Number of bits rounded up.
    fn ceil_size_in_bits() -> usize;

    /// Number of bits rounded down.
    fn floor_size_in_bits() -> usize;
}
