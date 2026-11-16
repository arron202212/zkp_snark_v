use crate::algebra::field_utils::bigint::BigInt;

/// Binary finite field common API.
///
/// Represents finite fields GF(2^n).
///
/// This trait documents the common interface for binary fields.
/// It does not require a particular internal representation.
pub trait BinaryField:
    Clone
    + Sized
    + PartialEq
{
    /// Binary field modulus.
    ///
    /// Usually represented as an irreducible polynomial.
    const MODULUS: u64;

    /// Number of bits in the field representation.
    const NUM_BITS: usize;

    /// Multiplicative generator.
    fn multiplicative_generator() -> Self;

    /// Characteristic of the field.
    fn field_char<const N: usize>() -> BigInt<N> {
        BigInt::<N>::from(2u64)
    }

    /// Addition.
    fn add_assign(&mut self, other: &Self);

    /// Subtraction.
    fn sub_assign(&mut self, other: &Self);

    /// Multiplication.
    fn mul_assign(&mut self, other: &Self);

    /// Exponentiation.
    fn pow_assign(&mut self, exp: u64);

    /// Square in place.
    fn square(&mut self);

    /// Inverse in place.
    fn invert(&mut self);

    /// Additive inverse.
    fn negate(&self) -> Self;

    /// Square root.
    fn sqrt(&self) -> Self;

    /// Squared value.
    fn squared(&self) -> Self {
        let mut result = self.clone();
        result.square();
        result
    }

    /// Multiplicative inverse.
    fn inverse(&self) -> Self {
        let mut result = self.clone();
        result.invert();
        result
    }

    /// Zero element.
    fn zero() -> Self;

    /// One element.
    fn one() -> Self;

    /// Random element.
    fn random_element() -> Self;

    /// Number of bits.
    fn ceil_size_in_bits() -> usize {
        Self::NUM_BITS
    }

    /// Number of bits.
    fn floor_size_in_bits() -> usize {
        Self::NUM_BITS
    }

    /// Extension degree.
    fn extension_degree() -> usize;

    /// Convert to little-endian u64 words.
    fn to_words(&self) -> Vec<u64>;

    /// Restore from little-endian u64 words.
    fn from_words(words: &[u64]) -> bool;

    /// Randomize current element.
    fn randomize(&mut self);

    /// Clear current element.
    fn clear(&mut self);

    /// Check whether element is zero.
    fn is_zero(&self) -> bool;

    /// Print field element.
    fn print(&self);
}

