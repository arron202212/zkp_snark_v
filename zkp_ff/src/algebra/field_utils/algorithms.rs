//! Exponentiation and Tonelli-Shanks utilities.
//!
//! Rust port of libff algorithms.tcc.
//!
//! Exponent limbs are little-endian:
//! exponent[0] contains the least significant 64 bits.

/// Minimal interface required by square-and-multiply.
pub trait FieldArithmetic: Clone {
    fn one() -> Self;

    fn is_zero(&self) -> bool;

    fn mul(&self, rhs: &Self) -> Self;

    /// Optimized squaring operation.
    fn squared(&self) -> Self;
}

/// Square-and-multiply using a little-endian limb array.
///
/// This follows libff::power(bigint<m>):
/// scan from the most significant set bit to the least significant bit.
pub fn power<F: FieldArithmetic>(base: &F, exponent: &[u64]) -> F {
    let mut result = F::one();
    let mut found_one = false;

    let max_bits = exponent
        .iter()
        .rposition(|&limb| limb != 0)
        .map(|i| i * 64 + (64 - exponent[i].leading_zeros() as usize))
        .unwrap_or(0);

    for bit_index in (0..max_bits).rev() {
        if found_one {
            result = result.squared();
        }

        if get_bit(exponent, bit_index) {
            found_one = true;
            result = result.mul(base);
        }
    }

    result
}

/// Square-and-multiply for an unsigned 64-bit exponent.
pub fn power_u64<F: FieldArithmetic>(base: &F, exponent: u64) -> F {
    power(base, &[exponent])
}

/// Square-and-multiply for an unsigned 128-bit exponent.
pub fn power_u128<F: FieldArithmetic>(base: &F, exponent: u128) -> F {
    power(base, &[exponent as u64, (exponent >> 64) as u64])
}

/// Square-and-multiply for a vector of 64-bit limbs.
///
/// The vector is interpreted as a little-endian integer:
/// exponent[0] is the least significant limb.
pub fn power_words<F: FieldArithmetic>(base: &F, exponent: &[u64]) -> F {
    power(base, exponent)
}

#[inline]
fn get_bit(words: &[u64], bit_index: usize) -> bool {
    let word_index = bit_index / 64;
    let bit_in_word = bit_index % 64;

    words
        .get(word_index)
        .map(|word| ((word >> bit_in_word) & 1) != 0)
        .unwrap_or(false)
}

/// Parameters required by the Tonelli-Shanks square-root algorithm.
#[derive(Clone, Debug)]
pub struct TonelliShanksParams<F> {
    /// Two-adicity of the multiplicative group order:
    /// p - 1 = 2^s * t, with t odd.
    pub s: usize,

    /// (t - 1) / 2, represented as little-endian u64 limbs.
    pub t_minus_1_over_2: Vec<u64>,

    /// A precomputed value nqr^t, where nqr is a quadratic non-residue.
    pub nqr_to_t: F,

    /// The odd component t, used to validate the parameters.
    pub t: Vec<u64>,
}

/// Errors returned by the safe Tonelli-Shanks implementation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TonelliShanksError {
    InvalidTwoAdicity,
    TIsNotOdd,
    ZeroNonResiduePower,
    NotASquare,
    InvalidIteration,
}

/// Interface required for the Tonelli-Shanks implementation.
///
/// `pow_words()` must interpret the exponent as little-endian u64 limbs.
pub trait TonelliShanksField: FieldArithmetic + PartialEq {
    fn zero() -> Self;

    fn pow_words(&self, exponent: &[u64]) -> Self;
}

/// Tonelli-Shanks square root.
///
/// This follows the structure of libff::cyclotomic_squared(), which
/// computes a square root using the precomputed field parameters.
///
/// Note: despite its C++ function name, this routine computes a square
/// root; it does not compute the square of its input.
pub fn cyclotomic_squared<F: TonelliShanksField>(
    value: &F,
    params: &TonelliShanksParams<F>,
) -> Result<F, TonelliShanksError> {
    if params.s == 0 {
        return Err(TonelliShanksError::InvalidTwoAdicity);
    }

    if params.t.is_empty() || params.t[0] & 1 == 0 {
        return Err(TonelliShanksError::TIsNotOdd);
    }

    if params.nqr_to_t.is_zero() {
        return Err(TonelliShanksError::ZeroNonResiduePower);
    }

    if value.is_zero() {
        return Ok(F::zero());
    }

    let one = F::one();

    let mut v = params.s;
    let mut z = params.nqr_to_t.clone();

    // w = value^((t - 1) / 2)
    let mut w = value.pow_words(&params.t_minus_1_over_2);

    // x = value * w
    let mut x = value.mul(&w);

    // b = value^t
    let mut b = x.mul(&w);

    // Euler criterion:
    // b^(2^(s-1)) must equal one for a quadratic residue.
    //
    // libff performs this check only in DEBUG builds. Rust returns
    // NotASquare instead of risking an infinite loop.
    let mut check = b.clone();

    for _ in 0..(v - 1) {
        check = check.squared();
    }

    if check != one {
        return Err(TonelliShanksError::NotASquare);
    }

    // Tonelli-Shanks iteration.
    while b != one {
        let mut m = 0usize;
        let mut b2m = b.clone();

        while b2m != one {
            b2m = b2m.squared();
            m += 1;

            if m >= v {
                return Err(TonelliShanksError::InvalidIteration);
            }
        }

        // j = v - m - 1
        let j = v
            .checked_sub(m + 1)
            .ok_or(TonelliShanksError::InvalidIteration)?;

        // w = z^(2^j)
        w = z.clone();

        for _ in 0..j {
            w = w.squared();
        }

        // z = w^2
        z = w.squared();

        // b = b * z
        b = b.mul(&z);

        // x = x * w
        x = x.mul(&w);

        v = m;
    }

    // Defensive verification of the result.
    if x.squared() != *value {
        return Err(TonelliShanksError::InvalidIteration);
    }

    Ok(x)
}
