//! Generic field utilities corresponding to libff/field_utils.tcc.
//!
//! Bit order and word order follow libff:
//! - bits are packed least-significant-bit first;
//! - words are stored little-endian;
//! - integer packing uses FieldT::floor_size_in_bits() bits per field element.

use core::fmt::Debug;

/// Identifies the category of a field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldType {
    Multiplicative,
    Additive,
}

/// The common interface required by field utility functions.
///
/// Implement this trait for Fp, Fp2, Fp3, Fp6, Fp12, etc.
/// `canonical_words()` must return the canonical (non-Montgomery)
/// representation, in little-endian u64 word order.
///
/// `from_canonical_words()` must interpret its input as a canonical
/// integer representation, not as Montgomery limbs.
pub trait FieldUtils: Clone + PartialEq + Debug {
    const FIELD_TYPE: FieldType;

    /// Number of bits used by the packed representation.
    const FLOOR_BITS: usize;

    /// Number of bits needed to represent an element.
    const CEIL_BITS: usize;

    /// For additive fields, this is the extension degree.
    const EXTENSION_DEGREE: usize;

    /// Two-adicity of the multiplicative group.
    const TWO_ADICITY: usize;

    fn zero() -> Self;
    fn one() -> Self;

    fn is_zero(&self) -> bool;

    fn add(&self, rhs: &Self) -> Self;
    fn sub(&self, rhs: &Self) -> Self;
    fn mul(&self, rhs: &Self) -> Self;

    /// Returns None if the element is not invertible.
    fn inverse(&self) -> Option<Self>;

    /// Canonical little-endian words.
    fn canonical_words(&self) -> Vec<u64>;

    /// Construct from canonical little-endian words.
    fn from_canonical_words(words: &[u64]) -> Option<Self>;

    /// Multiplicative generator. Used by coset_shift().
    fn multiplicative_generator() -> Self;

    /// A primitive 2^TWO_ADICITY-th root of unity.
    fn root_of_unity() -> Self;
}

/// Equivalent to libff::get_field_type().
pub fn get_field_type<F: FieldUtils>(_elem: &F) -> FieldType {
    F::FIELD_TYPE
}

/// Equivalent to libff::log_of_field_size_helper().
pub fn log_of_field_size<F: FieldUtils>(_elem: &F) -> usize {
    match F::FIELD_TYPE {
        FieldType::Multiplicative => F::CEIL_BITS,
        FieldType::Additive => F::EXTENSION_DEGREE,
    }
}

/// Equivalent to libff::soundness_log_of_field_size_helper().
pub fn soundness_log_of_field_size<F: FieldUtils>(_elem: &F) -> usize {
    match F::FIELD_TYPE {
        FieldType::Multiplicative => {
            assert!(F::CEIL_BITS > 0);
            F::CEIL_BITS - 1
        }
        FieldType::Additive => F::EXTENSION_DEGREE,
    }
}

/// Equivalent to libff::get_word_of_field_elem().
///
/// For multiplicative fields, libff reads the bigint representation.
/// For additive fields, libff reads to_words().
///
/// In this Rust interface, canonical_words() provides the corresponding
/// canonical representation for both categories.
pub fn get_word_of_field_elem<F: FieldUtils>(field_elem: &F, word: usize) -> u64 {
    field_elem.canonical_words()[word]
}

/// Equivalent to libff::coset_shift().
pub fn coset_shift<F: FieldUtils>() -> F {
    let generator = F::multiplicative_generator();
    generator.mul(&generator)
}

/// Equivalent to libff::get_root_of_unity() for finite fields.
///
/// `n` must be a power of two and log2(n) must not exceed the field's
/// two-adicity.
pub fn get_root_of_unity<F: FieldUtils>(n: usize) -> Result<F, &'static str> {
    if n == 0 || !n.is_power_of_two() {
        return Err("libff::get_root_of_unity: expected n to be a non-zero power of two");
    }

    let logn = n.trailing_zeros() as usize;

    if logn > F::TWO_ADICITY {
        return Err("libff::get_root_of_unity: expected logn <= FieldT::s");
    }

    let mut omega = F::root_of_unity();

    for _ in logn..F::TWO_ADICITY {
        omega = omega.mul(&omega);
    }

    Ok(omega)
}

/// Equivalent to libff::pack_int_vector_into_field_element_vector().
///
/// Each input usize contributes its lowest `w` bits.
/// Bits are packed little-endian into field elements, with
/// FLOOR_BITS bits per output element.
pub fn pack_int_vector_into_field_element_vector<F: FieldUtils>(
    values: &[usize],
    w: usize,
) -> Result<Vec<F>, &'static str> {
    if w == 0 {
        return Err("w must be greater than zero");
    }

    if w > usize::BITS as usize {
        return Err("w exceeds the bit width of usize");
    }

    if F::FLOOR_BITS == 0 {
        return Err("field packing capacity must be greater than zero");
    }

    let total_bits = values
        .len()
        .checked_mul(w)
        .ok_or("input bit length overflow")?;

    let chunk_bits = F::FLOOR_BITS;
    let output_len = div_ceil(total_bits, chunk_bits);

    let mut result = Vec::with_capacity(output_len);

    for i in 0..output_len {
        let mut words = vec![0u64; div_ceil(chunk_bits, 64)];

        for j in 0..chunk_bits {
            let global_bit = i * chunk_bits + j;

            if global_bit >= total_bits {
                break;
            }

            let word_index = global_bit / w;
            let bit_in_word = global_bit % w;
            let bit = (values[word_index] >> bit_in_word) & 1;

            if bit != 0 {
                set_bit(&mut words, j);
            }
        }

        let elem = F::from_canonical_words(&words)
            .ok_or("packed integer is not a valid canonical field element")?;

        result.push(elem);
    }

    Ok(result)
}

/// Equivalent to libff::pack_bit_vector_into_field_element_vector(v, chunk_bits).
pub fn pack_bit_vector_into_field_element_vector<F: FieldUtils>(
    bits: &[bool],
    chunk_bits: usize,
) -> Result<Vec<F>, &'static str> {
    if chunk_bits > F::FLOOR_BITS {
        return Err("chunk_bits exceeds field packing capacity");
    }

    if chunk_bits == 0 {
        return Err("chunk_bits must be greater than zero");
    }

    let output_len = div_ceil(bits.len(), chunk_bits);
    let mut result = Vec::with_capacity(output_len);

    for i in 0..output_len {
        let mut words = vec![0u64; div_ceil(chunk_bits, 64)];

        for j in 0..chunk_bits {
            let bit_index = i * chunk_bits + j;

            if bit_index < bits.len() && bits[bit_index] {
                set_bit(&mut words, j);
            }
        }

        let elem = F::from_canonical_words(&words)
            .ok_or("packed bits are not a valid canonical field element")?;

        result.push(elem);
    }

    Ok(result)
}

/// Equivalent to the overload using FieldT::floor_size_in_bits().
pub fn pack_bit_vector<F: FieldUtils>(bits: &[bool]) -> Result<Vec<F>, &'static str> {
    pack_bit_vector_into_field_element_vector::<F>(bits, F::FLOOR_BITS)
}

/// Equivalent to libff::convert_bit_vector_to_field_element_vector().
///
/// Each bool becomes a separate field element: false -> zero, true -> one.
pub fn convert_bit_vector_to_field_element_vector<F: FieldUtils>(bits: &[bool]) -> Vec<F> {
    bits.iter()
        .map(|&bit| if bit { F::one() } else { F::zero() })
        .collect()
}

/// Equivalent to libff::convert_field_element_to_bit_vector().
///
/// Returns CEIL_BITS bits in least-significant-bit-first order.
pub fn convert_field_element_to_bit_vector<F: FieldUtils>(elem: &F) -> Vec<bool> {
    let words = elem.canonical_words();
    let mut result = Vec::with_capacity(F::CEIL_BITS);

    for bit_index in 0..F::CEIL_BITS {
        let word_index = bit_index / 64;
        let bit_in_word = bit_index % 64;

        let bit = words
            .get(word_index)
            .map(|word| ((word >> bit_in_word) & 1) != 0)
            .unwrap_or(false);

        result.push(bit);
    }

    result
}

/// Equivalent to the overload with bitcount.
///
/// As in the C++ implementation, resize truncates excess bits or pads
/// with false values.
pub fn convert_field_element_to_bit_vector_with_count<F: FieldUtils>(
    elem: &F,
    bitcount: usize,
) -> Vec<bool> {
    let mut result = convert_field_element_to_bit_vector(elem);
    result.resize(bitcount, false);
    result
}

/// Equivalent to libff::convert_field_element_vector_to_bit_vector().
pub fn convert_field_element_vector_to_bit_vector<F: FieldUtils>(values: &[F]) -> Vec<bool> {
    let mut result = Vec::with_capacity(values.len() * F::CEIL_BITS);

    for elem in values {
        result.extend(convert_field_element_to_bit_vector(elem));
    }

    result
}

/// Equivalent to libff::convert_bit_vector_to_field_element().
///
/// The first bit is the least-significant bit:
/// value = b[0] + 2*b[1] + 4*b[2] + ...
pub fn convert_bit_vector_to_field_element<F: FieldUtils>(
    bits: &[bool],
) -> Result<F, &'static str> {
    if bits.len() > F::CEIL_BITS {
        return Err("bit vector exceeds field element bit capacity");
    }

    let mut result = F::zero();
    let mut coefficient = F::one();

    for &bit in bits {
        if bit {
            result = result.add(&coefficient);
        }

        coefficient = coefficient.add(&coefficient);
    }

    Ok(result)
}

/// Equivalent to libff::batch_invert().
///
/// All elements must be nonzero. The function uses one field inversion
/// and O(n) multiplications.
///
/// An empty vector is accepted and remains empty.
pub fn batch_invert<F: FieldUtils>(values: &mut [F]) {
    if values.is_empty() {
        return;
    }

    let mut prefix_products = Vec::with_capacity(values.len());
    let mut acc = F::one();

    for elem in values.iter() {
        assert!(
            !elem.is_zero(),
            "batch_invert requires every element to be nonzero"
        );

        prefix_products.push(acc.clone());
        acc = acc.mul(elem);
    }

    let mut acc_inverse = acc
        .inverse()
        .expect("product of nonzero field elements must be invertible");

    for i in (0..values.len()).rev() {
        let old_elem = values[i].clone();

        values[i] = acc_inverse.mul(&prefix_products[i]);
        acc_inverse = acc_inverse.mul(&old_elem);
    }
}

/// Integer ceiling division.
#[inline]
fn div_ceil(a: usize, b: usize) -> usize {
    debug_assert!(b > 0);

    a / b + usize::from(a % b != 0)
}

/// Set a bit in a little-endian u64 word vector.
#[inline]
fn set_bit(words: &mut [u64], bit_index: usize) {
    let word_index = bit_index / 64;
    let bit_in_word = bit_index % 64;

    words[word_index] |= 1u64 << bit_in_word;
}
