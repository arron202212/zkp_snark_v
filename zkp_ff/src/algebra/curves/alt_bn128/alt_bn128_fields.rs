/** @file
 *****************************************************************************
 * Declaration of interfaces for the ALT_BN128 finite fields.
 *
 * Rust port of:
 * libff/algebra/curves/alt_bn128/alt_bn128_fields.hpp
 *
 * Copyright (c) SCIPR Lab and contributors.
 * MIT license.
 *****************************************************************************/
use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_extension::fp2::Fp2;
use crate::algebra::fields::prime_extension::fp6_3over2::Fp6;
use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12;

// -----------------------------------------------------------------------------
// Field sizes
// -----------------------------------------------------------------------------

/// Scalar field bit count.
pub const ALT_BN128_R_BITCOUNT: usize = 254;

/// Base field bit count.
pub const ALT_BN128_Q_BITCOUNT: usize = 254;

/// Number of u64 limbs required by Fr.
pub const ALT_BN128_R_LIMBS: usize = (ALT_BN128_R_BITCOUNT + 63) / 64;

/// Number of u64 limbs required by Fq.
pub const ALT_BN128_Q_LIMBS: usize = (ALT_BN128_Q_BITCOUNT + 63) / 64;

// -----------------------------------------------------------------------------
// Field moduli
//
// Limbs are stored in little-endian order.
// These values use the existing Rust u64 Montgomery implementation.
// -----------------------------------------------------------------------------

/// BN254 scalar field modulus:
///
/// 21888242871839275222246405745257275088548364400416034343698204186575808495617
pub const ALT_BN128_MODULUS_R: [u64; ALT_BN128_R_LIMBS] = [
    0x43e1_f593_f000_0001,
    0x2833_e848_79b9_7091,
    0xb850_45b6_8181_585d,
    0x3064_4e72_e131_a029,
];

/// BN254 base field modulus:
///
/// 21888242871839275222246405745257275088696311157297823662689037894645226208583
pub const ALT_BN128_MODULUS_Q: [u64; ALT_BN128_Q_LIMBS] = [
    0x3c20_8c16_d87c_fd47,
    0x9781_6a91_6871_ca8d,
    0xb850_45b6_8181_585d,
    0x3064_4e72_e131_a029,
];

// -----------------------------------------------------------------------------
// Field types
// -----------------------------------------------------------------------------

/// BN254 scalar field Fr.
pub type AltBn128Fr = Fp<ALT_BN128_R_LIMBS>;

/// BN254 base field Fq.
pub type AltBn128Fq = Fp<ALT_BN128_Q_LIMBS>;

/// BN254 quadratic extension field Fq2.
pub type AltBn128Fq2 = Fp2;

/// BN254 sextic extension field Fq6.
pub type AltBn128Fq6 = Fp6;

/// BN254 degree-12 extension field Fq12.
pub type AltBn128Fq12 = Fp12;

/// BN254 target group GT.
pub type AltBn128GT = AltBn128Fq12;

// -----------------------------------------------------------------------------
// Initialization
// -----------------------------------------------------------------------------

/// Initialize all ALT_BN128 finite-field parameters.
///
/// This function initializes:
///
/// - Fr Montgomery parameters
/// - Fq Montgomery parameters
/// - Fq2 extension parameters
/// - Fq6 extension parameters
/// - Fq12 extension parameters
///
/// Implementation is provided by alt_bn128_init.rs.
pub fn init_alt_bn128_fields() {
    crate::algebra::curves::alt_bn128::alt_bn128_init::init_alt_bn128_fields_impl();
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alt_bn128_field_sizes() {
        assert_eq!(ALT_BN128_R_BITCOUNT, 254);
        assert_eq!(ALT_BN128_Q_BITCOUNT, 254);

        assert_eq!(ALT_BN128_R_LIMBS, 4);
        assert_eq!(ALT_BN128_Q_LIMBS, 4);
    }

    #[test]
    fn test_alt_bn128_modulus_r() {
        assert_eq!(
            ALT_BN128_MODULUS_R,
            [
                0x43e1_f593_f000_0001,
                0x2833_e848_79b9_7091,
                0xb850_45b6_8181_585d,
                0x3064_4e72_e131_a029,
            ]
        );
    }

    #[test]
    fn test_alt_bn128_modulus_q() {
        assert_eq!(
            ALT_BN128_MODULUS_Q,
            [
                0x3c20_8c16_d87c_fd47,
                0x9781_6a91_6871_ca8d,
                0xb850_45b6_8181_585d,
                0x3064_4e72_e131_a029,
            ]
        );
    }

    #[test]
    fn test_alt_bn128_moduli_are_distinct() {
        assert_ne!(ALT_BN128_MODULUS_R, ALT_BN128_MODULUS_Q);
    }
}
