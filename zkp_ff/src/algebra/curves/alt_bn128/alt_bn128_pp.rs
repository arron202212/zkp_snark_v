use crate::algebra::curves::alt_bn128::alt_bn128_g1::G1;
use crate::algebra::curves::alt_bn128::alt_bn128_g2::G2;
use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_extension::fp2::Fp2;
use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12;

use super::alt_bn128_init::init_alt_bn128_params;

use super::alt_bn128_pairing::{
    AltBn128G1Precomp, AltBn128G2Precomp, alt_bn128_ate_miller_loop, alt_bn128_ate_precompute_g1,
    alt_bn128_ate_precompute_g2, alt_bn128_final_exponentiation, alt_bn128_reduced_pairing,
};

/// BN254 / alt_bn128 scalar field.
pub type AltBn128Fr = Fp<4>;

/// BN254 / alt_bn128 base field.
pub type AltBn128Fq = Fp<4>;

/// BN254 quadratic extension field.
pub type AltBn128Fqe = Fp2<4>;

/// BN254 degree-12 extension field.
pub type AltBn128Fqk = Fp12;

/// BN254 target group.
pub type AltBn128GT = Fp12;

/// BN254 pairing parameter set.
///
/// This corresponds to libff's `alt_bn128_pp`.
pub struct AltBn128PP;

impl AltBn128PP {
    /// Whether this curve supports affine pairing.
    ///
    /// Matches libff:
    /// `static const bool has_affine_pairing = false`.
    pub const HAS_AFFINE_PAIRING: bool = false;

    /// Initialize all BN254 public parameters.
    ///
    /// This must be called before using curve and pairing operations
    /// that depend on initialized public parameters.
    pub fn init_public_params() {
        init_alt_bn128_params();
    }

    /// Final exponentiation.
    ///
    /// Converts an Fq12 Miller-loop result into an element of GT.
    pub fn final_exponentiation(elt: &AltBn128Fqk) -> AltBn128GT {
        alt_bn128_final_exponentiation(elt)
    }

    /// Precompute a G1 point for the Miller loop.
    pub fn precompute_g1(p: &AltBn128G1) -> AltBn128G1Precomp {
        alt_bn128_ate_precompute_g1(p)
    }

    /// Precompute a G2 point for the Miller loop.
    pub fn precompute_g2(q: &AltBn128G2) -> AltBn128G2Precomp {
        alt_bn128_ate_precompute_g2(q)
    }

    /// Miller loop using precomputed G1 and G2 points.
    ///
    /// The result is in Fq12 and has not undergone final exponentiation.
    pub fn miller_loop(prec_p: &AltBn128G1Precomp, prec_q: &AltBn128G2Precomp) -> AltBn128Fqk {
        alt_bn128_ate_miller_loop(prec_p, prec_q)
    }

    /// Double Miller loop.
    ///
    /// Computes:
    ///
    ///     Miller(P1, Q1) * Miller(P2, Q2)
    ///
    /// This implementation reuses the existing single Miller-loop
    /// implementation. It is mathematically equivalent to the
    /// simultaneous double Miller loop, but is not as optimized.
    pub fn double_miller_loop(
        prec_p1: &AltBn128G1Precomp,
        prec_q1: &AltBn128G2Precomp,
        prec_p2: &AltBn128G1Precomp,
        prec_q2: &AltBn128G2Precomp,
    ) -> AltBn128Fqk {
        let f1 = Self::miller_loop(prec_p1, prec_q1);
        let f2 = Self::miller_loop(prec_p2, prec_q2);

        f1.mul(&f2)
    }

    /// Raw pairing.
    ///
    /// Returns the Miller-loop result without final exponentiation.
    pub fn pairing(p: &AltBn128G1, q: &AltBn128G2) -> AltBn128Fqk {
        let prec_p = Self::precompute_g1(p);
        let prec_q = Self::precompute_g2(q);

        Self::miller_loop(&prec_p, &prec_q)
    }

    /// Reduced pairing.
    ///
    /// Returns the result after final exponentiation, in GT.
    pub fn reduced_pairing(p: &AltBn128G1, q: &AltBn128G2) -> AltBn128GT {
        alt_bn128_reduced_pairing(p, q)
    }
}

/// C++-style compatibility alias.
pub type AltBn128G1 = G1;

/// C++-style compatibility alias.
pub type AltBn128G2 = G2;

/// C++-style compatibility aliases for precomputation.
pub type AltBn128G1Precomputation = AltBn128G1Precomp;
pub type AltBn128G2Precomputation = AltBn128G2Precomp;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alt_bn128_pp_affine_pairing_flag() {
        assert!(!AltBn128PP::HAS_AFFINE_PAIRING);
    }

    #[test]
    fn test_alt_bn128_pp_type_aliases() {
        let _: Option<AltBn128Fr> = None;
        let _: Option<AltBn128Fq> = None;
        let _: Option<AltBn128Fqe> = None;
        let _: Option<AltBn128Fqk> = None;
        let _: Option<AltBn128GT> = None;
    }
}
