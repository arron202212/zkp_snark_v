use super::mnt4_g1::Mnt4G1;
use super::mnt4_g2::Mnt4G2;
use super::mnt4_init::{Mnt4InitializedParams, init_mnt4_params};
use super::mnt4_pairing::{
    Mnt4AteG1Precomp, Mnt4AteG2Precomp, Mnt4Fq4, Mnt4GT, mnt4_ate_double_miller_loop,
    mnt4_ate_miller_loop, mnt4_ate_precompute_g1, mnt4_ate_precompute_g2,
    mnt4_final_exponentiation, mnt4_pairing, mnt4_reduced_pairing,
};

/// MNT4 public parameters interface.
///
/// Corresponds to libff::mnt4_pp.
pub struct Mnt4PP;

/// Field and group aliases.
pub type FpType = crate::algebra::fields::prime_base::fp::Fp<5>;
pub type G1Type = Mnt4G1;
pub type G2Type = Mnt4G2;
pub type G1PrecompType = Mnt4AteG1Precomp;
pub type G2PrecompType = Mnt4AteG2Precomp;
pub type FqType = crate::algebra::fields::prime_base::fp::Fp<5>;
pub type FqeType = crate::algebra::fields::prime_extension::fp2::Fp2<5>;
pub type FqkType = Mnt4Fq4;
pub type GTType = Mnt4GT;

/// MNT4 supports affine pairing in the original libff implementation.
///
/// The affine Ate implementation is not yet exposed by the current
/// Rust mnt4_pairing module.
pub const HAS_AFFINE_PAIRING: bool = true;

impl Mnt4PP {
    /// Initialize MNT4 public parameters.
    ///
    /// C++:
    ///     mnt4_pp::init_public_params()
    ///
    /// Rust returns the initialized parameters instead of relying
    /// on C++-style global mutable state.
    pub fn init_public_params() -> Mnt4InitializedParams {
        init_mnt4_params()
    }

    /// Final exponentiation.
    ///
    /// C++:
    ///     mnt4_pp::final_exponentiation(elt)
    pub fn final_exponentiation(elt: &Mnt4Fq4) -> Mnt4GT {
        mnt4_final_exponentiation(elt)
    }

    /// Precompute G1.
    ///
    /// C++:
    ///     mnt4_pp::precompute_G1(P)
    pub fn precompute_g1(p: &G1Type) -> G1PrecompType {
        mnt4_ate_precompute_g1(p)
    }

    /// Precompute G2.
    ///
    /// C++:
    ///     mnt4_pp::precompute_G2(Q)
    pub fn precompute_g2(q: &G2Type) -> G2PrecompType {
        mnt4_ate_precompute_g2(q)
    }

    /// Miller loop using precomputed G1 and G2.
    ///
    /// C++:
    ///     mnt4_pp::miller_loop(prec_P, prec_Q)
    pub fn miller_loop(prec_p: &G1PrecompType, prec_q: &G2PrecompType) -> Mnt4Fq4 {
        mnt4_ate_miller_loop(prec_p, prec_q)
    }

    /// Double Miller loop.
    ///
    /// C++:
    ///     mnt4_pp::double_miller_loop(...)
    pub fn double_miller_loop(
        prec_p1: &G1PrecompType,
        prec_q1: &G2PrecompType,
        prec_p2: &G1PrecompType,
        prec_q2: &G2PrecompType,
    ) -> Mnt4Fq4 {
        mnt4_ate_double_miller_loop(prec_p1, prec_q1, prec_p2, prec_q2)
    }

    /// Pairing.
    ///
    /// C++:
    ///     mnt4_pp::pairing(P, Q)
    pub fn pairing(p: &G1Type, q: &G2Type) -> Mnt4Fq4 {
        mnt4_pairing(p, q)
    }

    /// Reduced pairing.
    ///
    /// C++:
    ///     mnt4_pp::reduced_pairing(P, Q)
    pub fn reduced_pairing(p: &G1Type, q: &G2Type) -> Mnt4Fq4 {
        mnt4_reduced_pairing(p, q)
    }

    /// Ate Miller loop using precomputed affine Ate data.
    ///
    /// This is the low-level Ate Miller loop already exposed
    /// by the current Rust pairing module.
    pub fn ate_miller_loop(prec_p: &G1PrecompType, prec_q: &G2PrecompType) -> Mnt4Fq4 {
        mnt4_ate_miller_loop(prec_p, prec_q)
    }
}

// use crate::algebra::curves::mnt::mnt4::mnt4_pp::Mnt4PP;

// let params = Mnt4PP::init_public_params();

// let p = Mnt4PP::precompute_g1(&g1);
// let q = Mnt4PP::precompute_g2(&g2);

// let miller = Mnt4PP::miller_loop(&p, &q);

// let result = Mnt4PP::final_exponentiation(&miller);
