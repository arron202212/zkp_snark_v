use crate::algebra::curves::alt_bn128::alt_bn128_g1::G1;
use crate::algebra::fields::prime_extension::fp2::Fp2;

// /// Miller Loop 中的一组直线系数。
// ///
// /// 对应 C++:
// /// struct alt_bn128_ate_ell_coeffs
// #[derive(Clone, Debug, PartialEq, Eq)]
// pub struct AteEllCoeffs {
//     pub ell_0: Fp2,
//     pub ell_vw: Fp2,
//     pub ell_vv: Fp2,
// }

// /// G2 配对预计算结果。
// ///
// /// 对应 C++:
// /// struct alt_bn128_ate_G2_precomp
// #[derive(Clone, Debug, PartialEq, Eq)]
// pub struct AteG2Precomp {
//     pub qx: Fp2,
//     pub qy: Fp2,
//     pub coeffs: Vec<AteEllCoeffs>,
// }

// /// G1 预计算
// pub fn ate_precompute_g1(p: &G1) -> AteG1Precomp;

// /// G2 预计算
// pub fn ate_precompute_g2(q: &G2) -> AteG2Precomp;

// /// Miller Loop
// pub fn ate_miller_loop(
//     p: &AteG1Precomp,
//     q: &AteG2Precomp,
// ) -> Fp12;

// /// Double Miller Loop
// pub fn ate_double_miller_loop(
//     p1: &AteG1Precomp,
//     q1: &AteG2Precomp,
//     p2: &AteG1Precomp,
//     q2: &AteG2Precomp,
// ) -> Fp12;

// /// Ate Pairing
// pub fn ate_pairing(p: &G1, q: &G2) -> Fp12;

// /// Final Exponentiation
// pub fn final_exponentiation(elt: &Fp12) -> Option<Fp12>;

// /// Ate Reduced Pairing
// pub fn ate_reduced_pairing(p: &G1, q: &G2) -> Option<Fp12>;

// /// 通用预计算接口
// pub fn precompute_g1(p: &G1) -> AteG1Precomp;

// pub fn precompute_g2(q: &G2) -> AteG2Precomp;

// /// 通用 Miller Loop
// pub fn miller_loop(
//     p: &AteG1Precomp,
//     q: &AteG2Precomp,
// ) -> Fp12;

// /// 通用 Double Miller Loop
// pub fn double_miller_loop(
//     p1: &AteG1Precomp,
//     q1: &AteG2Precomp,
//     p2: &AteG1Precomp,
//     q2: &AteG2Precomp,
// ) -> Fp12;

// /// 通用 Pairing
// pub fn pairing(p: &G1, q: &G2) -> Fp12;

// /// 通用 Reduced Pairing
// pub fn reduced_pairing(p: &G1, q: &G2) -> Option<Fp12>;

use crate::algebra::fields::prime_base::fp::Fp;
// use crate::algebra::fields::prime_extension::fp2::Fp2;
use crate::algebra::fields::prime_extension::fp6_2over3::Fp6;
use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AteG1Precomp {
    pub px: Fp,
    pub py: Fp,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AteEllCoeffs {
    pub ell_0: Fp2,
    pub ell_vw: Fp2,
    pub ell_vv: Fp2,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AteG2Precomp {
    pub qx: Fp2,
    pub qy: Fp2,
    pub coeffs: Vec<AteEllCoeffs>,
}

pub fn final_exponentiation_first_chunk(elt: &Fp12) -> Option<Fp12> {
    let a = Fp12::new(elt.c0.clone(), elt.c1.negate(), elt.params.clone());

    let b = elt.inverse()?;

    let c = a.mul(&b);
    let d = c.frobenius_map(2);

    Some(d.mul(&c))
}
