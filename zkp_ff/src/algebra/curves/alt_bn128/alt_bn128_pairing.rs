use crate::algebra::curves::alt_bn128::alt_bn128_g1::{
    // bn254_g1_params,
    G1,
};
use crate::algebra::curves::alt_bn128::alt_bn128_g2::{
    // bn254_g2_params,
    G2,
    // G2Params,
};
use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};
use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12;

use super::alt_bn128_init::{
    ALT_BN128_ATE_IS_LOOP_COUNT_NEG, ALT_BN128_ATE_LOOP_COUNT, ALT_BN128_FINAL_EXPONENT_IS_Z_NEG,
    ALT_BN128_FINAL_EXPONENT_Z, bn254_fp12_params,
};

/// G1 pairing precomputation.
///
/// G1 is represented in affine coordinates.
#[derive(Clone)]
pub struct AltBn128G1Precomp {
    pub px: Fp<4>,
    pub py: Fp<4>,
}

/// One line coefficient in the Ate Miller loop.
#[derive(Clone)]
pub struct AltBn128AteEllCoeffs {
    pub ell_0: Fp2<4>,
    pub ell_vw: Fp2<4>,
    pub ell_vv: Fp2<4>,
}

/// G2 pairing precomputation.
#[derive(Clone)]
pub struct AltBn128G2Precomp {
    pub qx: Fp2<4>,
    pub qy: Fp2<4>,
    pub coeffs: Vec<AltBn128AteEllCoeffs>,
}

/// Return the bit at `index` in a little-endian limb array.
fn get_bit(value: &[u64; 4], index: usize) -> bool {
    let limb = index / 64;
    let bit = index % 64;

    ((value[limb] >> bit) & 1) == 1
}

/// Return the highest set bit.
fn highest_bit(value: &[u64; 4]) -> Option<usize> {
    for i in (0..4).rev() {
        if value[i] != 0 {
            return Some(i * 64 + (63 - value[i].leading_zeros() as usize));
        }
    }

    None
}

/// Construct the BN254 twist non-residue:
///
/// xi = 9 + u
///
/// where u² = -1.
fn bn254_twist() -> Fp2<4> {
    let fq_params = bn254_fq_params();
    let fq2_params = bn254_fq2_params();

    Fp2::new(
        Fp::from_u64(9, fq_params.clone()),
        Fp::one(fq_params),
        fq2_params,
    )
}

fn bn254_fq_params() -> FieldParams<4> {
    FieldParams::new([
        0x3c20_8c16_d87c_fd47,
        0x9781_6a91_6871_ca8d,
        0xb850_45b6_8181_585d,
        0x3064_4e72_e131_a029,
    ])
}

fn bn254_fq2_params() -> Fp2Params<4> {
    let fq = bn254_fq_params();
    let one = Fp::<4>::one(fq.clone());
    let minus_one = one.negate();

    Fp2Params::new(fq, minus_one.clone(), [one, minus_one])
}

/// G1 precomputation.
pub fn alt_bn128_ate_precompute_g1(p: &G1) -> AltBn128G1Precomp {
    let (px, py) = p
        .affine_coordinates()
        .expect("G1 point at infinity has no affine coordinates");

    AltBn128G1Precomp { px, py }
}

/// Perform one flipped Miller doubling step.
///
/// The formulas follow libff's alt_bn128 pairing implementation.
fn doubling_step_for_flipped_miller_loop(
    current: &mut G2,
    two_inv: &Fp2<4>,
    twist: &Fp2<4>,
) -> AltBn128AteEllCoeffs {
    let x = current.x.clone();
    let y = current.y.clone();
    let z = current.z.clone();

    let a = two_inv.mul(&x.mul(&y));
    let b = y.square();
    let c = z.square();

    let three = Fp2::from_u64(3, current.params.fp2_params.clone());

    let d = three.mul(&c);
    let e = twist.mul(&d);
    let f = three.mul(&e);

    let g = two_inv.mul(&b.add(&f));

    let h = y.add(&z).square().sub(&b).sub(&c);

    let i = e.sub(&b);
    let j = x.square();

    let e_squared = e.square();

    current.x = a.mul(&b.sub(&f));
    current.y = g.square().sub(&three.mul(&e_squared));
    current.z = b.mul(&h);

    AltBn128AteEllCoeffs {
        ell_0: twist.mul(&i),
        ell_vw: h.negate(),
        ell_vv: three.mul(&j),
    }
}

/// Perform one flipped Miller mixed-addition step.
fn mixed_addition_step_for_flipped_miller_loop(
    current: &mut G2,
    base: &G2,
    twist: &Fp2<4>,
) -> AltBn128AteEllCoeffs {
    let x1 = current.x.clone();
    let y1 = current.y.clone();
    let z1 = current.z.clone();

    let x2 = base.x.clone();
    let y2 = base.y.clone();

    let d = x1.sub(&x2.mul(&z1));
    let e = y1.sub(&y2.mul(&z1));

    let f = d.square();
    let g = e.square();
    let h = d.mul(&f);
    let i = x1.mul(&f);

    let two = Fp2::from_u64(2, current.params.fp2_params.clone());

    let j = h.add(&z1.mul(&g)).sub(&two.mul(&i));

    current.x = d.mul(&j);
    current.y = e.mul(&i.sub(&j)).sub(&h.mul(&y1));
    current.z = z1.mul(&h);

    AltBn128AteEllCoeffs {
        ell_0: twist.mul(&e.mul(&x2).sub(&d.mul(&y2))),
        ell_vw: d.clone(),
        ell_vv: e.negate(),
    }
}

/// G2 precomputation for the Ate Miller loop.
pub fn alt_bn128_ate_precompute_g2(q: &G2) -> AltBn128G2Precomp {
    let q_affine = q
        .to_affine()
        .expect("G2 point at infinity has no affine coordinates");

    let qx = q_affine.x.clone();
    let qy = q_affine.y.clone();

    let mut current = q_affine.clone();

    let fq2_params = current.params.fp2_params.clone();

    let two = Fp2::from_u64(2, fq2_params.clone());
    let two_inv = two.inverse().expect("2 must be invertible in BN254 Fq2");

    let twist = bn254_twist();

    let mut coeffs = Vec::new();

    if let Some(msb) = highest_bit(&ALT_BN128_ATE_LOOP_COUNT) {
        // Skip the leading one bit.
        for bit_index in (0..msb).rev() {
            coeffs.push(doubling_step_for_flipped_miller_loop(
                &mut current,
                &two_inv,
                &twist,
            ));

            if get_bit(&ALT_BN128_ATE_LOOP_COUNT, bit_index) {
                coeffs.push(mixed_addition_step_for_flipped_miller_loop(
                    &mut current,
                    &q_affine,
                    &twist,
                ));
            }
        }
    }

    if ALT_BN128_ATE_IS_LOOP_COUNT_NEG {
        current.y = current.y.negate();
    }

    // Frobenius images Q1 and Q2.
    let q1 = q_affine.mul_by_q();
    let mut q2 = q1.mul_by_q();

    q2.y = q2.y.negate();

    coeffs.push(mixed_addition_step_for_flipped_miller_loop(
        &mut current,
        &q1,
        &twist,
    ));

    coeffs.push(mixed_addition_step_for_flipped_miller_loop(
        &mut current,
        &q2,
        &twist,
    ));

    AltBn128G2Precomp { qx, qy, coeffs }
}

/// Consume one sparse line from the Miller loop.
fn mul_by_line(f: &Fp12, coeff: &AltBn128AteEllCoeffs, p: &AltBn128G1Precomp) -> Fp12 {
    let ell_vw = coeff.ell_vw.mul_by_fp(&p.py);
    let ell_vv = coeff.ell_vv.mul_by_fp(&p.px);

    f.mul_by_024(&coeff.ell_0, &ell_vw, &ell_vv)
}

/// Ate Miller loop.
pub fn alt_bn128_ate_miller_loop(prec_p: &AltBn128G1Precomp, prec_q: &AltBn128G2Precomp) -> Fp12 {
    let params = bn254_fp12_params();

    let mut f = Fp12::one(params);

    let mut coeff_index = 0usize;

    if let Some(msb) = highest_bit(&ALT_BN128_ATE_LOOP_COUNT) {
        for bit_index in (0..msb).rev() {
            f = f.squared();

            f = mul_by_line(&f, &prec_q.coeffs[coeff_index], prec_p);

            coeff_index += 1;

            if get_bit(&ALT_BN128_ATE_LOOP_COUNT, bit_index) {
                f = mul_by_line(&f, &prec_q.coeffs[coeff_index], prec_p);

                coeff_index += 1;
            }
        }
    }

    if ALT_BN128_ATE_IS_LOOP_COUNT_NEG {
        f = f.inverse().expect("Miller result must be invertible");
    }

    // Final two Frobenius lines.
    for _ in 0..2 {
        f = mul_by_line(&f, &prec_q.coeffs[coeff_index], prec_p);

        coeff_index += 1;
    }

    f
}

/// First chunk of the final exponentiation:
///
/// (q^6 - 1)(q^2 + 1)
pub fn alt_bn128_final_exponentiation_first_chunk(elt: &Fp12) -> Fp12 {
    let params = elt.params.clone();

    // Conjugate: (c0, c1) -> (c0, -c1).
    let conjugate = Fp12::new(elt.c0.clone(), elt.c1.negate(), params);

    let inverse = elt
        .inverse()
        .expect("final exponentiation input must be nonzero");

    let c = conjugate.mul(&inverse);
    let d = c.frobenius_map(2);

    d.mul(&c)
}

/// Exponentiate by -z in the cyclotomic subgroup.
fn exp_by_neg_z(elt: &Fp12) -> Fp12 {
    let mut result = elt.cyclotomic_exp(&ALT_BN128_FINAL_EXPONENT_Z);

    if !ALT_BN128_FINAL_EXPONENT_IS_Z_NEG {
        result = result.unitary_inverse();
    }

    result
}

/// Last chunk of the final exponentiation.
///
/// This follows the libff BN254 addition chain.
pub fn alt_bn128_final_exponentiation_last_chunk(elt: &Fp12) -> Fp12 {
    let a = exp_by_neg_z(elt);

    let b = a.cyclotomic_squared();
    let c = b.cyclotomic_squared();
    let d = c.mul(&b);

    let e = exp_by_neg_z(&d);
    let f = e.cyclotomic_squared();
    let g = exp_by_neg_z(&f);

    let h = d.unitary_inverse();
    let i = g.unitary_inverse();

    let j = i.mul(&e);
    let k = j.mul(&h);
    let l = k.mul(&b);
    let m = k.mul(&e);

    let n = m.mul(elt);
    let o = l.frobenius_map(1);
    let p = o.mul(&n);

    let q = k.frobenius_map(2);
    let r = q.mul(&p);

    let s = elt.unitary_inverse();
    let t = s.mul(&l);

    let u = t.frobenius_map(3);
    u.mul(&r)
}

/// Full final exponentiation.
pub fn alt_bn128_final_exponentiation(elt: &Fp12) -> Fp12 {
    let first = alt_bn128_final_exponentiation_first_chunk(elt);

    alt_bn128_final_exponentiation_last_chunk(&first)
}

/// Ate pairing without final exponentiation.
pub fn alt_bn128_ate_pairing(p: &G1, q: &G2) -> Fp12 {
    let prec_p = alt_bn128_ate_precompute_g1(p);
    let prec_q = alt_bn128_ate_precompute_g2(q);

    alt_bn128_ate_miller_loop(&prec_p, &prec_q)
}

/// Reduced Tate/Ate pairing.
pub fn alt_bn128_reduced_pairing(p: &G1, q: &G2) -> Fp12 {
    let miller = alt_bn128_ate_pairing(p, q);

    alt_bn128_final_exponentiation(&miller)
}
