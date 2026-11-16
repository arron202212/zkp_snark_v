use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};
use crate::algebra::fields::prime_extension::fp4::{Fp4, Fp4Params};

use super::mnt4_g1::{Mnt4G1, Mnt4G1Params};
use super::mnt4_g2::{Mnt4G2, mnt4_fq_params, mnt4_fq2_params};

/// MNT4-298 base field.
pub type Mnt4Fq = Fp<5>;

/// MNT4-298 quadratic extension.
pub type Mnt4Fq2 = Fp2<5>;

/// MNT4-298 quartic extension.
pub type Mnt4Fq4 = Fp4<5>;

/// Target group.
pub type Mnt4GT = Mnt4Fq4;

/// Fq4 parameters.
pub type Mnt4Fq4Params = Fp4Params<5>;

// ============================================================================
// Ate loop parameters
// ============================================================================

/// Ate loop count.
///
/// Little-endian limbs.
const MNT4_ATE_LOOP_COUNT: [u64; 3] = [
    0x0dc9_a1b6_7166_0000,
    0x4660_9756_bec2_a33f,
    0x0000_0000_001e_ef55,
];

const MNT4_ATE_LOOP_BITS: usize = 192;

// ============================================================================
// Final exponent
// ============================================================================

/// Full final exponent:
/// ```text
///     (q^4 - 1) / r
/// ```
/// Little-endian u64 limbs.
const MNT4_FINAL_EXPONENT: [u64; 14] = [
    0xe7e6_9541_c598_0000,
    0x065e_41b0_1242_2619,
    0x36c7_a580_1f55_ac8e,
    0xa570_04bb_9a1e_7a60,
    0xb351_384d_8485_e987,
    0x11ce_e431_be83_48de,
    0x4e61_f58b_68ea_590b,
    0x8507_4066_8186_8763,
    0xba73_85b8_d20f_ea4e,
    0x95e8_7dcb_6c97_234c,
    0x42e0_0292_64c0_324a,
    0x35ac_ca5a_07b2_394f,
    0xefe2_16b3_7afb_6d30,
    0x343c_7ac3_174c_87a1,
];

// ============================================================================
// Fq4 initialization
// ============================================================================

/// Official libff MNT4 Fq4 Frobenius coefficients.
///
/// Each entry is a canonical Fq integer in little-endian limbs.
const MNT4_FQ4_FROBENIUS_C1: [[u64; 5]; 4] = [
    [1, 0, 0, 0, 0],
    [
        0x94dd5d7def6980c4,
        0x8cd9fae5c1f7bdcf,
        0x8d534beb17daf751,
        0x9916dfdcc2fd1f96,
        0x0000000f73779fe0,
    ],
    [
        0xc90cd65a71660000,
        0x41a9e35e51200e12,
        0xcaeec9635d1330ea,
        0xa266249da7b0548e,
        0x000003bcf7bcd473,
    ],
    [
        0x342f78dc81fc7f3d,
        0xb4cfe8788f285043,
        0x3d9b7d7845383998,
        0x094f44c0e4b334f8,
        0x000003ad84453493,
    ],
];

/// Convert the canonical Frobenius coefficients into Fp elements.
///
/// Fp::new() performs the conversion to Montgomery representation.
fn mnt4_fq4_frobenius_coeffs(fq_params: &FieldParams<5>) -> [Fp<5>; 4] {
    std::array::from_fn(|i| Fp::<5>::new(MNT4_FQ4_FROBENIUS_C1[i], fq_params.clone()))
}

/// Construct MNT4 Fq4 parameters.
///
/// Fq2 = Fq[u] / (u² - 17)
/// Fq4 = Fq2[v] / (v² - u)
pub fn mnt4_fq4_params() -> Mnt4Fq4Params {
    let fq_params = mnt4_fq_params();
    let fq2_params = mnt4_fq2_params();

    let non_residue = Mnt4Fq::from_u64(17, fq_params.clone());

    let frobenius_coeffs_c1 = mnt4_fq4_frobenius_coeffs(&fq_params);

    Mnt4Fq4Params::new(fq2_params, non_residue, frobenius_coeffs_c1)
}

// ============================================================================
// Precomputation structures
// ============================================================================

#[derive(Clone, Debug)]
pub struct Mnt4AteG1Precomp {
    pub px: Mnt4Fq,
    pub py: Mnt4Fq,

    pub px_twist: Mnt4Fq2,
    pub py_twist: Mnt4Fq2,

    pub fq2_params: Fp2Params<5>,

    pub is_zero: bool,
}

#[derive(Clone, Debug)]
pub struct Mnt4AteDblCoeffs {
    pub c_h: Mnt4Fq2,
    pub c_4c: Mnt4Fq2,
    pub c_j: Mnt4Fq2,
    pub c_l: Mnt4Fq2,
}

#[derive(Clone, Debug)]
pub struct Mnt4AteAddCoeffs {
    pub c_l1: Mnt4Fq2,
    pub c_rz: Mnt4Fq2,
}

#[derive(Clone, Debug)]
pub struct Mnt4AteG2Precomp {
    pub qx: Mnt4Fq2,
    pub qy: Mnt4Fq2,
    pub qy2: Mnt4Fq2,

    pub qx_over_twist: Mnt4Fq2,
    pub qy_over_twist: Mnt4Fq2,

    pub dbl_coeffs: Vec<Mnt4AteDblCoeffs>,
    pub add_coeffs: Vec<Mnt4AteAddCoeffs>,

    pub fq4_params: Mnt4Fq4Params,

    pub is_zero: bool,
}

// ============================================================================
// Extended projective G2 representation
// ============================================================================

/// Extended projective representation used by the flipped Miller loop.
///
/// Invariant:
///  ```text
///     T = Z²
///  ```
#[derive(Clone, Debug)]
struct ExtendedMnt4G2 {
    x: Mnt4Fq2,
    y: Mnt4Fq2,
    z: Mnt4Fq2,
    t: Mnt4Fq2,
}

impl ExtendedMnt4G2 {
    fn new(q: &Mnt4G2) -> Self {
        let affine = q.to_affine_element();

        let one = Mnt4Fq2::one(affine.params.fq2_params.clone());

        Self {
            x: affine.x,
            y: affine.y,
            z: one.clone(),
            t: one,
        }
    }

    fn check_invariant(&self) -> bool {
        self.t == self.z.square()
    }
}

// ============================================================================
// Small helpers
// ============================================================================

fn fq2_two(params: &Fp2Params<5>) -> Mnt4Fq2 {
    Mnt4Fq2::from_u64(2, params.clone())
}

fn fq2_four(params: &Fp2Params<5>) -> Mnt4Fq2 {
    let two = fq2_two(params);
    two.add(&two)
}

fn fq2_eight(params: &Fp2Params<5>) -> Mnt4Fq2 {
    let four = fq2_four(params);
    four.add(&four)
}

fn loop_bit(index: usize) -> bool {
    let limb = index / 64;
    let bit = index % 64;

    (MNT4_ATE_LOOP_COUNT[limb] & (1u64 << bit)) != 0
}

// ============================================================================
// Flipped Miller loop: doubling step
// ============================================================================

fn doubling_step_for_flipped_miller_loop(current: &mut ExtendedMnt4G2) -> Mnt4AteDblCoeffs {
    let x = current.x.clone();
    let y = current.y.clone();
    let z = current.z.clone();
    let t = current.t.clone();

    let params = &x.params;

    let a = t.square();
    let b = x.square();
    let c = y.square();
    let d = c.square();

    let e = x.add(&c).square().sub(&b).sub(&d);

    let three_b = b.add(&b).add(&b);

    // MNT4 twist curve coefficient a = 34.
    let curve_a = Mnt4Fq2::new(
        Mnt4Fq::from_u64(34, params.fp_params.clone()),
        Mnt4Fq::zero(params.fp_params.clone()),
        params.clone(),
    );

    let f = three_b.add(&curve_a.mul(&a));
    let g = f.square();

    let four_e = e.add(&e).add(&e).add(&e);

    let x3 = four_e.negate().add(&g);

    let two_e_minus_x3 = e.add(&e).sub(&x3);

    let y3 = fq2_eight(params)
        .mul(&d)
        .negate()
        .add(&f.mul(&two_e_minus_x3));

    let z3 = y.add(&z).square().sub(&c).sub(&z.square());

    let t3 = z3.square();

    let c_h = z3.add(&t).square().sub(&t3).sub(&a);

    let c_4c = c.add(&c).add(&c).add(&c);

    let c_j = f.add(&t).square().sub(&g).sub(&a);

    let c_l = f.add(&x).square().sub(&g).sub(&b);

    current.x = x3;
    current.y = y3;
    current.z = z3;
    current.t = t3;

    debug_assert!(current.check_invariant());

    Mnt4AteDblCoeffs {
        c_h,
        c_4c,
        c_j,
        c_l,
    }
}

// ============================================================================
// Flipped Miller loop: mixed addition step
// ============================================================================

fn mixed_addition_step_for_flipped_miller_loop(
    base_x: &Mnt4Fq2,
    base_y: &Mnt4Fq2,
    base_y_squared: &Mnt4Fq2,
    current: &mut ExtendedMnt4G2,
) -> Mnt4AteAddCoeffs {
    let x = current.x.clone();
    let y = current.y.clone();
    let z = current.z.clone();
    let t = current.t.clone();

    let b = base_x.mul(&t);

    let d = base_y.add(&z).square().sub(base_y_squared).sub(&t).mul(&t);

    let h = b.sub(&x);
    let i = h.square();

    let e = i.add(&i).add(&i).add(&i);

    let j = h.mul(&e);
    let v = x.mul(&e);

    let l1 = d.sub(&y.add(&y));

    let x3 = l1.square().sub(&j).sub(&v.add(&v));

    let y3 = l1.mul(&v.sub(&x3)).sub(&y.add(&y).mul(&j));

    let z3 = z.add(&h).square().sub(&t).sub(&i);

    let t3 = z3.square();

    current.x = x3;
    current.y = y3;
    current.z = z3.clone();
    current.t = t3;

    debug_assert!(current.check_invariant());

    Mnt4AteAddCoeffs { c_l1: l1, c_rz: z3 }
}

// ============================================================================
// G1 precomputation
// ============================================================================

pub fn mnt4_ate_precompute_g1(p: &Mnt4G1) -> Mnt4AteG1Precomp {
    let fq2_params = mnt4_fq2_params();

    let affine = p.to_affine_point().expect("G1 affine conversion failed");

    let zero = Mnt4Fq::zero(p.params.fq_params.clone());

    let twist = Mnt4Fq2::new(
        zero,
        Mnt4Fq::one(p.params.fq_params.clone()),
        fq2_params.clone(),
    );

    let px_twist = twist.mul_by_fp(&affine.x);
    let py_twist = twist.mul_by_fp(&affine.y);

    Mnt4AteG1Precomp {
        px: affine.x,
        py: affine.y,
        px_twist,
        py_twist,
        fq2_params,
        is_zero: p.is_zero(),
    }
}

// ============================================================================
// G2 precomputation
// ============================================================================

pub fn mnt4_ate_precompute_g2(q: &Mnt4G2) -> Mnt4AteG2Precomp {
    let affine = q.to_affine_element();

    let params = &affine.params.fq2_params;

    let twist = Mnt4Fq2::new(
        Mnt4Fq::zero(affine.params.fq_params.clone()),
        Mnt4Fq::one(affine.params.fq_params.clone()),
        params.clone(),
    );

    let twist_inv = twist.inverse().expect("MNT4 twist must be invertible");

    let qx = affine.x.clone();
    let qy = affine.y.clone();
    let qy2 = qy.square();

    let qx_over_twist = qx.mul(&twist_inv);
    let qy_over_twist = qy.mul(&twist_inv);

    let fq4_params = mnt4_fq4_params();

    let mut result = Mnt4AteG2Precomp {
        qx: qx.clone(),
        qy: qy.clone(),
        qy2: qy2.clone(),

        qx_over_twist,
        qy_over_twist,

        dbl_coeffs: Vec::new(),
        add_coeffs: Vec::new(),

        fq4_params,

        is_zero: q.is_zero(),
    };

    if result.is_zero {
        return result;
    }

    let mut r = ExtendedMnt4G2::new(&affine);

    let mut found_one = false;

    for i in (0..MNT4_ATE_LOOP_BITS).rev() {
        let bit = loop_bit(i);

        if !found_one {
            found_one |= bit;
            continue;
        }

        let dc = doubling_step_for_flipped_miller_loop(&mut r);

        result.dbl_coeffs.push(dc);

        if bit {
            let ac = mixed_addition_step_for_flipped_miller_loop(&qx, &qy, &qy2, &mut r);

            result.add_coeffs.push(ac);
        }
    }

    result
}

// ============================================================================
// Miller loop
// ============================================================================

pub fn mnt4_ate_miller_loop(prec_p: &Mnt4AteG1Precomp, prec_q: &Mnt4AteG2Precomp) -> Mnt4Fq4 {
    let params = prec_q.fq4_params.clone();

    if prec_p.is_zero || prec_q.is_zero {
        return Mnt4Fq4::one(params);
    }

    let zero = Mnt4Fq::zero(prec_p.px.params.clone());

    let l1_coeff =
        Mnt4Fq2::new(prec_p.px.clone(), zero, prec_p.fq2_params.clone()).sub(&prec_q.qx_over_twist);

    let mut f = Mnt4Fq4::one(params.clone());

    let mut dbl_idx = 0usize;
    let mut add_idx = 0usize;

    let mut found_one = false;

    for i in (0..MNT4_ATE_LOOP_BITS).rev() {
        let bit = loop_bit(i);

        if !found_one {
            found_one |= bit;
            continue;
        }

        let dc = &prec_q.dbl_coeffs[dbl_idx];
        dbl_idx += 1;

        let g_rr_c0 = dc
            .c_4c
            .negate()
            .sub(&dc.c_j.mul(&prec_p.px_twist))
            .add(&dc.c_l);

        let g_rr_c1 = dc.c_h.mul(&prec_p.py_twist);

        let g_rr = Mnt4Fq4::new(g_rr_c0, g_rr_c1, params.clone());

        f = f.square().mul(&g_rr);

        if bit {
            let ac = &prec_q.add_coeffs[add_idx];
            add_idx += 1;

            let g_rq_c0 = ac.c_rz.mul(&prec_p.py_twist);

            let g_rq_c1 = prec_q
                .qy_over_twist
                .mul(&ac.c_rz)
                .add(&l1_coeff.mul(&ac.c_l1))
                .negate();

            let g_rq = Mnt4Fq4::new(g_rq_c0, g_rq_c1, params.clone());

            f = f.mul(&g_rq);
        }
    }

    debug_assert_eq!(dbl_idx, prec_q.dbl_coeffs.len());

    debug_assert_eq!(add_idx, prec_q.add_coeffs.len());

    f
}

// ============================================================================
// Double Miller loop
// ============================================================================

/// Compute:
/// ```text
///     Miller(P1, Q1) * Miller(P2, Q2)
/// ```
/// This version reuses the single Miller loop.
/// It can later be replaced by the fused implementation from libff.
pub fn mnt4_ate_double_miller_loop(
    prec_p1: &Mnt4AteG1Precomp,
    prec_q1: &Mnt4AteG2Precomp,
    prec_p2: &Mnt4AteG1Precomp,
    prec_q2: &Mnt4AteG2Precomp,
) -> Mnt4Fq4 {
    let f1 = mnt4_ate_miller_loop(prec_p1, prec_q1);
    let f2 = mnt4_ate_miller_loop(prec_p2, prec_q2);

    f1.mul(&f2)
}

// ============================================================================
// Final exponentiation
// ============================================================================

/// Full final exponentiation.
///
/// Computes:
/// ```text
///     elt ^ ((q^4 - 1) / r)
/// ```
/// Correctness-first implementation.
/// It does not yet use libff's optimized Frobenius/cyclotomic decomposition.
pub fn mnt4_final_exponentiation(elt: &Mnt4Fq4) -> Mnt4GT {
    elt.pow(&MNT4_FINAL_EXPONENT)
}

// ============================================================================
// Public pairing API
// ============================================================================

/// Miller loop for a pair of points.
pub fn mnt4_miller_loop(p: &Mnt4G1, q: &Mnt4G2) -> Mnt4Fq4 {
    let prec_p = mnt4_ate_precompute_g1(p);
    let prec_q = mnt4_ate_precompute_g2(q);

    mnt4_ate_miller_loop(&prec_p, &prec_q)
}

/// Reduced pairing.
///```text
///     e(P, Q) = Miller(P, Q)^((q^4 - 1) / r)
///```
pub fn mnt4_reduced_pairing(p: &Mnt4G1, q: &Mnt4G2) -> Mnt4GT {
    let f = mnt4_miller_loop(p, q);

    mnt4_final_exponentiation(&f)
}

/// Alias matching the usual pairing API.
pub fn mnt4_pairing(p: &Mnt4G1, q: &Mnt4G2) -> Mnt4GT {
    mnt4_reduced_pairing(p, q)
}

/// Compute the product of two reduced pairings.
pub fn mnt4_double_reduced_pairing(p1: &Mnt4G1, q1: &Mnt4G2, p2: &Mnt4G1, q2: &Mnt4G2) -> Mnt4GT {
    let prec_p1 = mnt4_ate_precompute_g1(p1);
    let prec_q1 = mnt4_ate_precompute_g2(q1);

    let prec_p2 = mnt4_ate_precompute_g1(p2);
    let prec_q2 = mnt4_ate_precompute_g2(q2);

    let f = mnt4_ate_double_miller_loop(&prec_p1, &prec_q1, &prec_p2, &prec_q2);

    mnt4_final_exponentiation(&f)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn g1_generator() -> Mnt4G1 {
        let params = std::sync::Arc::new(Mnt4G1Params::new(mnt4_fq_params()));

        Mnt4G1::one(params)
    }

    fn g2_generator() -> Mnt4G2 {
        Mnt4G2::generator()
    }

    #[test]
    fn test_mnt4_fq4_params() {
        let params = mnt4_fq4_params();

        let one = Mnt4Fq4::one(params.clone());
        let zero = Mnt4Fq4::zero(params.clone());

        assert!(one.is_one());
        assert!(zero.is_zero());

        assert_eq!(params.non_residue.to_canonical()[0], 17);

        assert_eq!(params.frobenius_coeffs_c1[0].to_canonical()[0], 1);

        let q_minus_one = Mnt4Fq::new(
            [
                0xc90cd65a71660000,
                0x41a9e35e51200e12,
                0xcaeec9635d1330ea,
                0xa266249da7b0548e,
                0x000003bcf7bcd473,
            ],
            mnt4_fq_params(),
        );

        assert_eq!(params.frobenius_coeffs_c1[2], q_minus_one);
    }

    #[test]
    fn test_mnt4_g1_precomputation() {
        let p = g1_generator();

        let prec = mnt4_ate_precompute_g1(&p);

        assert!(!prec.is_zero);
    }

    #[test]
    fn test_mnt4_g2_precomputation() {
        let q = g2_generator();

        let prec = mnt4_ate_precompute_g2(&q);

        assert!(!prec.is_zero);
        assert!(!prec.dbl_coeffs.is_empty());
    }

    #[test]
    fn test_mnt4_miller_loop_identity() {
        let p = g1_generator();
        let q = g2_generator();

        let p_zero = Mnt4G1::zero(std::sync::Arc::clone(&p.params));

        let q_zero = Mnt4G2::zero(q.params.clone());

        let result1 = mnt4_miller_loop(&p_zero, &q);
        let result2 = mnt4_miller_loop(&p, &q_zero);

        assert!(result1.is_one());
        assert!(result2.is_one());
    }

    #[test]
    fn test_mnt4_pairing_identity() {
        let p = g1_generator();
        let q = g2_generator();

        let p_zero = Mnt4G1::zero(std::sync::Arc::clone(&p.params));

        let q_zero = Mnt4G2::zero(q.params.clone());

        assert!(mnt4_pairing(&p_zero, &q).is_one());
        assert!(mnt4_pairing(&p, &q_zero).is_one());
    }

    #[test]
    fn test_final_exponentiation_one() {
        let params = mnt4_fq4_params();
        let one = Mnt4Fq4::one(params);

        assert!(mnt4_final_exponentiation(&one).is_one());
    }

    #[test]
    fn test_mnt4_pairing_nonzero() {
        let p = g1_generator();
        let q = g2_generator();

        let e = mnt4_pairing(&p, &q);

        assert!(
            !e.is_one(),
            "pairing of non-zero generators must not be one"
        );
    }

    fn assert_fp4_eq(a: &Mnt4Fq4, b: &Mnt4Fq4) {
        assert_eq!(a.c0, b.c0, "Fp4 c0 mismatch");
        assert_eq!(a.c1, b.c1, "Fp4 c1 mismatch");
    }

    #[test]
    fn test_mnt4_pairing_bilinearity_g1() {
        let p = g1_generator();
        let q = g2_generator();

        let two_p = p.mul_scalar([2, 0, 0, 0, 0]);

        let e1 = mnt4_pairing(&p, &q);
        let e2 = mnt4_pairing(&two_p, &q);

        let expected = e1.square();

        assert_fp4_eq(&e2, &expected);
    }
    #[test]
    fn test_mnt4_pairing_bilinearity_g2() {
        let p = g1_generator();
        let q = g2_generator();

        let two_q = q.scalar_mul(&[2]);

        let e1 = mnt4_pairing(&p, &q);
        let e2 = mnt4_pairing(&p, &two_q);

        let expected = e1.square();

        assert_fp4_eq(&e2, &expected);
    }

    #[test]
    fn test_mnt4_fq4_extension_relation() {
        let params = mnt4_fq4_params();
        let fq2_params = params.fp2_params.clone();
        let fq_params = fq2_params.fp_params.clone();

        let zero = Mnt4Fq2::zero(fq2_params.clone());
        let one = Mnt4Fq2::one(fq2_params.clone());

        // u = 0 + 1*u
        let u = Mnt4Fq2::new(
            Mnt4Fq::zero(fq_params.clone()),
            Mnt4Fq::one(fq_params),
            fq2_params,
        );

        // v = 0 + 1*v
        let v = Mnt4Fq4::new(zero.clone(), one, params.clone());

        let expected = Mnt4Fq4::new(u, zero, params);

        assert!(
            v.square().equals(&expected),
            "MNT4 Fq4 must satisfy v^2 = u"
        );
    }

    #[test]
    fn test_mnt4_fq4_frobenius() {
        let params = mnt4_fq4_params();
        let fq2_params = params.fp2_params.clone();
        let fq_params = fq2_params.fp_params.clone();

        let a = Mnt4Fq2::new(
            Mnt4Fq::from_u64(2, fq_params.clone()),
            Mnt4Fq::from_u64(3, fq_params.clone()),
            fq2_params.clone(),
        );

        let b = Mnt4Fq2::new(
            Mnt4Fq::from_u64(5, fq_params.clone()),
            Mnt4Fq::from_u64(7, fq_params),
            fq2_params,
        );

        let x = Mnt4Fq4::new(a, b, params);

        let q = [
            0xc90c_d65a_7166_0001,
            0x41a9_e35e_5120_0e12,
            0xcaee_c963_5d13_30ea,
            0xa266_249d_a7b0_548e,
            0x0000_03bc_f7bc_d473,
        ];

        let expected = x.pow(&q);
        let actual = x.frobenius_map(1);

        assert!(
            actual.equals(&expected),
            "Fq4 Frobenius map must equal exponentiation by q"
        );

        assert!(
            x.frobenius_map(4).equals(&x),
            "Fq4 Frobenius^4 must be identity"
        );
    }
}
