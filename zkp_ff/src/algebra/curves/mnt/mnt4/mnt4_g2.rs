use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};

pub type Mnt4Fq = Fp<5>;
pub type Mnt4Fq2 = Fp2<5>;

/// MNT4 G2 parameters.
#[derive(Clone, Debug)]
pub struct Mnt4G2Params {
    pub fq_params: FieldParams<5>,
    pub fq2_params: Fp2Params<5>,

    /// Twist curve coefficient a.
    ///
    /// twist_a = (34, 0)
    pub twist_coeff_a: Mnt4Fq2,

    /// Twist curve coefficient b.
    ///
    /// twist_b = (0, 17 * curve_b)
    pub twist_coeff_b: Mnt4Fq2,

    /// Optimized multiplication by twist coefficient a.
    pub twist_mul_by_a_c0: Mnt4Fq,
    pub twist_mul_by_a_c1: Mnt4Fq,

    /// Optimized multiplication by twist coefficient b.
    pub twist_mul_by_b_c0: Mnt4Fq,
    pub twist_mul_by_b_c1: Mnt4Fq,

    /// Frobenius multiplication coefficients.
    pub twist_mul_by_q_x: Mnt4Fq,
    pub twist_mul_by_q_y: Mnt4Fq,

    /// G2 cofactor.
    pub cofactor: [u64; 5],
}

/// MNT4 G2 point.
///
/// Uses ordinary projective coordinates:
///  ```text
///     x = X / Z
///     y = Y / Z
///  ```
/// This is NOT Jacobian representation.
#[derive(Clone, Debug)]
pub struct Mnt4G2 {
    pub x: Mnt4Fq2,
    pub y: Mnt4Fq2,
    pub z: Mnt4Fq2,

    pub params: Mnt4G2Params,
}

impl PartialEq for Mnt4G2 {
    fn eq(&self, other: &Self) -> bool {
        self.equals(other)
    }
}

impl Eq for Mnt4G2 {}

// ============================================================================
// Field constants
// ============================================================================

/// Construct a canonical Fq element from five little-endian limbs.
fn fq(value: [u64; 5], params: &FieldParams<5>) -> Mnt4Fq {
    Mnt4Fq::new(value, params.clone())
}

/// Construct an Fq2 element from two canonical Fq components.
fn fq2(c0: [u64; 5], c1: [u64; 5], params: &Fp2Params<5>) -> Mnt4Fq2 {
    Mnt4Fq2::new(
        fq(c0, &params.fp_params),
        fq(c1, &params.fp_params),
        params.clone(),
    )
}

fn fq_zero(params: &FieldParams<5>) -> Mnt4Fq {
    Mnt4Fq::zero(params.clone())
}

fn fq_one(params: &FieldParams<5>) -> Mnt4Fq {
    Mnt4Fq::one(params.clone())
}

fn fq2_zero(params: &Fp2Params<5>) -> Mnt4Fq2 {
    Mnt4Fq2::zero(params.clone())
}

fn fq2_one(params: &Fp2Params<5>) -> Mnt4Fq2 {
    Mnt4Fq2::one(params.clone())
}

// ============================================================================
// MNT4 G2 initialization
// ============================================================================

/// MNT4 base-field modulus.
///
/// Decimal:
///
/// 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
pub fn mnt4_fq_params() -> FieldParams<5> {
    FieldParams::new([
        0xc90c_d65a_7166_0001,
        0x41a9_e35e_5120_0e12,
        0xcaee_c963_5d13_30ea,
        0xa266_249d_a7b0_548e,
        0x0000_03bc_f7bc_d473,
    ])
}

/// MNT4 Fq2 parameters.
///
/// Fq2 = Fq[u] / (u^2 - 17)
///
/// Frobenius:
///
/// (a + bu)^q = a - bu
pub fn mnt4_fq2_params() -> Fp2Params<5> {
    let fq_params = mnt4_fq_params();

    let non_residue = Mnt4Fq::from_u64(17, fq_params.clone());

    let one = fq_one(&fq_params);
    let minus_one = one.negate();

    Fp2Params::new(fq_params, non_residue, [one, minus_one])
}

/// Initialize MNT4 G2 parameters.
pub fn mnt4_g2_params() -> Mnt4G2Params {
    let fq_params = mnt4_fq_params();
    let fq2_params = mnt4_fq2_params();

    let zero = fq_zero(&fq_params);

    // ------------------------------------------------------------------------
    // Twist coefficient a
    //
    // twist_a = (34, 0)
    // ------------------------------------------------------------------------

    let twist_coeff_a = Mnt4Fq2::new(
        Mnt4Fq::from_u64(34, fq_params.clone()),
        zero.clone(),
        fq2_params.clone(),
    );

    // ------------------------------------------------------------------------
    // Twist coefficient b
    //
    // twist_b = (0, 17 * curve_b)
    //
    // 17 * curve_b mod q
    // ------------------------------------------------------------------------

    let twist_coeff_b = fq2(
        [
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ],
        [
            0x6a3b_6efb_12b3_0d16,
            0x5c24_b1b1_a912_ec4d,
            0x1ecd_2ff0_038a_b06a,
            0xe6e2_23b5_5400_63ee,
            0x0000_0087_788d_2a0e,
        ],
        &fq2_params,
    );

    // ------------------------------------------------------------------------
    // Optimized multiplication constants.
    //
    // mul_by_a(elt):
    //   (34 * elt.c0, 34 * elt.c1)
    //
    // mul_by_b(elt):
    //   (b_c0 * elt.c1, b_c1 * elt.c0)
    // ------------------------------------------------------------------------

    let twist_mul_by_a_c0 = Mnt4Fq::from_u64(34, fq_params.clone());

    let twist_mul_by_a_c1 = Mnt4Fq::from_u64(34, fq_params.clone());

    // curve_b * 17^2 mod q
    let twist_mul_by_b_c0 = fq(
        [
            0x7bd8_b1f7_5b17_de74,
            0x9b1c_060f_9801_94fe,
            0x75c2_9c29_820f_553b,
            0x1038_15cf_44a5_f9b2,
            0x0000_0185_11e6_2216,
        ],
        &fq_params,
    );

    // curve_b * 17 mod q
    let twist_mul_by_b_c1 = fq(
        [
            0x6a3b_6efb_12b3_0d16,
            0x5c24_b1b1_a912_ec4d,
            0x1ecd_2ff0_038a_b06a,
            0xe6e2_23b5_5400_63ee,
            0x0000_0087_788d_2a0e,
        ],
        &fq_params,
    );

    // ------------------------------------------------------------------------
    // Frobenius multiplication constants.
    // ------------------------------------------------------------------------

    // q - 1
    let twist_mul_by_q_x = fq(
        [
            0xc90c_d65a_7166_0000,
            0x41a9_e35e_5120_0e12,
            0xcaee_c963_5d13_30ea,
            0xa266_249d_a7b0_548e,
            0x0000_03bc_f7bc_d473,
        ],
        &fq_params,
    );

    let twist_mul_by_q_y = fq(
        [
            0x94dd_5d7d_ef69_80c4,
            0x8cd9_fae5_c1f7_bdcf,
            0x8d53_4beb_17da_f751,
            0x9916_dfdc_c2fd_1f96,
            0x0000_000f_7377_9fe0,
        ],
        &fq_params,
    );

    // ------------------------------------------------------------------------
    // G2 cofactor.
    // ------------------------------------------------------------------------

    let cofactor = [
        0xd6d6_7810_e2cc_0001,
        0x880a_7ab5_0fe2_b151,
        0xcaee_c963_5d32_203f,
        0xa266_249d_a7b0_548e,
        0x0000_03bc_f7bc_d473,
    ];

    Mnt4G2Params {
        fq_params,
        fq2_params,
        twist_coeff_a,
        twist_coeff_b,
        twist_mul_by_a_c0,
        twist_mul_by_a_c1,
        twist_mul_by_b_c0,
        twist_mul_by_b_c1,
        twist_mul_by_q_x,
        twist_mul_by_q_y,
        cofactor,
    }
}

// ============================================================================
// G2 implementation
// ============================================================================

impl Mnt4G2 {
    pub fn new(x: Mnt4Fq2, y: Mnt4Fq2, z: Mnt4Fq2, params: Mnt4G2Params) -> Self {
        Self { x, y, z, params }
    }

    /// Point at infinity.
    ///
    /// Canonical representation:
    ///
    /// (0, 1, 0)
    pub fn zero(params: Mnt4G2Params) -> Self {
        let x = fq2_zero(&params.fq2_params);
        let y = fq2_one(&params.fq2_params);
        let z = fq2_zero(&params.fq2_params);

        Self::new(x, y, z, params)
    }

    /// Generator.
    pub fn one(params: Mnt4G2Params) -> Self {
        let x = fq2(
            [
                0xba79_03df_6c09_a9a4,
                0x2cb1_4f01_a931_e72d,
                0x9001_f205_151e_12a7,
                0x6605_71ff_542f_2ef8,
                0x0000_0371_7804_91c5,
            ],
            [
                0x968b_e32c_0ae0_a989,
                0x0330_2bb6_c02c_712c,
                0x697c_851f_002f_5763,
                0xda16_5def_8380_81af,
                0x0000_004b_a59a_3f72,
            ],
            &params.fq2_params,
        );

        let y = fq2(
            [
                0x0d95_22bc_a4e7_9f22,
                0xec98_f0f6_10a5_aafd,
                0xd31e_5c4b_3b2e_0b60,
                0xaad8_68a1_c47d_6605,
                0x0000_004b_471f_33ff,
            ],
            [
                0x84fb_61a3_cbf0_e0f1,
                0x3b51_68ed_8d75_c7c4,
                0xcb7d_982f_78ec_9cfc,
                0xa503_1f3f_81a5_c100,
                0x0000_0355_d05a_1c69,
            ],
            &params.fq2_params,
        );

        let z = fq2_one(&params.fq2_params);

        Self::new(x, y, z, params)
    }

    pub fn generator() -> Self {
        Self::one(mnt4_g2_params())
    }

    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.z.is_zero()
    }

    pub fn is_special(&self) -> bool {
        self.is_zero() || self.z.is_one()
    }

    /// Convert to affine coordinates.
    ///```text
    /// x = X/Z
    /// y = Y/Z
    ///```
    pub fn to_affine(&mut self) {
        if self.is_zero() {
            self.x = fq2_zero(&self.params.fq2_params);
            self.y = fq2_one(&self.params.fq2_params);
            self.z = fq2_zero(&self.params.fq2_params);
            return;
        }

        if self.z.is_one() {
            return;
        }

        let z_inv = self
            .z
            .inverse()
            .expect("nonzero projective Z must be invertible");

        self.x = self.x.mul(&z_inv);
        self.y = self.y.mul(&z_inv);
        self.z = fq2_one(&self.params.fq2_params);
    }

    pub fn to_affine_element(&self) -> Self {
        let mut result = self.clone();
        result.to_affine();
        result
    }

    /// Compare projective points.
    ///
    /// X1/Z1 == X2/Z2
    /// Y1/Z1 == Y2/Z2
    pub fn equals(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero();
        }

        if other.is_zero() {
            return false;
        }

        let x1z2 = self.x.mul(&other.z);
        let x2z1 = other.x.mul(&self.z);

        if x1z2 != x2z1 {
            return false;
        }

        let y1z2 = self.y.mul(&other.z);
        let y2z1 = other.y.mul(&self.z);

        y1z2 == y2z1
    }

    // ------------------------------------------------------------------------
    // Optimized multiplication by twist coefficients
    // ------------------------------------------------------------------------

    /// Multiply an Fq2 element by twist coefficient a.
    ///
    /// (a0, a1) -> (34*a0, 34*a1)
    pub fn mul_by_a(&self, elt: &Mnt4Fq2) -> Mnt4Fq2 {
        Mnt4Fq2::new(
            elt.c0.mul(&self.params.twist_mul_by_a_c0),
            elt.c1.mul(&self.params.twist_mul_by_a_c1),
            self.params.fq2_params.clone(),
        )
    }

    /// Multiply an Fq2 element by twist coefficient b.
    ///
    /// (a0, a1) -> (b_c0*a1, b_c1*a0)
    pub fn mul_by_b(&self, elt: &Mnt4Fq2) -> Mnt4Fq2 {
        Mnt4Fq2::new(
            self.params.twist_mul_by_b_c0.mul(&elt.c1),
            self.params.twist_mul_by_b_c1.mul(&elt.c0),
            self.params.fq2_params.clone(),
        )
    }

    // ------------------------------------------------------------------------
    // Negation
    // ------------------------------------------------------------------------

    pub fn negate(&self) -> Self {
        if self.is_zero() {
            return Self::zero(self.params.clone());
        }

        Self::new(
            self.x.clone(),
            self.y.negate(),
            self.z.clone(),
            self.params.clone(),
        )
    }

    // ------------------------------------------------------------------------
    // Doubling
    // ------------------------------------------------------------------------

    /// Point doubling.
    ///
    /// This follows libff's ordinary-projective doubling formula.
    pub fn double(&self) -> Self {
        if self.is_zero() || self.y.is_zero() {
            return Self::zero(self.params.clone());
        }

        let xx = self.x.square();
        let zz = self.z.square();

        // w = a*Z^2 + 3*X^2
        let three = Mnt4Fq2::from_u64(3, self.params.fq2_params.clone());

        let w = self
            .mul_by_a(&zz)
            .add(&xx.mul_by_fp(&Mnt4Fq::from_u64(3, self.params.fq_params.clone())));

        // s = 2*Y*Z
        let two = Mnt4Fq2::from_u64(2, self.params.fq2_params.clone());

        let s = self.y.mul(&self.z).mul(&two);
        let ss = s.square();
        let sss = s.mul(&ss);

        // R = Y*s
        let r = self.y.mul(&s);
        let rr = r.square();

        // B = (X+R)^2 - XX - RR
        let b = self.x.add(&r).square().sub(&xx).sub(&rr);

        // h = w^2 - 2*B
        let h = w.square().sub(&b.mul(&two));

        // X3 = h*s
        let x3 = h.mul(&s);

        // Y3 = w*(B-h) - 2*RR
        let y3 = w.mul(&b.sub(&h)).sub(&rr.mul(&two));

        // Z3 = s^3
        let z3 = sss;

        let _ = three;

        Self::new(x3, y3, z3, self.params.clone())
    }

    pub fn dbl(&self) -> Self {
        self.double()
    }

    // ------------------------------------------------------------------------
    // Addition
    // ------------------------------------------------------------------------

    /// Complete ordinary-projective addition.
    ///
    /// Formula follows mnt4_g2.cpp.
    pub fn add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        if self.equals(other) {
            return self.double();
        }

        // Same affine x and opposite affine y.
        //
        // Avoid the degenerate all-zero output of the generic formula.
        if self.x.mul(&other.z) == other.x.mul(&self.z)
            && self.y.mul(&other.z) == other.y.negate().mul(&self.z)
        {
            return Self::zero(self.params.clone());
        }

        let y1z2 = self.y.mul(&other.z);
        let x1z2 = self.x.mul(&other.z);
        let z1z2 = self.z.mul(&other.z);

        let u = other.y.mul(&self.z).sub(&y1z2);
        let uu = u.square();

        let v = other.x.mul(&self.z).sub(&x1z2);
        let vv = v.square();

        let vvv = v.mul(&vv);
        let r = vv.mul(&x1z2);

        let two = Mnt4Fq2::from_u64(2, self.params.fq2_params.clone());

        let a = uu.mul(&z1z2).sub(&vvv).sub(&r.mul(&two));

        let x3 = v.mul(&a);

        let y3 = u.mul(&r.sub(&a)).sub(&vvv.mul(&y1z2));

        let z3 = vvv.mul(&z1z2);

        Self::new(x3, y3, z3, self.params.clone())
    }

    // ------------------------------------------------------------------------
    // Mixed addition
    // ------------------------------------------------------------------------

    /// Mixed addition.
    ///
    /// The second point must be affine (Z = 1), or infinity.
    pub fn mixed_add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        assert!(
            other.is_special(),
            "mixed_add requires the second point to be affine"
        );

        if self.equals(other) {
            return self.double();
        }

        if self.x == other.x.mul(&self.z) && self.y == other.y.negate().mul(&self.z) {
            return Self::zero(self.params.clone());
        }

        let u = self.z.mul(&other.y).sub(&self.y);
        let uu = u.square();

        let v = self.z.mul(&other.x).sub(&self.x);
        let vv = v.square();

        let vvv = v.mul(&vv);
        let r = vv.mul(&self.x);

        let two = Mnt4Fq2::from_u64(2, self.params.fq2_params.clone());

        let a = uu.mul(&self.z).sub(&vvv).sub(&r.mul(&two));

        let x3 = v.mul(&a);

        let y3 = u.mul(&r.sub(&a)).sub(&vvv.mul(&self.y));

        let z3 = vvv.mul(&self.z);

        Self::new(x3, y3, z3, self.params.clone())
    }

    // ------------------------------------------------------------------------
    // Scalar multiplication
    // ------------------------------------------------------------------------

    /// Multiply by a little-endian integer.
    ///
    /// scalar[0] contains the least significant 64 bits.
    pub fn scalar_mul(&self, scalar: &[u64]) -> Self {
        let mut result = Self::zero(self.params.clone());
        let mut base = self.clone();

        for limb in scalar {
            let mut word = *limb;

            for _ in 0..64 {
                if (word & 1) != 0 {
                    result = result.add(&base);
                }

                word >>= 1;
                base = base.double();
            }
        }

        result
    }

    /// Multiply by the G2 cofactor.
    pub fn mul_by_cofactor(&self) -> Self {
        self.scalar_mul(&self.params.cofactor)
    }

    // ------------------------------------------------------------------------
    // Frobenius / multiplication by q
    // ------------------------------------------------------------------------

    /// Multiply the point by the base-field characteristic q.
    ///
    /// This is the twist Frobenius map used by libff.
    pub fn mul_by_q(&self) -> Self {
        if self.is_zero() {
            return Self::zero(self.params.clone());
        }

        let x_frob = self.x.frobenius_map(1);
        let y_frob = self.y.frobenius_map(1);
        let z_frob = self.z.frobenius_map(1);

        let x = x_frob.mul_by_fp(&self.params.twist_mul_by_q_x);

        let y = y_frob.mul_by_fp(&self.params.twist_mul_by_q_y);

        Self::new(x, y, z_frob, self.params.clone())
    }

    // ------------------------------------------------------------------------
    // Curve validation
    // ------------------------------------------------------------------------

    /// Check the twist curve equation.
    ///
    /// Z * (Y^2 - b*Z^2) = X * (X^2 + a*Z^2)
    pub fn is_well_formed(&self) -> bool {
        if self.is_zero() {
            return true;
        }

        let z2 = self.z.square();

        let y2 = self.y.square();
        let lhs = self
            .y
            .square()
            .sub(&self.params.twist_coeff_b.mul(&z2))
            .mul(&self.z);

        let x2 = self.x.square();

        let rhs = self.x.mul(&x2.add(&self.params.twist_coeff_a.mul(&z2)));

        let _ = y2;

        lhs == rhs
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_g2_zero() {
        let params = mnt4_g2_params();
        let zero = Mnt4G2::zero(params);

        assert!(zero.is_zero());
        assert!(zero.is_well_formed());
    }

    #[test]
    fn test_g2_generator_is_well_formed() {
        let g = Mnt4G2::generator();

        assert!(!g.is_zero());
        assert!(g.is_well_formed());
    }

    #[test]
    fn test_g2_add_zero() {
        let g = Mnt4G2::generator();
        let zero = Mnt4G2::zero(g.params.clone());

        assert_eq!(g.add(&zero), g);
        assert_eq!(zero.add(&g), g);
    }

    #[test]
    fn test_g2_add_inverse() {
        let g = Mnt4G2::generator();
        let neg_g = g.negate();

        let result = g.add(&neg_g);

        assert!(result.is_zero());
    }

    #[test]
    fn test_g2_double_equals_add_self() {
        let g = Mnt4G2::generator();

        assert_eq!(g.double(), g.add(&g));
    }

    #[test]
    fn test_g2_affine_conversion() {
        let g = Mnt4G2::generator();

        let affine = g.to_affine_element();

        assert!(affine.is_special());
        assert!(affine.is_well_formed());
        assert_eq!(g, affine);
    }

    #[test]
    fn test_g2_scalar_mul_zero() {
        let g = Mnt4G2::generator();

        let result = g.scalar_mul(&[0, 0, 0, 0, 0]);

        assert!(result.is_zero());
    }

    #[test]
    fn test_g2_scalar_mul_one() {
        let g = Mnt4G2::generator();

        let result = g.scalar_mul(&[1, 0, 0, 0, 0]);

        assert_eq!(result, g);
    }

    #[test]
    fn test_g2_scalar_mul_two() {
        let g = Mnt4G2::generator();

        let result = g.scalar_mul(&[2, 0, 0, 0, 0]);

        assert_eq!(result, g.double());
    }

    #[test]
    fn test_g2_mixed_add() {
        let g = Mnt4G2::generator();

        let affine = g.to_affine_element();

        let result = g.mixed_add(&affine);

        assert_eq!(result, g.double());
    }
}
