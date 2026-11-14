use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::Limbs;
use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};

/// BN254 G2 参数。
///
/// G2 使用 Jacobian 坐标：
/// ```text
///     x = X / Z²
///     y = Y / Z³
/// ```
/// 曲线方程：
/// ```text
///     Y² = X³ + b * Z⁶
/// ```
#[derive(Clone, Debug)]
pub struct G2Params {
    pub fp2_params: Fp2Params,
    pub b: Fp2,

    /// G2 generator coordinates.
    pub generator_x: Fp2,
    pub generator_y: Fp2,

    /// q-Frobenius twist coefficients.
    pub twist_mul_by_q_x: Fp2,
    pub twist_mul_by_q_y: Fp2,

    /// Scalar field modulus r.
    pub scalar_modulus: [u64; 4],

    /// G2 cofactor.
    pub cofactor: [u64; 4],
}

/// BN254 G2 point in Jacobian coordinates.
#[derive(Clone, Debug)]
pub struct G2 {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
    pub params: G2Params,
}

impl G2Params {
    pub fn new(
        fp2_params: Fp2Params,
        b: Fp2,
        generator_x: Fp2,
        generator_y: Fp2,
        twist_mul_by_q_x: Fp2,
        twist_mul_by_q_y: Fp2,
        scalar_modulus: [u64; 4],
        cofactor: [u64; 4],
    ) -> Self {
        Self {
            fp2_params,
            b,
            generator_x,
            generator_y,
            twist_mul_by_q_x,
            twist_mul_by_q_y,
            scalar_modulus,
            cofactor,
        }
    }

    pub fn zero(&self) -> G2 {
        G2 {
            x: Fp2::zero(self.fp2_params.clone()),
            y: Fp2::one(self.fp2_params.clone()),
            z: Fp2::zero(self.fp2_params.clone()),
            params: self.clone(),
        }
    }

    pub fn one(&self) -> G2 {
        G2 {
            x: self.generator_x.clone(),
            y: self.generator_y.clone(),
            z: Fp2::one(self.fp2_params.clone()),
            params: self.clone(),
        }
    }
}

impl PartialEq for G2 {
    fn eq(&self, other: &Self) -> bool {
        if self.params.fp2_params.fp_params.modulus != other.params.fp2_params.fp_params.modulus {
            return false;
        }

        if self.is_zero() || other.is_zero() {
            return self.is_zero() && other.is_zero();
        }

        let z1_squared = self.z.square();
        let z2_squared = other.z.square();

        let u1 = self.x.mul(&z2_squared);
        let u2 = other.x.mul(&z1_squared);

        if u1 != u2 {
            return false;
        }

        let z1_cubed = self.z.mul(&z1_squared);
        let z2_cubed = other.z.mul(&z2_squared);

        let s1 = self.y.mul(&z2_cubed);
        let s2 = other.y.mul(&z1_cubed);

        s1 == s2
    }
}

impl Eq for G2 {}

impl G2 {
    pub fn new(x: Fp2, y: Fp2, z: Fp2, params: G2Params) -> Self {
        Self { x, y, z, params }
    }

    pub fn zero(params: G2Params) -> Self {
        params.zero()
    }

    pub fn one(params: G2Params) -> Self {
        params.one()
    }

    pub fn is_zero(&self) -> bool {
        self.z.is_zero()
    }

    /// Check whether Z = 1.
    pub fn is_special(&self) -> bool {
        self.z.is_one()
    }

    fn check_field(&self, other: &Self) {
        assert_eq!(
            self.params.fp2_params.fp_params.modulus, other.params.fp2_params.fp_params.modulus,
            "G2 points belong to different base fields"
        );

        assert_eq!(
            self.params.fp2_params.non_residue, other.params.fp2_params.non_residue,
            "G2 points belong to different extension fields"
        );
    }

    fn two(&self) -> Fp2 {
        Fp2::from_u64(2, self.params.fp2_params.clone())
    }

    fn three(&self) -> Fp2 {
        Fp2::from_u64(3, self.params.fp2_params.clone())
    }

    fn four(&self) -> Fp2 {
        Fp2::from_u64(4, self.params.fp2_params.clone())
    }

    fn eight(&self) -> Fp2 {
        Fp2::from_u64(8, self.params.fp2_params.clone())
    }

    // ------------------------------------------------------------
    // Coordinate conversion
    // ------------------------------------------------------------

    /// Convert this point to affine coordinates in place.
    ///
    /// Returns false for the point at infinity.
    pub fn to_affine_coordinates(&mut self) -> bool {
        if self.is_zero() {
            return false;
        }

        if self.is_special() {
            return true;
        }

        let Some(z_inv) = self.z.inverse() else {
            return false;
        };

        let z_inv2 = z_inv.square();
        let z_inv3 = z_inv2.mul(&z_inv);

        self.x = self.x.mul(&z_inv2);
        self.y = self.y.mul(&z_inv3);
        self.z = Fp2::one(self.params.fp2_params.clone());

        true
    }

    /// Return an affine copy.
    pub fn to_affine(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let mut result = self.clone();

        if !result.to_affine_coordinates() {
            return None;
        }

        Some(result)
    }

    /// Convert to the special representation.
    pub fn to_special(&mut self) {
        self.to_affine_coordinates();
    }

    // ------------------------------------------------------------
    // Negation
    // ------------------------------------------------------------

    pub fn negate(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }

        Self {
            x: self.x.clone(),
            y: self.y.negate(),
            z: self.z.clone(),
            params: self.params.clone(),
        }
    }

    // ------------------------------------------------------------
    // Point doubling
    // ------------------------------------------------------------

    /// Jacobian point doubling.
    ///
    /// A = X²
    /// B = Y²
    /// C = B²
    /// D = 2 * ((X+B)² - A - C)
    /// E = 3*A
    /// F = E²
    ///
    /// X3 = F - 2*D
    /// Y3 = E*(D-X3) - 8*C
    /// Z3 = 2*Y*Z
    pub fn dbl(&self) -> Self {
        if self.is_zero() || self.y.is_zero() {
            return self.params.zero();
        }

        let a = self.x.square();
        let b = self.y.square();
        let c = b.square();

        let x_plus_b = self.x.add(&b);

        let d = x_plus_b.square().sub(&a).sub(&c).mul(&self.two());

        let e = a.mul(&self.three());
        let f = e.square();

        let x3 = f.sub(&d.mul(&self.two()));

        let y3 = e.mul(&d.sub(&x3)).sub(&c.mul(&self.eight()));

        let z3 = self.y.mul(&self.z).mul(&self.two());

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: self.params.clone(),
        }
    }

    // ------------------------------------------------------------
    // Full Jacobian point addition
    // ------------------------------------------------------------

    pub fn add(&self, other: &Self) -> Self {
        self.check_field(other);

        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        let z1z1 = self.z.square();
        let z2z2 = other.z.square();

        let u1 = self.x.mul(&z2z2);
        let u2 = other.x.mul(&z1z1);

        let z1_cubed = self.z.mul(&z1z1);
        let z2_cubed = other.z.mul(&z2z2);

        let s1 = self.y.mul(&z2_cubed);
        let s2 = other.y.mul(&z1_cubed);

        if u1 == u2 {
            if s1 == s2 {
                return self.dbl();
            }

            return self.params.zero();
        }

        let h = u2.sub(&u1);
        let s2_minus_s1 = s2.sub(&s1);

        let two_h = h.mul(&self.two());
        let i = two_h.square();

        let j = h.mul(&i);
        let r = s2_minus_s1.mul(&self.two());
        let v = u1.mul(&i);

        let x3 = r.square().sub(&j).sub(&v.mul(&self.two()));

        let y3 = r.mul(&v.sub(&x3)).sub(&s1.mul(&j).mul(&self.two()));

        let z3 = self.z.add(&other.z).square().sub(&z1z1).sub(&z2z2).mul(&h);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: self.params.clone(),
        }
    }

    // ------------------------------------------------------------
    // Mixed addition
    // ------------------------------------------------------------

    /// Mixed addition.
    ///
    /// `other` must be an affine/special point (Z = 1).
    pub fn mixed_add(&self, other: &Self) -> Self {
        self.check_field(other);

        assert!(
            other.is_zero() || other.is_special(),
            "mixed_add requires an affine second operand"
        );

        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        let z1z1 = self.z.square();

        let u1 = self.x.clone();
        let u2 = other.x.mul(&z1z1);

        let z1_cubed = self.z.mul(&z1z1);

        let s1 = self.y.clone();
        let s2 = other.y.mul(&z1_cubed);

        if u1 == u2 {
            if s1 == s2 {
                return self.dbl();
            }

            return self.params.zero();
        }

        let h = u2.sub(&self.x);
        let hh = h.square();

        let i = hh.mul(&self.four());
        let j = h.mul(&i);

        let r = s2.sub(&self.y).mul(&self.two());
        let v = self.x.mul(&i);

        let x3 = r.square().sub(&j).sub(&v.mul(&self.two()));

        let y3 = r.mul(&v.sub(&x3)).sub(&self.y.mul(&j).mul(&self.two()));

        let z3 = self.z.add(&h).square().sub(&z1z1).sub(&hh);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: self.params.clone(),
        }
    }

    // ------------------------------------------------------------
    // Scalar multiplication
    // ------------------------------------------------------------

    /// Multiply by a little-endian 256-bit scalar.
    ///
    /// This is the straightforward double-and-add implementation.
    pub fn scalar_mul(&self, scalar: &[u64; 4]) -> Self {
        let mut result = self.params.zero();
        let mut base = self.clone();

        for limb in scalar.iter() {
            let mut word = *limb;

            for _ in 0..64 {
                if word & 1 == 1 {
                    result = result.add(&base);
                }

                word >>= 1;
                base = base.dbl();
            }
        }

        result
    }

    /// Multiply by the base-field modulus q using the twist
    /// Frobenius coefficients.
    ///
    /// For a point on the BN254 twist:
    ///
    /// X' = Frobenius(X) * twist_mul_by_q_x
    /// Y' = Frobenius(Y) * twist_mul_by_q_y
    /// Z' = Frobenius(Z)
    pub fn mul_by_q(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }

        Self {
            x: self.x.frobenius_map(1).mul(&self.params.twist_mul_by_q_x),

            y: self.y.frobenius_map(1).mul(&self.params.twist_mul_by_q_y),

            z: self.z.frobenius_map(1),

            params: self.params.clone(),
        }
    }

    /// Multiply by the G2 cofactor.
    pub fn mul_by_cofactor(&self) -> Self {
        self.scalar_mul(&self.params.cofactor)
    }

    /// Check membership in the prime-order subgroup.
    ///
    /// P belongs to G2[r] iff [r]P = O.
    pub fn is_in_safe_subgroup(&self) -> bool {
        self.is_zero() || self.scalar_mul(&self.params.scalar_modulus).is_zero()
    }

    // ------------------------------------------------------------
    // Curve validation
    // ------------------------------------------------------------

    /// Check:
    ///
    /// Y² = X³ + b*Z⁶
    pub fn is_well_formed(&self) -> bool {
        if self.is_zero() {
            return true;
        }

        let y2 = self.y.square();
        let x3 = self.x.square().mul(&self.x);

        let z2 = self.z.square();
        let z6 = z2.square().mul(&z2);

        let rhs = x3.add(&self.params.b.mul(&z6));

        y2 == rhs
    }

    // ------------------------------------------------------------
    // Batch normalization
    // ------------------------------------------------------------

    /// Normalize all nonzero points using one field inversion.
    ///
    /// Infinity points are preserved.
    pub fn batch_to_special_all_non_zeros(points: &mut [Self]) {
        if points.is_empty() {
            return;
        }

        let mut prefixes: Vec<Fp2> = Vec::with_capacity(points.len());
        let mut product = Fp2::one(points[0].params.fp2_params.clone());

        for point in points.iter() {
            prefixes.push(product.clone());

            if !point.is_zero() {
                product = product.mul(&point.z);
            }
        }

        let Some(mut inverse_product) = product.inverse() else {
            return;
        };

        for i in (0..points.len()).rev() {
            if points[i].is_zero() {
                continue;
            }

            let z_inverse = inverse_product.mul(&prefixes[i]);

            inverse_product = inverse_product.mul(&points[i].z);

            let z_inverse2 = z_inverse.square();
            let z_inverse3 = z_inverse2.mul(&z_inverse);

            points[i].x = points[i].x.mul(&z_inverse2);
            points[i].y = points[i].y.mul(&z_inverse3);
            points[i].z = Fp2::one(points[i].params.fp2_params.clone());
        }
    }

    // ------------------------------------------------------------
    // Compressed serialization
    // ------------------------------------------------------------

    /// Compressed encoding:
    ///
    /// [flag, X.c0 (4 words), X.c1 (4 words), sign]
    ///
    /// flag:
    ///   0 = infinity
    ///   1 = ordinary point
    ///
    /// sign is the least significant bit of Y.c0.
    pub fn to_compressed_words(&self) -> Vec<u64> {
        if self.is_zero() {
            return vec![0u64; 10];
        }

        let affine = self
            .to_affine()
            .expect("nonzero point must have affine coordinates");

        let mut words = Vec::with_capacity(10);

        words.push(1);
        words.extend_from_slice(&affine.x.to_words());

        let y_words = affine.y.to_words();
        let sign = y_words[0] & 1;

        words.push(sign);

        words
    }

    /// Decode a compressed point.
    pub fn from_compressed_words(words: &[u64], params: G2Params) -> Option<Self> {
        if words.len() != 10 {
            return None;
        }

        if words[0] == 0 {
            return Some(params.zero());
        }

        if words[0] != 1 || words[9] > 1 {
            return None;
        }

        let x = Fp2::from_words(&words[1..9], params.fp2_params.clone())?;

        let x3 = x.square().mul(&x);
        let rhs = x3.add(&params.b);

        let mut y = rhs.sqrt()?;

        let y_sign = y.to_words()[0] & 1;

        if y_sign != words[9] {
            y = y.negate();
        }

        let point = Self {
            x,
            y,
            z: Fp2::one(params.fp2_params.clone()),
            params,
        };

        if !point.is_well_formed() {
            return None;
        }

        Some(point)
    }

    // ------------------------------------------------------------
    // Uncompressed serialization
    // ------------------------------------------------------------

    /// Uncompressed encoding:
    ///
    /// [flag, X (8 words), Y (8 words)]
    pub fn to_uncompressed_words(&self) -> Vec<u64> {
        if self.is_zero() {
            return vec![0u64; 17];
        }

        let affine = self
            .to_affine()
            .expect("nonzero point must have affine coordinates");

        let mut words = Vec::with_capacity(17);

        words.push(1);
        words.extend_from_slice(&affine.x.to_words());
        words.extend_from_slice(&affine.y.to_words());

        words
    }

    pub fn from_uncompressed_words(words: &[u64], params: G2Params) -> Option<Self> {
        if words.len() != 17 {
            return None;
        }

        if words[0] == 0 {
            return Some(params.zero());
        }

        if words[0] != 1 {
            return None;
        }

        let x = Fp2::from_words(&words[1..9], params.fp2_params.clone())?;

        let y = Fp2::from_words(&words[9..17], params.fp2_params.clone())?;

        let point = Self {
            x,
            y,
            z: Fp2::one(params.fp2_params.clone()),
            params,
        };

        if !point.is_well_formed() {
            return None;
        }

        Some(point)
    }
}

// ============================================================================
// BN254 G2 parameters
// ============================================================================

/// BN254 G2 generator and parameters.
///  ```text
/// Fq2 = Fq[u] / (u² + 1)
///
/// Twist coefficient:
///
///     b = 3 / (9 + u)
///
/// Fp2 coefficient order:
///
///     c0 = real component
///     c1 = imaginary component
/// ```
pub fn bn254_g2_params() -> G2Params {
    let fp_params = crate::algebra::fields::prime_base::fp_montgomery::FieldParams::new([
        0x3c20_8c16_d87c_fd47,
        0x9781_6a91_6871_ca8d,
        0xb850_45b6_8181_585d,
        0x3064_4e72_e131_a029,
    ]);

    let one_fp = Fp::one(fp_params.clone());
    let minus_one_fp = one_fp.negate();

    let fp2_params = Fp2Params::new(
        fp_params.clone(),
        minus_one_fp.clone(),
        [one_fp, minus_one_fp],
    );

    // ------------------------------------------------------------
    // Twist coefficient b = 3 / (9 + u)
    // ------------------------------------------------------------

    let nine = Fp::from_u64(9, fp_params.clone());
    let three = Fp::from_u64(3, fp_params.clone());
    let one = Fp::one(fp_params.clone());

    let twist_denominator = Fp2::new(nine, one, fp2_params.clone());

    let twist_b = Fp2::new(three, Fp::zero(fp_params.clone()), fp2_params.clone());

    let b = twist_b.mul(
        &twist_denominator
            .inverse()
            .expect("9 + u must be invertible"),
    );

    // ------------------------------------------------------------
    // G2 generator X
    // ------------------------------------------------------------

    let gx0 = Fp::from_words(
        &[
            0x46de_bd5c_d992_f6ed,
            0x6743_22d4_f75e_dadd,
            0x426a_0066_5e5c_4479,
            0x1800_deef_121f_1e76,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 G2 generator X.c0");

    let gx1 = Fp::from_words(
        &[
            0x97e4_85b7_aef3_12c2,
            0xf1aa_4933_35a9_e712,
            0x7260_bfb7_31fb_5d25,
            0x198e_9393_920d_483a,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 G2 generator X.c1");

    let generator_x = Fp2::new(gx0, gx1, fp2_params.clone());

    // ------------------------------------------------------------
    // G2 generator Y
    // ------------------------------------------------------------

    let gy0 = Fp::from_words(
        &[
            0x4ce6_cc01_66fa_7daa,
            0xe3d1_e769_0c43_d37b,
            0x4aab_7180_8dcb_408f,
            0x12c8_5ea5_db8c_6deb,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 G2 generator Y.c0");

    let gy1 = Fp::from_words(
        &[
            0x55ac_dadc_d122_975b,
            0xbc4b_3133_70b3_8ef3,
            0xec9e_99ad_690c_3395,
            0x0906_89d0_585f_f075,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 G2 generator Y.c1");

    let generator_y = Fp2::new(gy0, gy1, fp2_params.clone());

    // ------------------------------------------------------------
    // q-Frobenius twist coefficients
    //
    // Values are from libff alt_bn128_init.cpp.
    // ------------------------------------------------------------

    let qx0 = Fp::from_words(
        &[
            0x99e3_9557_176f_553d,
            0xb78c_c310_c2c3_330c,
            0x4c0b_ec3c_f559_b143,
            0x2fb3_4798_4f79_11f7,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 twist_mul_by_q_x.c0");

    let qx1 = Fp::from_words(
        &[
            0x1665_d51c_640f_cba2,
            0x32ae_2a1d_0b7c_9dce,
            0x4ba4_cc8b_d75a_0794,
            0x16c9_e550_61eb_ae20,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 twist_mul_by_q_x.c1");

    let twist_mul_by_q_x = Fp2::new(qx0, qx1, fp2_params.clone());

    let qy0 = Fp::from_words(
        &[
            0xdc54_0146_71a0_135a,
            0xdbaa_e0ed_a9c9_5998,
            0xdc5e_c698_b6e2_f9b9,
            0x063c_f305_489a_f5dc,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 twist_mul_by_q_y.c0");

    let qy1 = Fp::from_words(
        &[
            0x82d3_7f63_2623_b0e3,
            0x2180_7dc9_8fa2_5bd2,
            0x0704_b5a7_ec79_6f2b,
            0x07c0_3cbc_ac41_049a,
        ],
        fp_params.clone(),
    )
    .expect("valid BN254 twist_mul_by_q_y.c1");

    let twist_mul_by_q_y = Fp2::new(qy0, qy1, fp2_params.clone());

    // ------------------------------------------------------------
    // Scalar field modulus r
    // ------------------------------------------------------------

    let scalar_modulus: Limbs = [
        0x43e1_f593_f000_0001,
        0x2833_e848_79b9_7091,
        0xb850_45b6_8181_585d,
        0x3064_4e72_e131_a029,
    ];

    // ------------------------------------------------------------
    // BN254 G2 cofactor
    // ------------------------------------------------------------

    let cofactor: Limbs = [
        0x345f_2299_c0f9_fa8d,
        0x06ce_ecda_572a_2489,
        0xb850_45b6_8181_585e,
        0x3064_4e72_e131_a029,
    ];

    G2Params::new(
        fp2_params,
        b,
        generator_x,
        generator_y,
        twist_mul_by_q_x,
        twist_mul_by_q_y,
        scalar_modulus,
        cofactor,
    )
}

/// Return the BN254 G2 generator.
pub fn bn254_g2_generator() -> G2 {
    bn254_g2_params().one()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn generator() -> G2 {
        bn254_g2_generator()
    }

    #[test]
    fn test_g2_zero() {
        let params = bn254_g2_params();
        let zero = params.zero();

        assert!(zero.is_zero());
        assert!(zero.is_well_formed());
    }

    #[test]
    fn test_g2_generator() {
        let g = generator();

        assert!(!g.is_zero());
        assert!(g.is_special());
        assert!(g.is_well_formed());
    }

    #[test]
    fn test_g2_add_zero() {
        let g = generator();
        let zero = g.params.zero();

        assert_eq!(g.add(&zero), g);
        assert_eq!(zero.add(&g), g);
    }

    #[test]
    fn test_g2_negate() {
        let g = generator();
        let negative = g.negate();

        assert!(g.add(&negative).is_zero());
    }

    #[test]
    fn test_g2_double_equals_add() {
        let g = generator();

        assert_eq!(g.dbl(), g.add(&g));
    }

    #[test]
    fn test_g2_mixed_add() {
        let g = generator();

        let doubled1 = g.mixed_add(&g);
        let doubled2 = g.dbl();

        assert_eq!(doubled1, doubled2);
    }

    #[test]
    fn test_g2_add_associativity() {
        let g = generator();

        let two_g = g.dbl();
        let three_g = two_g.add(&g);

        let left = g.add(&g).add(&g);
        let right = g.add(&g.add(&g));

        assert_eq!(left, right);
        assert_eq!(left, three_g);
    }

    #[test]
    fn test_g2_scalar_zero() {
        let g = generator();

        let result = g.scalar_mul(&[0, 0, 0, 0]);

        assert!(result.is_zero());
    }

    #[test]
    fn test_g2_scalar_one() {
        let g = generator();

        let result = g.scalar_mul(&[1, 0, 0, 0]);

        assert_eq!(result, g);
    }

    #[test]
    fn test_g2_scalar_two() {
        let g = generator();

        let result = g.scalar_mul(&[2, 0, 0, 0]);

        assert_eq!(result, g.dbl());
    }

    #[test]
    fn test_g2_scalar_order() {
        let g = generator();

        let result = g.scalar_mul(&g.params.scalar_modulus);

        assert!(result.is_zero());
    }

    #[test]
    fn test_g2_subgroup() {
        let g = generator();

        assert!(g.is_in_safe_subgroup());
    }

    #[test]
    fn test_g2_mul_by_cofactor() {
        let g = generator();

        let result = g.mul_by_cofactor();

        assert!(result.is_well_formed());
    }

    #[test]
    fn test_g2_affine_conversion() {
        let g = generator();

        let two_g = g.dbl();
        let affine = two_g.to_affine().unwrap();

        assert!(affine.is_special());
        assert!(affine.is_well_formed());
        assert_eq!(affine, two_g);
    }

    #[test]
    fn test_g2_compressed_serialization() {
        let g = generator();

        let words = g.to_compressed_words();

        assert_eq!(words.len(), 10);

        let decoded = G2::from_compressed_words(&words, g.params.clone()).unwrap();

        assert_eq!(decoded, g);
    }

    #[test]
    fn test_g2_uncompressed_serialization() {
        let g = generator();

        let words = g.to_uncompressed_words();

        assert_eq!(words.len(), 17);

        let decoded = G2::from_uncompressed_words(&words, g.params.clone()).unwrap();

        assert_eq!(decoded, g);
    }

    #[test]
    fn test_g2_infinity_serialization() {
        let params = bn254_g2_params();
        let zero = params.zero();

        let compressed = zero.to_compressed_words();

        let decoded = G2::from_compressed_words(&compressed, params.clone()).unwrap();

        assert!(decoded.is_zero());

        let uncompressed = zero.to_uncompressed_words();

        let decoded = G2::from_uncompressed_words(&uncompressed, params).unwrap();

        assert!(decoded.is_zero());
    }

    #[test]
    fn test_g2_batch_normalization() {
        let g = generator();

        let two_g = g.dbl();
        let three_g = two_g.add(&g);

        let mut points = vec![g.clone(), two_g, three_g, g.params.zero()];

        G2::batch_to_special_all_non_zeros(&mut points);

        assert!(points[0].is_special());
        assert!(points[1].is_special());
        assert!(points[2].is_special());
        assert!(points[3].is_zero());

        assert!(points[0].is_well_formed());
        assert!(points[1].is_well_formed());
        assert!(points[2].is_well_formed());
    }

    #[test]
    fn test_generator_is_well_formed() {
        let g = bn254_g2_generator();

        assert!(g.is_well_formed());
        assert!(g.is_in_safe_subgroup());
    }

    #[test]
    fn test_generator_order() {
        let g = bn254_g2_generator();

        let r = g.params.scalar_modulus;

        assert!(g.scalar_mul(&r).is_zero());
    }

    #[test]
    fn test_generator_cofactor() {
        let g = bn254_g2_generator();

        let h = g.params.cofactor;
        let result = g.scalar_mul(&h);

        assert!(result.is_well_formed());
    }

    #[test]
    fn test_projective_equality() {
        let g = bn254_g2_generator();

        let lambda = Fp2::from_u64(7, g.params.fp2_params.clone());

        let lambda2 = lambda.square();
        let lambda3 = lambda2.mul(&lambda);

        let equivalent = G2 {
            x: g.x.mul(&lambda2),
            y: g.y.mul(&lambda3),
            z: lambda,
            params: g.params.clone(),
        };

        assert_eq!(g, equivalent);
    }

    #[test]
    fn test_g2_mul_by_q_is_well_formed() {
        let g = generator();

        let result = g.mul_by_q();

        assert!(result.is_well_formed());
    }

    #[test]
    fn test_g2_mul_by_q_matches_scalar_q() {
        let g = generator();

        let q: Limbs = [
            0x3c20_8c16_d87c_fd47,
            0x9781_6a91_6871_ca8d,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ];

        let expected = g.scalar_mul(&q);
        let actual = g.mul_by_q();

        assert_eq!(actual, expected);
    }
}
