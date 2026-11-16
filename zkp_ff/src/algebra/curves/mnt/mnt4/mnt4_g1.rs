use std::sync::Arc;

use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::{FieldParams, LimbsN};

/// MNT4-298 base field.
pub type Mnt4Fq = Fp<5>;

/// MNT4-298 G1 curve parameters.
#[derive(Clone, Debug)]
pub struct Mnt4G1Params {
    pub fq_params: FieldParams<5>,
    pub coeff_a: Mnt4Fq,
    pub coeff_b: Mnt4Fq,
    pub generator_x: Mnt4Fq,
    pub generator_y: Mnt4Fq,
}

impl Mnt4G1Params {
    /// Initialize MNT4-298 G1 parameters.
    ///
    /// The supplied FieldParams must use the MNT4-298 Fq modulus.
    pub fn new(fq_params: FieldParams<5>) -> Self {
        const B: LimbsN<5> = [
            0x5d4b_0ef5_28ec_0fd5,
            0x7b98_0f4e_9cd2_1a51,
            0xc3ed_d2a2_070a_085c,
            0x1558_5ea4_d523_234f,
            0x0000_0354_5a27_6394,
        ];

        const GENERATOR_X: LimbsN<5> = [
            0xa216_b5b9_c8d2_af46,
            0x60c4_0794_92b9_48de,
            0xaee8_6aba_8f73_d690,
            0xba85_213f_e6ca_3875,
            0x0000_007a_2caf_82a1,
        ];

        const GENERATOR_Y: LimbsN<5> = [
            0x9cba_48e7_10a4_8ab2,
            0x6778_c1af_d96a_71e2,
            0x89d0_148d_cc98_62d3,
            0x8267_2f7f_159f_ec2e,
            0x0000_02db_6194_61cc,
        ];

        let coeff_a = Mnt4Fq::from_u64(2, fq_params.clone());
        let coeff_b = Mnt4Fq::new(B, fq_params.clone());

        let generator_x = Mnt4Fq::new(GENERATOR_X, fq_params.clone());

        let generator_y = Mnt4Fq::new(GENERATOR_Y, fq_params.clone());

        Self {
            fq_params,
            coeff_a,
            coeff_b,
            generator_x,
            generator_y,
        }
    }

    /// Return the MNT4-298 G1 generator.
    pub fn generator(self: &Arc<Self>) -> Mnt4G1 {
        Mnt4G1 {
            x: self.generator_x.clone(),
            y: self.generator_y.clone(),
            z: Mnt4Fq::one(self.fq_params.clone()),
            params: Arc::clone(self),
        }
    }

    /// Return the point at infinity.
    pub fn zero(self: &Arc<Self>) -> Mnt4G1 {
        Mnt4G1 {
            x: Mnt4Fq::zero(self.fq_params.clone()),
            y: Mnt4Fq::one(self.fq_params.clone()),
            z: Mnt4Fq::zero(self.fq_params.clone()),
            params: Arc::clone(self),
        }
    }
}

/// MNT4-298 G1.
///
/// Projective coordinates:
///
/// ```text
/// x_affine = X / Z
/// y_affine = Y / Z
/// ```
///
/// This is NOT Jacobian coordinates.
#[derive(Clone, Debug)]
pub struct Mnt4G1 {
    pub x: Mnt4Fq,
    pub y: Mnt4Fq,
    pub z: Mnt4Fq,

    pub params: Arc<Mnt4G1Params>,
}

impl PartialEq for Mnt4G1 {
    fn eq(&self, other: &Self) -> bool {
        self.check_curve(other);

        if self.is_zero() {
            return other.is_zero();
        }

        if other.is_zero() {
            return false;
        }

        // X1 / Z1 == X2 / Z2
        if self.x.mul(&other.z) != other.x.mul(&self.z) {
            return false;
        }

        // Y1 / Z1 == Y2 / Z2
        self.y.mul(&other.z) == other.y.mul(&self.z)
    }
}

impl Eq for Mnt4G1 {}

impl Mnt4G1 {
    // ------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------

    /// Construct a projective point.
    pub fn new(x: Mnt4Fq, y: Mnt4Fq, z: Mnt4Fq, params: Arc<Mnt4G1Params>) -> Self {
        assert_eq!(
            x.params.modulus, params.fq_params.modulus,
            "X belongs to a different field"
        );
        assert_eq!(
            y.params.modulus, params.fq_params.modulus,
            "Y belongs to a different field"
        );
        assert_eq!(
            z.params.modulus, params.fq_params.modulus,
            "Z belongs to a different field"
        );

        Self { x, y, z, params }
    }

    /// Point at infinity.
    ///
    /// O = (0, 1, 0)
    pub fn zero(params: Arc<Mnt4G1Params>) -> Self {
        params.zero()
    }

    /// Standard G1 generator.
    pub fn one(params: Arc<Mnt4G1Params>) -> Self {
        params.generator()
    }

    fn field_zero(&self) -> Mnt4Fq {
        Mnt4Fq::zero(self.params.fq_params.clone())
    }

    fn field_one(&self) -> Mnt4Fq {
        Mnt4Fq::one(self.params.fq_params.clone())
    }

    fn check_curve(&self, other: &Self) {
        assert_eq!(
            self.params.fq_params.modulus, other.params.fq_params.modulus,
            "Points belong to different base fields"
        );

        assert!(
            self.params.coeff_a == other.params.coeff_a
                && self.params.coeff_b == other.params.coeff_b,
            "Points belong to different curves"
        );
    }

    // ------------------------------------------------------------------
    // Basic predicates
    // ------------------------------------------------------------------

    /// Match libff:
    ///```text
    ///     X == 0 && Z == 0
    ///```
    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.z.is_zero()
    }

    /// A point is special if it is infinity or Z == 1.
    pub fn is_special(&self) -> bool {
        self.is_zero() || self.z == self.field_one()
    }

    // ------------------------------------------------------------------
    // Affine conversion
    // ------------------------------------------------------------------

    /// Convert to affine coordinates in place.
    ///
    /// Returns false if Z has no inverse.
    pub fn to_affine(&mut self) -> bool {
        if self.is_zero() {
            self.x = self.field_zero();
            self.y = self.field_one();
            self.z = self.field_zero();
            return true;
        }

        let Some(z_inv) = self.z.inverse() else {
            return false;
        };

        self.x = self.x.mul(&z_inv);
        self.y = self.y.mul(&z_inv);
        self.z = self.field_one();

        true
    }

    /// Return an affine copy.
    pub fn to_affine_point(&self) -> Option<Self> {
        let mut point = self.clone();

        if point.to_affine() { Some(point) } else { None }
    }

    /// Convert to special representation.
    pub fn to_special(&mut self) -> bool {
        self.to_affine()
    }

    // ------------------------------------------------------------------
    // Curve validation
    // ------------------------------------------------------------------

    /// Check:
    ///```text
    ///     Z * (Y² - b*Z²) == X * (X² + a*Z²)
    ///```
    /// This is the projective form of:
    ////// The curve equation is:
    ///
    /// ```text
    /// y^2 = x^3 + a*x + b
    /// ```

    pub fn is_well_formed(&self) -> bool {
        if self.is_zero() {
            return true;
        }

        let x2 = self.x.square();
        let y2 = self.y.square();
        let z2 = self.z.square();

        let lhs = self.z.mul(&y2.sub(&self.params.coeff_b.mul(&z2)));

        let rhs = self.x.mul(&x2.add(&self.params.coeff_a.mul(&z2)));

        lhs == rhs
    }

    // ------------------------------------------------------------------
    // Negation
    // ------------------------------------------------------------------

    /// Return -P = (X, -Y, Z).
    pub fn negate(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.negate(),
            z: self.z.clone(),
            params: Arc::clone(&self.params),
        }
    }

    // ------------------------------------------------------------------
    // Point doubling
    // ------------------------------------------------------------------

    /// Double a point.
    ///
    /// Corresponds to libff mnt4_G1::dbl().
    pub fn double(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }

        let xx = self.x.square();
        let zz = self.z.square();

        // w = a*ZZ + 3*XX
        let three_xx = xx.add(&xx).add(&xx);

        let w = self.params.coeff_a.mul(&zz).add(&three_xx);

        // s = 2*Y*Z
        let y1z1 = self.y.mul(&self.z);
        let s = y1z1.add(&y1z1);

        let ss = s.square();
        let sss = s.mul(&ss);

        // R = Y*s
        let r = self.y.mul(&s);
        let rr = r.square();

        // B = (X+R)^2 - XX - RR
        let b = self.x.add(&r).square().sub(&xx).sub(&rr);

        // h = w² - 2*B
        let h = w.square().sub(&b.add(&b));

        // X3 = h*s
        let x3 = h.mul(&s);

        // Y3 = w*(B-h) - 2*RR
        let y3 = w.mul(&b.sub(&h)).sub(&rr.add(&rr));

        // Z3 = sss
        let z3 = sss;

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: Arc::clone(&self.params),
        }
    }

    // ------------------------------------------------------------------
    // General point addition
    // ------------------------------------------------------------------

    /// Add two projective points.
    ///
    /// Corresponds to libff mnt4_G1::add().
    pub fn add(&self, other: &Self) -> Self {
        self.check_curve(other);

        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        // Equal points: use doubling.
        if self == other {
            return self.double();
        }

        // Y1Z2 = Y1 * Z2
        let y1z2 = self.y.mul(&other.z);

        // X1Z2 = X1 * Z2
        let x1z2 = self.x.mul(&other.z);

        // Z1Z2 = Z1 * Z2
        let z1z2 = self.z.mul(&other.z);

        // u = Y2*Z1 - Y1*Z2
        let u = other.y.mul(&self.z).sub(&y1z2);

        let uu = u.square();

        // v = X2*Z1 - X1*Z2
        let v = other.x.mul(&self.z).sub(&x1z2);

        let vv = v.square();
        let vvv = v.mul(&vv);

        // R = vv * X1Z2
        let r = vv.mul(&x1z2);

        // A = uu*Z1Z2 - vvv - 2*R
        let a = uu.mul(&z1z2).sub(&vvv).sub(&r).sub(&r);

        // X3 = v*A
        let x3 = v.mul(&a);

        // Y3 = u*(R-A) - vvv*Y1Z2
        let y3 = u.mul(&r.sub(&a)).sub(&vvv.mul(&y1z2));

        // Z3 = vvv*Z1Z2
        let z3 = vvv.mul(&z1z2);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: Arc::clone(&self.params),
        }
    }

    // ------------------------------------------------------------------
    // Mixed addition
    // ------------------------------------------------------------------

    /// Add a point whose coordinates are affine (Z == 1),
    /// or infinity.
    ///
    /// Corresponds to libff mnt4_G1::mixed_add().
    pub fn mixed_add(&self, other: &Self) -> Self {
        self.check_curve(other);

        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        debug_assert!(
            other.is_special(),
            "mixed_add expects an affine/special second operand"
        );

        // Since other.Z == 1:
        //
        // X1Z2 = X1
        // Y1Z2 = Y1
        let x1z2 = self.x.clone();
        let y1z2 = self.y.clone();

        let x2z1 = self.z.mul(&other.x);
        let y2z1 = self.z.mul(&other.y);

        // Equal points: double.
        if x1z2 == x2z1 && y1z2 == y2z1 {
            return self.double();
        }

        // u = Y2*Z1 - Y1
        let u = y2z1.sub(&self.y);
        let uu = u.square();

        // v = X2*Z1 - X1
        let v = x2z1.sub(&self.x);
        let vv = v.square();
        let vvv = v.mul(&vv);

        // R = vv*X1
        let r = vv.mul(&self.x);

        // A = uu*Z1 - vvv - 2*R
        let a = uu.mul(&self.z).sub(&vvv).sub(&r).sub(&r);

        // X3 = v*A
        let x3 = v.mul(&a);

        // Y3 = u*(R-A) - vvv*Y1
        let y3 = u.mul(&r.sub(&a)).sub(&vvv.mul(&self.y));

        // Z3 = vvv*Z1
        let z3 = vvv.mul(&self.z);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: Arc::clone(&self.params),
        }
    }

    // ------------------------------------------------------------------
    // Subtraction
    // ------------------------------------------------------------------

    /// P - Q = P + (-Q).
    pub fn sub(&self, other: &Self) -> Self {
        self.check_curve(other);
        self.add(&other.negate())
    }

    // ------------------------------------------------------------------
    // Cofactor
    // ------------------------------------------------------------------

    /// MNT4-298 G1 cofactor is 1.
    pub fn mul_by_cofactor(&self) -> Self {
        self.clone()
    }

    // ------------------------------------------------------------------
    // Scalar multiplication
    // ------------------------------------------------------------------

    /// Multiply by a little-endian 5-limb integer.
    ///
    /// This is a straightforward double-and-add implementation.
    /// The scalar is interpreted as an unsigned integer, not reduced
    /// modulo the scalar field.
    pub fn mul_scalar(&self, scalar: LimbsN<5>) -> Self {
        let mut result = Self::zero(Arc::clone(&self.params));
        let mut base = self.clone();

        for limb in scalar.iter() {
            let mut value = *limb;

            for _ in 0..64 {
                if value & 1 == 1 {
                    result = result.add(&base);
                }

                base = base.double();
                value >>= 1;
            }
        }

        result
    }

    // ------------------------------------------------------------------
    // Batch affine conversion
    // ------------------------------------------------------------------

    /// Convert all non-infinity points to affine coordinates.
    ///
    /// Uses one field inversion for the entire batch.
    ///
    /// Infinity points are preserved unchanged.
    ///
    /// Returns false if a non-infinity point has a non-invertible Z.
    pub fn batch_to_special_all_non_zeros(points: &mut [Self]) -> bool {
        if points.is_empty() {
            return true;
        }

        let mut indices = Vec::new();
        let mut prefix_products: Vec<Mnt4Fq> = Vec::new();

        let mut product: Option<Mnt4Fq> = None;

        for (index, point) in points.iter().enumerate() {
            if point.is_zero() {
                continue;
            }

            if point.z.is_zero() {
                return false;
            }

            indices.push(index);

            // Store product of all previous non-zero Z values.
            prefix_products.push(match &product {
                Some(value) => value.clone(),
                None => point.field_one(),
            });

            product = Some(match product {
                Some(value) => value.mul(&point.z),
                None => point.z.clone(),
            });
        }

        if indices.is_empty() {
            return true;
        }

        let Some(mut inverse_product) = product
            .expect("non-empty indices imply a product")
            .inverse()
        else {
            return false;
        };

        // Walk backwards:
        //
        // inv(Z_i) = inverse(product(Z_0..Z_n)) *
        //            product(Z_0..Z_(i-1)) *
        //            product(Z_(i+1)..Z_n)
        for position in (0..indices.len()).rev() {
            let index = indices[position];

            let z_inv = inverse_product.mul(&prefix_products[position]);

            // Remove this Z from the accumulated inverse product.
            inverse_product = inverse_product.mul(&points[index].z);

            points[index].x = points[index].x.mul(&z_inv);
            points[index].y = points[index].y.mul(&z_inv);
            points[index].z = points[index].field_one();
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> Arc<Mnt4G1Params> {
        let fq_modulus: LimbsN<5> = [
            0xc90c_d65a_7166_0001,
            0x41a9_e35e_5120_0e12,
            0xcaee_c963_5d13_30ea,
            0xa266_249d_a7b0_548e,
            0x0000_03bc_f7bc_d473,
        ];

        Arc::new(Mnt4G1Params::new(FieldParams::<5>::new(fq_modulus)))
    }

    #[test]
    fn test_generator_is_well_formed() {
        let params = setup();
        let g = Mnt4G1::one(params);

        assert!(g.is_well_formed());
        assert!(!g.is_zero());
    }

    #[test]
    fn test_add_zero() {
        let params = setup();

        let g = Mnt4G1::one(Arc::clone(&params));
        let zero = Mnt4G1::zero(params);

        assert_eq!(g.add(&zero), g);
    }

    #[test]
    fn test_double_equals_add_self() {
        let params = setup();

        let g = Mnt4G1::one(params);

        assert_eq!(g.double(), g.add(&g));
    }

    #[test]
    fn test_add_inverse() {
        let params = setup();

        let g = Mnt4G1::one(params);
        let neg_g = g.negate();

        let result = g.add(&neg_g);

        assert!(result.is_zero());
    }

    #[test]
    fn test_scalar_multiplication() {
        let params = setup();

        let g = Mnt4G1::one(params);

        let two = g.mul_scalar([2, 0, 0, 0, 0]);
        let three = g.mul_scalar([3, 0, 0, 0, 0]);
        let five = g.mul_scalar([5, 0, 0, 0, 0]);

        assert_eq!(two.add(&three), five);
        assert!(five.is_well_formed());
    }

    #[test]
    fn test_affine_conversion() {
        let params = setup();

        let g = Mnt4G1::one(params);

        // 2G 通常不是 Z=1 的表示。
        let mut p = g.double();

        assert!(p.is_well_formed());
        assert!(p.to_affine());
        assert!(p.is_special());
        assert!(p.is_well_formed());
    }
}
