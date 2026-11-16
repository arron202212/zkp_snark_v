use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

/// Fp2 parameters.
///
/// Fp2 = Fp[u] / (u^2 - non_residue)
///
/// An element is represented as c0 + c1 * u.
///
/// N is the number of u64 limbs in the base field.
///
/// The default N = 4 preserves compatibility with existing BN254 code.
#[derive(Clone, Debug)]
pub struct Fp2Params<const N: usize = 4> {
    /// Base-field parameters.
    pub fp_params: FieldParams<N>,

    /// u^2 = non_residue.
    pub non_residue: Fp<N>,

    /// Frobenius coefficients for the c1 component.
    ///
    /// For the quadratic extension:
    ///
    /// frobenius_coeffs_c1[0] = 1
    /// frobenius_coeffs_c1[1] = u^p / u
    ///
    /// The index is reduced modulo 2.
    pub frobenius_coeffs_c1: [Fp<N>; 2],
}

impl<const N: usize> Fp2Params<N> {
    pub fn new(
        fp_params: FieldParams<N>,
        non_residue: Fp<N>,
        frobenius_coeffs_c1: [Fp<N>; 2],
    ) -> Self {
        assert_eq!(
            fp_params.modulus, non_residue.params.modulus,
            "non_residue belongs to a different base field"
        );

        for coeff in &frobenius_coeffs_c1 {
            assert_eq!(
                fp_params.modulus, coeff.params.modulus,
                "Frobenius coefficient belongs to a different base field"
            );
        }

        Self {
            fp_params,
            non_residue,
            frobenius_coeffs_c1,
        }
    }

    pub fn zero(&self) -> Fp2<N> {
        Fp2 {
            c0: Fp::<N>::zero(self.fp_params.clone()),
            c1: Fp::<N>::zero(self.fp_params.clone()),
            params: self.clone(),
        }
    }

    pub fn one(&self) -> Fp2<N> {
        Fp2 {
            c0: Fp::<N>::one(self.fp_params.clone()),
            c1: Fp::<N>::zero(self.fp_params.clone()),
            params: self.clone(),
        }
    }

    /// Construct a base-field element embedded into Fp2.
    pub fn from_u64(&self, value: u64) -> Fp2<N> {
        Fp2 {
            c0: Fp::<N>::from_u64(value, self.fp_params.clone()),
            c1: Fp::<N>::zero(self.fp_params.clone()),
            params: self.clone(),
        }
    }
}

/// Represents:
///
/// c0 + c1 * u
///
/// where:
///
/// u^2 = params.non_residue
#[derive(Clone, Debug)]
pub struct Fp2<const N: usize = 4> {
    pub c0: Fp<N>,
    pub c1: Fp<N>,
    pub params: Fp2Params<N>,
}

impl<const N: usize> PartialEq for Fp2<N> {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0
            && self.c1 == other.c1
            && self.params.fp_params.modulus == other.params.fp_params.modulus
            && self.params.non_residue == other.params.non_residue
            && self.params.frobenius_coeffs_c1 == other.params.frobenius_coeffs_c1
    }
}

impl<const N: usize> Eq for Fp2<N> {}

impl<const N: usize> Fp2<N> {
    pub fn new(c0: Fp<N>, c1: Fp<N>, params: Fp2Params<N>) -> Self {
        assert_eq!(
            c0.params.modulus, params.fp_params.modulus,
            "c0 belongs to a different base field"
        );

        assert_eq!(
            c1.params.modulus, params.fp_params.modulus,
            "c1 belongs to a different base field"
        );

        Self { c0, c1, params }
    }

    pub fn zero(params: Fp2Params<N>) -> Self {
        params.zero()
    }

    pub fn one(params: Fp2Params<N>) -> Self {
        params.one()
    }

    pub fn from_u64(value: u64, params: Fp2Params<N>) -> Self {
        params.from_u64(value)
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.c0 == Fp::<N>::one(self.params.fp_params.clone()) && self.c1.is_zero()
    }

    fn check_field(&self, other: &Self) {
        assert_eq!(
            self.params.fp_params.modulus, other.params.fp_params.modulus,
            "Different base fields"
        );

        assert_eq!(
            self.params.non_residue, other.params.non_residue,
            "Different quadratic extensions"
        );
    }

    fn with_coefficients(&self, c0: Fp<N>, c1: Fp<N>) -> Self {
        Self {
            c0,
            c1,
            params: self.params.clone(),
        }
    }

    // ------------------------------------------------------------
    // Addition and subtraction
    // ------------------------------------------------------------

    pub fn add(&self, other: &Self) -> Self {
        self.check_field(other);

        self.with_coefficients(self.c0.add(&other.c0), self.c1.add(&other.c1))
    }

    pub fn sub(&self, other: &Self) -> Self {
        self.check_field(other);

        self.with_coefficients(self.c0.sub(&other.c0), self.c1.sub(&other.c1))
    }

    pub fn negate(&self) -> Self {
        self.with_coefficients(self.c0.negate(), self.c1.negate())
    }

    // ------------------------------------------------------------
    // Scalar multiplication
    // ------------------------------------------------------------

    /// Multiply both coefficients by an Fp element.
    pub fn mul_by_fp(&self, scalar: &Fp<N>) -> Self {
        assert_eq!(
            scalar.params.modulus, self.params.fp_params.modulus,
            "Scalar belongs to a different base field"
        );

        self.with_coefficients(self.c0.mul(scalar), self.c1.mul(scalar))
    }

    // ------------------------------------------------------------
    // Multiplication
    // ------------------------------------------------------------

    /// Karatsuba multiplication.
    ///
    /// (a + bu)(A + Bu)
    ///
    /// c0 = aA + non_residue * bB
    /// c1 = (a+b)(A+B) - aA - bB
    pub fn mul(&self, other: &Self) -> Self {
        self.check_field(other);

        let a_a = self.c0.mul(&other.c0);
        let b_b = self.c1.mul(&other.c1);

        let c0 = a_a.add(&self.params.non_residue.mul(&b_b));

        let a_plus_b = self.c0.add(&self.c1);
        let aa_plus_bb = other.c0.add(&other.c1);

        let c1 = a_plus_b.mul(&aa_plus_bb).sub(&a_a).sub(&b_b);

        self.with_coefficients(c0, c1)
    }

    /// Squaring using the Karatsuba formula.
    pub fn square(&self) -> Self {
        let a2 = self.c0.square();
        let b2 = self.c1.square();

        let c0 = a2.add(&self.params.non_residue.mul(&b2));

        let sum = self.c0.add(&self.c1);

        let c1 = sum.square().sub(&a2).sub(&b2);

        self.with_coefficients(c0, c1)
    }

    pub fn squared(&self) -> Self {
        self.square()
    }

    pub fn square_in_place(&mut self) {
        *self = self.square();
    }

    /// Alternative complex squaring formula.
    ///
    /// c0 = (a+b)(a+non_residue*b)
    ///      - ab - non_residue*ab
    ///
    /// c1 = 2ab
    pub fn square_complex(&self) -> Self {
        let ab = self.c0.mul(&self.c1);

        let a_plus_b = self.c0.add(&self.c1);

        let a_plus_nr_b = self.c0.add(&self.params.non_residue.mul(&self.c1));

        let c0 = a_plus_b
            .mul(&a_plus_nr_b)
            .sub(&ab)
            .sub(&self.params.non_residue.mul(&ab));

        let two = Fp::<N>::from_u64(2, self.params.fp_params.clone());

        let c1 = ab.mul(&two);

        self.with_coefficients(c0, c1)
    }

    // ------------------------------------------------------------
    // Exponentiation
    // ------------------------------------------------------------

    /// Exponentiation with a little-endian limb array.
    ///
    /// exponent[0] contains the least significant 64 bits.
    ///
    /// Supports arbitrary exponent lengths.
    pub fn pow(&self, exponent: &[u64]) -> Self {
        let mut result = self.params.one();
        let mut base = self.clone();

        for limb in exponent {
            let mut word = *limb;

            for _ in 0..64 {
                if (word & 1) != 0 {
                    result = result.mul(&base);
                }

                word >>= 1;
                base = base.square();
            }
        }

        result
    }

    // ------------------------------------------------------------
    // Inversion
    // ------------------------------------------------------------

    /// Multiplicative inverse.
    ///
    /// For x = a + bu:
    ///
    /// x^-1 = (a - bu) / (a^2 - non_residue*b^2)
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let a2 = self.c0.square();
        let b2 = self.c1.square();

        let denominator = a2.sub(&self.params.non_residue.mul(&b2));

        let denominator_inverse = denominator.inverse()?;

        let c0 = self.c0.mul(&denominator_inverse);

        let c1 = self.c1.mul(&denominator_inverse).negate();

        Some(self.with_coefficients(c0, c1))
    }

    /// In-place style inversion.
    ///
    /// Returns false if the element is zero.
    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(inverse) => {
                *self = inverse;
                true
            }
            None => false,
        }
    }

    // ------------------------------------------------------------
    // Frobenius
    // ------------------------------------------------------------

    /// Apply the Frobenius map.
    ///
    /// (c0 + c1*u)^p
    ///     = c0 + c1 * Frobenius_coeffs_c1[power % 2] * u
    ///
    /// This follows the coefficient-based quadratic-extension
    /// Frobenius mapping used by libff.
    pub fn frobenius_map(&self, power: usize) -> Self {
        let coeff = &self.params.frobenius_coeffs_c1[power % 2];

        self.with_coefficients(self.c0.clone(), self.c1.mul(coeff))
    }

    // ------------------------------------------------------------
    // Serialization
    // ------------------------------------------------------------

    /// Serialize as 2*N canonical little-endian u64 words.
    ///
    /// Layout:
    ///
    /// [c0.limbs..., c1.limbs...]
    ///
    /// Values are canonical, not Montgomery representations.
    pub fn to_words(&self) -> Vec<u64> {
        let mut words = Vec::with_capacity(2 * N);

        words.extend_from_slice(&self.c0.to_words());
        words.extend_from_slice(&self.c1.to_words());

        words
    }

    /// Deserialize 2*N canonical little-endian words.
    pub fn from_words(words: &[u64], params: Fp2Params<N>) -> Option<Self> {
        if words.len() != 2 * N {
            return None;
        }

        let c0 = Fp::<N>::from_words(&words[..N], params.fp_params.clone())?;

        let c1 = Fp::<N>::from_words(&words[N..], params.fp_params.clone())?;

        Some(Self::new(c0, c1, params))
    }

    // ------------------------------------------------------------
    // Square root
    // ------------------------------------------------------------

    /// Algebraic square root in Fp2.
    ///
    /// For:
    ///
    /// z = a + b*u
    /// norm(z) = a^2 - non_residue*b^2
    ///
    /// If z = (x + y*u)^2, then:
    ///
    /// x^2 = (a + sqrt(norm(z))) / 2
    /// y   = b / (2*x)
    ///
    /// The result is verified before returning.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(self.params.zero());
        }

        let fp_params = self.params.fp_params.clone();

        let two = Fp::<N>::from_u64(2, fp_params.clone());
        let inv_two = two.inverse()?;

        let nr = &self.params.non_residue;

        // Special case: b == 0.
        //
        // Try (sqrt(a), 0).
        // If that fails, try (0, sqrt(a / nr)).
        if self.c1.is_zero() {
            if let Some(x) = self.c0.sqrt() {
                let candidate = Self::new(x, Fp::<N>::zero(fp_params.clone()), self.params.clone());

                if candidate.square() == *self {
                    return Some(candidate);
                }
            }

            let nr_inv = nr.inverse()?;
            let y2 = self.c0.mul(&nr_inv);

            if let Some(y) = y2.sqrt() {
                let candidate = Self::new(Fp::<N>::zero(fp_params), y, self.params.clone());

                if candidate.square() == *self {
                    return Some(candidate);
                }
            }

            return None;
        }

        // norm = a^2 - nr*b^2
        let a2 = self.c0.mul(&self.c0);
        let b2 = self.c1.mul(&self.c1);
        let nr_b2 = nr.mul(&b2);

        let norm = a2.sub(&nr_b2);

        // A square in Fp2 must have a square norm in Fp.
        let alpha = norm.sqrt()?;

        // Try both square roots of norm:
        //
        // alpha and -alpha.
        for alpha_candidate in [alpha.clone(), alpha.negate()] {
            // x^2 = (a + alpha) / 2
            let x2 = self.c0.add(&alpha_candidate).mul(&inv_two);

            let Some(x) = x2.sqrt() else {
                continue;
            };

            if x.is_zero() {
                continue;
            }

            // y = b / (2*x)
            let denominator = two.mul(&x);

            let Some(denominator_inv) = denominator.inverse() else {
                continue;
            };

            let y = self.c1.mul(&denominator_inv);

            let candidate = Self::new(x, y, self.params.clone());

            // Mandatory verification.
            if candidate.square() == *self {
                return Some(candidate);
            }
        }

        None
    }
}

// ============================================================================
// BN254 Fq2 parameters
// ============================================================================

/// BN254 base field Fq.
///
/// BN254 uses four u64 limbs.
pub fn bn254_fq_params() -> FieldParams<4> {
    FieldParams::new([
        0x3c20_8c16_d87c_fd47,
        0x9781_6a91_6871_ca8d,
        0xb850_45b6_8181_585d,
        0x3064_4e72_e131_a029,
    ])
}

/// BN254 Fq2 parameters.
///
/// The extension is:
///
/// Fq2 = Fq[u] / (u^2 + 1)
///
/// Therefore:
///
/// non_residue = -1
///
/// Since the BN254 base-field modulus p is congruent to 3 mod 4:
///
/// u^p = -u
///
/// Hence the Frobenius coefficients are [1, -1].
pub fn bn254_fq2_params() -> Fp2Params<4> {
    let fp_params = bn254_fq_params();

    let one = Fp::<4>::one(fp_params.clone());
    let minus_one = one.negate();

    Fp2Params::new(fp_params, minus_one.clone(), [one, minus_one])
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> Fp2Params<4> {
        bn254_fq2_params()
    }

    fn element(a: u64, b: u64) -> Fp2<4> {
        let p = params();

        Fp2::new(
            Fp::<4>::from_u64(a, p.fp_params.clone()),
            Fp::<4>::from_u64(b, p.fp_params.clone()),
            p,
        )
    }

    #[test]
    fn test_zero_and_one() {
        let p = params();

        assert!(p.zero().is_zero());
        assert!(p.one().is_one());
    }

    #[test]
    fn test_add() {
        let a = element(3, 5);
        let b = element(7, 11);

        let result = a.add(&b);

        assert_eq!(result, element(10, 16));
    }

    #[test]
    fn test_sub() {
        let a = element(10, 16);
        let b = element(7, 11);

        assert_eq!(a.sub(&b), element(3, 5));
    }

    #[test]
    fn test_negate() {
        let a = element(3, 5);

        assert!(a.add(&a.negate()).is_zero());
    }

    #[test]
    fn test_mul() {
        // u^2 = -1
        //
        // (3 + 5u)(7 + 11u)
        // = 21 + 33u + 35u + 55u^2
        // = -34 + 68u

        let a = element(3, 5);
        let b = element(7, 11);

        let result = a.mul(&b);

        let p = params();

        let minus_34 = Fp::<4>::from_u64(34, p.fp_params.clone()).negate();

        let sixty_eight = Fp::<4>::from_u64(68, p.fp_params.clone());

        let expected = Fp2::new(minus_34, sixty_eight, p);

        assert_eq!(result, expected);
    }

    #[test]
    fn test_square() {
        let a = element(3, 5);

        assert_eq!(a.square(), a.mul(&a));
    }

    #[test]
    fn test_square_complex() {
        let a = element(3, 5);

        assert_eq!(a.square_complex(), a.square());
    }

    #[test]
    fn test_mul_by_fp() {
        let a = element(3, 5);

        let scalar = Fp::<4>::from_u64(7, a.params.fp_params.clone());

        assert_eq!(a.mul_by_fp(&scalar), element(21, 35));
    }

    #[test]
    fn test_inverse() {
        let a = element(3, 5);

        let inverse = a.inverse().unwrap();

        assert_eq!(a.mul(&inverse), a.params.one());
    }

    #[test]
    fn test_zero_inverse() {
        let zero = params().zero();

        assert!(zero.inverse().is_none());
    }

    #[test]
    fn test_pow() {
        let a = element(3, 5);

        let a2 = a.square();
        let a4 = a2.square();

        assert_eq!(a.pow(&[4]), a4);
    }

    #[test]
    fn test_frobenius_zero() {
        let a = element(3, 5);

        assert_eq!(a.frobenius_map(0), a);
    }

    #[test]
    fn test_frobenius_one() {
        // For BN254 Fq2 with u^2 = -1:
        //
        // (a + bu)^p = a - bu

        let a = element(3, 5);

        let expected = a.with_coefficients(a.c0.clone(), a.c1.negate());

        assert_eq!(a.frobenius_map(1), expected);
    }

    #[test]
    fn test_frobenius_two() {
        let a = element(3, 5);

        assert_eq!(a.frobenius_map(2), a);
    }

    #[test]
    fn test_serialization() {
        let a = element(123, 456);

        let words = a.to_words();

        assert_eq!(words.len(), 8);

        let b = Fp2::from_words(&words, params()).unwrap();

        assert_eq!(a, b);
    }

    #[test]
    fn test_invalid_serialization_length() {
        let words = [1u64, 2, 3, 4];

        assert!(Fp2::from_words(&words, params()).is_none());
    }

    #[test]
    fn test_sqrt_square() {
        let a = element(3, 5);

        let square = a.square();

        let root = square.sqrt().unwrap();

        assert_eq!(root.square(), square);
    }

    #[test]
    fn test_sqrt_zero() {
        let zero = params().zero();

        assert_eq!(zero.sqrt().unwrap(), zero);
    }

    #[test]
    fn test_sqrt_non_residue() {
        // In Fq2, (0, 1) is not necessarily a non-residue.
        //
        // This test checks that sqrt either returns a valid root
        // or correctly reports that no root exists.
        let a = element(3, 5);

        let result = a.sqrt();

        if let Some(root) = result {
            assert_eq!(root.square(), a);
        }
    }

    #[test]
    fn test_mul_inverse_identity() {
        let a = element(17, 29);

        let inverse = a.inverse().unwrap();

        assert_eq!(a.mul(&inverse), params().one());
    }
}
