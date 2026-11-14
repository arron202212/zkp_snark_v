use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

/// Parameters for the cubic extension:
/// ```text
///     Fp3 = Fp[u] / (u^3 - non_residue)
/// ```
/// An element is represented as:
///```text
///     c0 + c1*u + c2*u^2
/// ```
#[derive(Clone, Debug)]
pub struct Fp3Params {
    pub fp_params: FieldParams,

    /// u^3 = non_residue.
    pub non_residue: Fp,

    /// Frobenius coefficients for c1.
    pub frobenius_coeffs_c1: [Fp; 3],

    /// Frobenius coefficients for c2.
    pub frobenius_coeffs_c2: [Fp; 3],
}

impl Fp3Params {
    pub fn new(
        fp_params: FieldParams,
        non_residue: Fp,
        frobenius_coeffs_c1: [Fp; 3],
        frobenius_coeffs_c2: [Fp; 3],
    ) -> Self {
        assert_eq!(
            fp_params.modulus, non_residue.params.modulus,
            "non_residue belongs to a different base field"
        );

        for coeff in &frobenius_coeffs_c1 {
            assert_eq!(
                fp_params.modulus, coeff.params.modulus,
                "Frobenius c1 coefficient belongs to a different base field"
            );
        }

        for coeff in &frobenius_coeffs_c2 {
            assert_eq!(
                fp_params.modulus, coeff.params.modulus,
                "Frobenius c2 coefficient belongs to a different base field"
            );
        }

        Self {
            fp_params,
            non_residue,
            frobenius_coeffs_c1,
            frobenius_coeffs_c2,
        }
    }

    pub fn zero(&self) -> Fp3 {
        Fp3 {
            c0: Fp::zero(self.fp_params.clone()),
            c1: Fp::zero(self.fp_params.clone()),
            c2: Fp::zero(self.fp_params.clone()),
            params: self.clone(),
        }
    }

    pub fn one(&self) -> Fp3 {
        Fp3 {
            c0: Fp::one(self.fp_params.clone()),
            c1: Fp::zero(self.fp_params.clone()),
            c2: Fp::zero(self.fp_params.clone()),
            params: self.clone(),
        }
    }

    /// Embed an Fp element into Fp3.
    pub fn from_u64(&self, value: u64) -> Fp3 {
        Fp3 {
            c0: Fp::from_u64(value, self.fp_params.clone()),
            c1: Fp::zero(self.fp_params.clone()),
            c2: Fp::zero(self.fp_params.clone()),
            params: self.clone(),
        }
    }
}

/// Element of the cubic extension field Fp3.
/// ```text
///     c0 + c1*u + c2*u^2
/// ```
/// with:
/// ```text
///     u^3 = non_residue
/// ```
#[derive(Clone, Debug)]
pub struct Fp3 {
    pub c0: Fp,
    pub c1: Fp,
    pub c2: Fp,
    pub params: Fp3Params,
}

impl PartialEq for Fp3 {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0
            && self.c1 == other.c1
            && self.c2 == other.c2
            && self.params.fp_params.modulus == other.params.fp_params.modulus
            && self.params.non_residue == other.params.non_residue
            && self.params.frobenius_coeffs_c1 == other.params.frobenius_coeffs_c1
            && self.params.frobenius_coeffs_c2 == other.params.frobenius_coeffs_c2
    }
}

impl Eq for Fp3 {}

impl Fp3 {
    pub fn new(c0: Fp, c1: Fp, c2: Fp, params: Fp3Params) -> Self {
        assert_eq!(
            c0.params.modulus, params.fp_params.modulus,
            "c0 belongs to a different base field"
        );

        assert_eq!(
            c1.params.modulus, params.fp_params.modulus,
            "c1 belongs to a different base field"
        );

        assert_eq!(
            c2.params.modulus, params.fp_params.modulus,
            "c2 belongs to a different base field"
        );

        Self { c0, c1, c2, params }
    }

    pub fn zero(params: Fp3Params) -> Self {
        params.zero()
    }

    pub fn one(params: Fp3Params) -> Self {
        params.one()
    }

    pub fn from_u64(value: u64, params: Fp3Params) -> Self {
        params.from_u64(value)
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.c0 == Fp::one(self.params.fp_params.clone()) && self.c1.is_zero() && self.c2.is_zero()
    }

    fn check_field(&self, other: &Self) {
        assert_eq!(
            self.params.fp_params.modulus, other.params.fp_params.modulus,
            "Different base fields"
        );

        assert_eq!(
            self.params.non_residue, other.params.non_residue,
            "Different cubic extensions"
        );
    }

    fn with_coefficients(&self, c0: Fp, c1: Fp, c2: Fp) -> Self {
        Self {
            c0,
            c1,
            c2,
            params: self.params.clone(),
        }
    }

    // ------------------------------------------------------------
    // Addition and subtraction
    // ------------------------------------------------------------

    /// Component-wise addition.
    pub fn add(&self, other: &Self) -> Self {
        self.check_field(other);

        self.with_coefficients(
            self.c0.add(&other.c0),
            self.c1.add(&other.c1),
            self.c2.add(&other.c2),
        )
    }

    /// Component-wise subtraction.
    pub fn sub(&self, other: &Self) -> Self {
        self.check_field(other);

        self.with_coefficients(
            self.c0.sub(&other.c0),
            self.c1.sub(&other.c1),
            self.c2.sub(&other.c2),
        )
    }

    /// Component-wise negation.
    pub fn negate(&self) -> Self {
        self.with_coefficients(self.c0.negate(), self.c1.negate(), self.c2.negate())
    }

    // ------------------------------------------------------------
    // Scalar multiplication
    // ------------------------------------------------------------

    /// Multiply all three coefficients by an Fp element.
    pub fn mul_by_fp(&self, scalar: &Fp) -> Self {
        assert_eq!(
            scalar.params.modulus, self.params.fp_params.modulus,
            "Scalar belongs to a different base field"
        );

        self.with_coefficients(
            self.c0.mul(scalar),
            self.c1.mul(scalar),
            self.c2.mul(scalar),
        )
    }

    // ------------------------------------------------------------
    // Multiplication
    // ------------------------------------------------------------

    /// Karatsuba multiplication.
    ///
    /// This follows the formula in libff fp3.tcc:
    ///
    /// aA = a * A
    /// bB = b * B
    /// cC = c * C
    ///
    /// result.c0 =
    ///     aA + nr * ((b+c)(B+C) - bB - cC)
    ///
    /// result.c1 =
    ///     (a+b)(A+B) - aA - bB + nr*cC
    ///
    /// result.c2 =
    ///     (a+c)(A+C) - aA + bB - cC
    pub fn mul(&self, other: &Self) -> Self {
        self.check_field(other);

        let a = &self.c0;
        let b = &self.c1;
        let c = &self.c2;

        let aa = &other.c0;
        let bb = &other.c1;
        let cc = &other.c2;

        let a_a = a.mul(aa);
        let b_b = b.mul(bb);
        let c_c = c.mul(cc);

        let b_plus_c = b.add(c);
        let bb_plus_cc = bb.add(cc);

        let c0 = a_a.add(
            &self
                .params
                .non_residue
                .mul(&b_plus_c.mul(&bb_plus_cc).sub(&b_b).sub(&c_c)),
        );

        let a_plus_b = a.add(b);
        let aa_plus_bb = aa.add(bb);

        let c1 = a_plus_b
            .mul(&aa_plus_bb)
            .sub(&a_a)
            .sub(&b_b)
            .add(&self.params.non_residue.mul(&c_c));

        let a_plus_c = a.add(c);
        let aa_plus_cc = aa.add(cc);

        let c2 = a_plus_c.mul(&aa_plus_cc).sub(&a_a).add(&b_b).sub(&c_c);

        self.with_coefficients(c0, c1, c2)
    }

    // ------------------------------------------------------------
    // Squaring
    // ------------------------------------------------------------

    /// CH-SQR2 squaring formula from libff fp3.tcc.
    ///
    /// s0 = a^2
    /// s1 = 2ab
    /// s2 = (a-b+c)^2
    /// s3 = 2bc
    /// s4 = c^2
    ///
    /// result:
    ///
    /// c0 = s0 + nr*s3
    /// c1 = s1 + nr*s4
    /// c2 = s1 + s2 + s3 - s0 - s4
    pub fn square(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;
        let c = &self.c2;

        let s0 = a.square();

        let ab = a.mul(b);
        let s1 = ab.add(&ab);

        let s2 = a.sub(b).add(c).square();

        let bc = b.mul(c);
        let s3 = bc.add(&bc);

        let s4 = c.square();

        let c0 = s0.add(&self.params.non_residue.mul(&s3));

        let c1 = s1.add(&self.params.non_residue.mul(&s4));

        let c2 = s1.add(&s2).add(&s3).sub(&s0).sub(&s4);

        self.with_coefficients(c0, c1, c2)
    }

    /// Alias matching the libff squared() method.
    pub fn squared(&self) -> Self {
        self.square()
    }

    // ------------------------------------------------------------
    // Exponentiation
    // ------------------------------------------------------------

    /// Exponentiation using little-endian limbs.
    ///
    /// Supports arbitrary exponent lengths.
    ///
    /// For Fp3 Tonelli-Shanks, the exponent p^3 - 1
    /// can require up to 768 bits when p is 256 bits.
    pub fn pow(&self, exponent: &[u64]) -> Self {
        let mut result = self.params.one();
        let mut base = self.clone();

        for limb in exponent {
            let mut word = *limb;

            for _ in 0..64 {
                if word & 1 == 1 {
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
    /// This follows Algorithm 17 from the libff implementation.
    ///
    /// t0 = a^2
    /// t1 = b^2
    /// t2 = c^2
    /// t3 = ab
    /// t4 = ac
    /// t5 = bc
    ///
    /// c0 = t0 - nr*t5
    /// c1 = nr*t2 - t3
    /// c2 = t1 - t4
    ///
    /// t6 = (a*c0 + nr*(c*c1 + b*c2))^-1
    ///
    /// result = (t6*c0, t6*c1, t6*c2)
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let a = &self.c0;
        let b = &self.c1;
        let c = &self.c2;

        let t0 = a.square();
        let t1 = b.square();
        let t2 = c.square();

        let t3 = a.mul(b);
        let t4 = a.mul(c);
        let t5 = b.mul(c);

        let nr = &self.params.non_residue;

        let c0 = t0.sub(&nr.mul(&t5));

        let c1 = nr.mul(&t2).sub(&t3);

        let c2 = t1.sub(&t4);

        let denominator = a.mul(&c0).add(&nr.mul(&c.mul(&c1).add(&b.mul(&c2))));

        let denominator_inverse = denominator.inverse()?;

        let result = self.with_coefficients(
            denominator_inverse.mul(&c0),
            denominator_inverse.mul(&c1),
            denominator_inverse.mul(&c2),
        );

        // Defensive verification of the inverse formula.
        debug_assert_eq!(self.mul(&result), self.params.one());

        Some(result)
    }

    /// In-place style inversion.
    ///
    /// Returns false if the element is zero.
    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(value) => {
                *self = value;
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
    /// (c0 + c1*u + c2*u^2)^(p^power)
    ///
    /// = c0
    ///   + Frobenius_coeffs_c1[power % 3] * c1 * u
    ///   + Frobenius_coeffs_c2[power % 3] * c2 * u^2
    pub fn frobenius_map(&self, power: usize) -> Self {
        let index = power % 3;

        self.with_coefficients(
            self.c0.clone(),
            self.params.frobenius_coeffs_c1[index].mul(&self.c1),
            self.params.frobenius_coeffs_c2[index].mul(&self.c2),
        )
    }

    // ------------------------------------------------------------
    // Serialization
    // ------------------------------------------------------------

    /// Serialize in coefficient order:
    ///
    /// c0 words || c1 words || c2 words
    ///
    /// Each Fp coefficient contributes four canonical u64 words.
    pub fn to_words(&self) -> Vec<u64> {
        let mut words = Vec::with_capacity(12);

        words.extend_from_slice(&self.c0.to_words());
        words.extend_from_slice(&self.c1.to_words());
        words.extend_from_slice(&self.c2.to_words());

        words
    }

    /// Deserialize twelve canonical little-endian words.
    pub fn from_words(words: &[u64], params: Fp3Params) -> Option<Self> {
        if words.len() != 12 {
            return None;
        }

        let c0 = Fp::from_words(&words[0..4], params.fp_params.clone())?;

        let c1 = Fp::from_words(&words[4..8], params.fp_params.clone())?;

        let c2 = Fp::from_words(&words[8..12], params.fp_params.clone())?;

        Some(Self::new(c0, c1, c2, params))
    }

    // ------------------------------------------------------------
    // Square root
    // ------------------------------------------------------------

    /// Tonelli-Shanks square root in Fp3.
    ///
    /// The multiplicative group of Fp3 has order:
    /// ```text
    ///     p^3 - 1
    /// ```
    /// This implementation uses up to twelve 64-bit limbs.
    ///
    /// Returns None if the element is not a quadratic residue.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(self.params.zero());
        }

        let p = self.params.fp_params.modulus;

        // Compute p^3 - 1.
        let p_squared = mul_words(&p, &p);
        let mut order = mul_words(&p_squared, &p);

        sub_one_words(&mut order);

        // Factor:
        //
        // p^3 - 1 = q * 2^s
        //
        // where q is odd.
        let mut q = order.clone();
        let mut s = 0usize;

        while q[0] & 1 == 0 {
            shr_one_words(&mut q);
            s += 1;
        }

        if s == 0 {
            return None;
        }

        let one = self.params.one();

        // Euler criterion:
        //
        // self^((p^3 - 1)/2) = 1
        //
        // for a nonzero quadratic residue.
        let legendre_exp = shr_words(&order);

        let legendre = self.pow(&legendre_exp);

        if legendre != one {
            return None;
        }

        // Find a quadratic non-residue z in Fp3.
        let z = self.find_non_residue(&legendre_exp)?;

        // x = self^((q+1)/2)
        let mut q_plus_one = q.clone();
        add_one_words(&mut q_plus_one);

        let x_exp = shr_words(&q_plus_one);

        let mut x = self.pow(&x_exp);

        // t = self^q
        let mut t = self.pow(&q);

        // c = z^q
        let mut c = z.pow(&q);

        let mut m = s;

        while !t.is_one() {
            let mut i = 0usize;
            let mut t2 = t.square();

            while i < m && !t2.is_one() {
                t2 = t2.square();
                i += 1;
            }

            if i == m {
                return None;
            }

            // b = c^(2^(m-i-1))
            let mut b = c.clone();

            for _ in 0..(m - i - 1) {
                b = b.square();
            }

            let b2 = b.square();

            x = x.mul(&b);
            t = t.mul(&b2);
            c = b2;

            m = i;
        }

        // Verify the candidate root.
        if x.square() == *self { Some(x) } else { None }
    }

    fn find_non_residue(&self, legendre_exp: &[u64]) -> Option<Self> {
        let one = self.params.one();

        // Try candidates:
        //
        //     z = k + u
        //
        // with k = 0, 1, 2, ...
        //
        // A candidate is a quadratic non-residue if:
        //
        //     z^((p^3 - 1)/2) != 1
        //
        // This search assumes the supplied cubic extension
        // parameters define a valid field.
        for k in 0u64..4096 {
            let candidate = Self::new(
                Fp::from_u64(k, self.params.fp_params.clone()),
                Fp::one(self.params.fp_params.clone()),
                Fp::zero(self.params.fp_params.clone()),
                self.params.clone(),
            );

            if candidate.pow(legendre_exp) != one {
                return Some(candidate);
            }
        }

        None
    }
}

// ============================================================================
// Wide integer helpers
// ============================================================================

/// Multiply arbitrary little-endian limb arrays.
///
/// The output length is a.len() + b.len().
///
/// Each inner multiply-accumulate fits into u128:
///
/// (2^64 - 1)^2 + (2^64 - 1) + (2^64 - 1) = 2^128 - 1.
fn mul_words(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; a.len() + b.len()];

    for i in 0..a.len() {
        let mut carry = 0u128;

        for j in 0..b.len() {
            let index = i + j;

            let uv = (a[i] as u128) * (b[j] as u128) + (out[index] as u128) + carry;

            out[index] = uv as u64;
            carry = uv >> 64;
        }

        let mut index = i + b.len();
        let mut c = carry as u64;

        while c != 0 && index < out.len() {
            let (value, overflow) = out[index].overflowing_add(c);

            out[index] = value;
            c = overflow as u64;
            index += 1;
        }

        debug_assert_eq!(c, 0);
    }

    out
}

fn sub_one_words(a: &mut [u64]) {
    let mut borrow = true;

    for limb in a.iter_mut() {
        if !borrow {
            break;
        }

        let (value, underflow) = limb.overflowing_sub(1);

        *limb = value;
        borrow = underflow;
    }

    debug_assert!(!borrow);
}

fn add_one_words(a: &mut [u64]) {
    let mut carry = true;

    for limb in a.iter_mut() {
        if !carry {
            break;
        }

        let (value, overflow) = limb.overflowing_add(1);

        *limb = value;
        carry = overflow;
    }

    debug_assert!(!carry);
}

fn shr_one_words(a: &mut [u64]) {
    let mut carry = 0u64;

    for i in (0..a.len()).rev() {
        let next_carry = a[i] << 63;

        a[i] = (a[i] >> 1) | carry;

        carry = next_carry;
    }
}

fn shr_words(a: &[u64]) -> Vec<u64> {
    let mut result = a.to_vec();
    shr_one_words(&mut result);
    result
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // A small prime field Fp with p = 7.
    //
    // x^3 - 2 is irreducible over F7:
    // 2 is not a cube in F7.
    //
    // Therefore F7[u]/(u^3 - 2) is a field with 7^3 elements.
    fn test_params() -> Fp3Params {
        let fp_params = FieldParams::new([7, 0, 0, 0]);

        let one = Fp::one(fp_params.clone());

        let _minus_one = one.negate();

        let nr = Fp::from_u64(2, fp_params.clone());

        // For p = 7 and u^3 = 2:
        //
        // u^7 = 2*u
        // u^14 = 4*u^2
        //
        // Frobenius:
        //
        // (a + bu + cu^2)^7
        // = a + 2*b*u + 4*c*u^2
        //
        // Its square:
        //
        // (a + bu + cu^2)^49
        // = a + 4*b*u + 2*c*u^2
        //
        // The third power is identity.
        Fp3Params::new(
            fp_params,
            nr,
            [
                one.clone(),
                Fp::from_u64(2, FieldParams::new([7, 0, 0, 0])),
                Fp::from_u64(4, FieldParams::new([7, 0, 0, 0])),
            ],
            [
                one,
                Fp::from_u64(4, FieldParams::new([7, 0, 0, 0])),
                Fp::from_u64(2, FieldParams::new([7, 0, 0, 0])),
            ],
        )
    }

    fn element(a: u64, b: u64, c: u64) -> Fp3 {
        let p = test_params();

        Fp3::new(
            Fp::from_u64(a, p.fp_params.clone()),
            Fp::from_u64(b, p.fp_params.clone()),
            Fp::from_u64(c, p.fp_params.clone()),
            p,
        )
    }

    #[test]
    fn test_zero_and_one() {
        let p = test_params();

        assert!(p.zero().is_zero());
        assert!(p.one().is_one());
    }

    #[test]
    fn test_add() {
        let a = element(1, 2, 3);
        let b = element(3, 4, 5);

        assert_eq!(a.add(&b), element(4, 6, 1),);
    }

    #[test]
    fn test_sub() {
        let a = element(4, 6, 1);
        let b = element(3, 4, 5);

        assert_eq!(a.sub(&b), element(1, 2, 3),);
    }

    #[test]
    fn test_negate() {
        let a = element(1, 2, 3);

        assert!(a.add(&a.negate()).is_zero());
    }

    #[test]
    fn test_mul() {
        // u^3 = 2 over F7.
        //
        // (1 + 2u + 3u²)(3 + 4u + 5u²)
        //
        // c0 = 1*3 + 2*(2*5 + 3*4)
        //    = 47 mod 7 = 5
        //
        // c1 = 1*4 + 2*3 + 2*(3*5)
        //    = 40 mod 7 = 5
        //
        // c2 = 1*5 + 2*4 + 3*3
        //    = 22 mod 7 = 1

        let a = element(1, 2, 3);
        let b = element(3, 4, 5);

        assert_eq!(a.mul(&b), element(5, 5, 1),);
    }

    #[test]
    fn test_square_matches_mul() {
        let a = element(1, 2, 3);

        assert_eq!(a.square(), a.mul(&a),);
    }

    #[test]
    fn test_squared_alias() {
        let a = element(1, 2, 3);

        assert_eq!(a.squared(), a.square(),);
    }

    #[test]
    fn test_scalar_mul() {
        let a = element(1, 2, 3);

        let scalar = Fp::from_u64(3, a.params.fp_params.clone());

        assert_eq!(a.mul_by_fp(&scalar), element(3, 6, 2),);
    }

    #[test]
    fn test_inverse() {
        let a = element(1, 2, 3);

        let inverse = a.inverse().unwrap();

        assert_eq!(a.mul(&inverse), test_params().one(),);
    }

    #[test]
    fn test_zero_inverse() {
        let zero = test_params().zero();

        assert!(zero.inverse().is_none());
    }

    #[test]
    fn test_invert() {
        let mut a = element(1, 2, 3);

        assert!(a.invert());

        assert_eq!(a.mul(&element(1, 2, 3)), test_params().one(),);
    }

    #[test]
    fn test_pow() {
        let a = element(1, 2, 3);

        assert_eq!(a.pow(&[2]), a.square(),);
    }

    #[test]
    fn test_frobenius_zero() {
        let a = element(1, 2, 3);

        assert_eq!(a.frobenius_map(0), a,);
    }

    #[test]
    fn test_frobenius_one() {
        let a = element(1, 2, 3);

        assert_eq!(a.frobenius_map(1), element(1, 4, 5),);
    }

    #[test]
    fn test_frobenius_two() {
        let a = element(1, 2, 3);

        assert_eq!(a.frobenius_map(2), element(1, 1, 6),);
    }

    #[test]
    fn test_frobenius_three_is_identity() {
        let a = element(1, 2, 3);

        assert_eq!(a.frobenius_map(3), a,);
    }

    #[test]
    fn test_serialization() {
        let a = element(1, 2, 3);

        let words = a.to_words();

        assert_eq!(words.len(), 12);

        let b = Fp3::from_words(&words, test_params()).unwrap();

        assert_eq!(a, b);
    }

    #[test]
    fn test_invalid_serialization_length() {
        let words = [1u64, 2, 3, 4];

        assert!(Fp3::from_words(&words, test_params(),).is_none());
    }

    #[test]
    fn test_sqrt_zero() {
        let zero = test_params().zero();

        assert_eq!(zero.sqrt().unwrap(), zero,);
    }

    #[test]
    fn test_sqrt_square() {
        let a = element(1, 2, 3);

        let square = a.square();

        let root = square.sqrt().unwrap();

        assert_eq!(root.square(), square,);
    }

    #[test]
    fn test_wide_multiplication() {
        let a = [u64::MAX, u64::MAX, 0, 0];
        let b = [2u64, 0, 0, 0];

        let result = mul_words(&a, &b);

        assert_eq!(result.len(), 8);
        assert_eq!(result[0], u64::MAX - 1);
        assert_eq!(result[1], u64::MAX);
        assert_eq!(result[2], 1);
    }
}
