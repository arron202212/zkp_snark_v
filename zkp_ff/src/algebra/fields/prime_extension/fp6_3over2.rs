use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};

/// Fp6 = Fp2[w] / (w^3 - non_residue).
///
/// Frobenius coefficients:
///   c0' = Frobenius(c0)
///   c1' = coeff_c1[k] * Frobenius(c1)
///   c2' = coeff_c2[k] * Frobenius(c2)
#[derive(Clone)]
pub struct Fp6Params {
    pub fp2_params: Fp2Params,
    pub non_residue: Fp2,
    pub frobenius_coeffs_c1: [Fp2; 6],
    pub frobenius_coeffs_c2: [Fp2; 6],
}

#[derive(Clone)]
pub struct Fp6 {
    pub c0: Fp2,
    pub c1: Fp2,
    pub c2: Fp2,
    pub params: Fp6Params,
}

impl Fp6 {
    pub fn new(c0: Fp2, c1: Fp2, c2: Fp2, params: Fp6Params) -> Self {
        Self { c0, c1, c2, params }
    }

    pub fn zero(params: Fp6Params) -> Self {
        let z = Fp2::zero(params.fp2_params.clone());

        Self::new(z.clone(), z.clone(), z, params)
    }

    pub fn one(params: Fp6Params) -> Self {
        let one = Fp2::one(params.fp2_params.clone());
        let zero = Fp2::zero(params.fp2_params.clone());

        Self::new(one, zero.clone(), zero, params)
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero() && self.c2.is_zero()
    }

    pub fn equals(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }

    /// Multiply an Fp2 element by w^3 = non_residue.
    pub fn mul_by_non_residue(&self, x: &Fp2) -> Fp2 {
        self.params.non_residue.mul(x)
    }

    pub fn add(&self, other: &Self) -> Self {
        Self::new(
            self.c0.add(&other.c0),
            self.c1.add(&other.c1),
            self.c2.add(&other.c2),
            self.params.clone(),
        )
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self::new(
            self.c0.sub(&other.c0),
            self.c1.sub(&other.c1),
            self.c2.sub(&other.c2),
            self.params.clone(),
        )
    }

    pub fn negate(&self) -> Self {
        Self::new(
            self.c0.negate(),
            self.c1.negate(),
            self.c2.negate(),
            self.params.clone(),
        )
    }

    pub fn mul_by_fp(&self, scalar: &Fp) -> Self {
        Self::new(
            self.c0.mul_by_fp(scalar),
            self.c1.mul_by_fp(scalar),
            self.c2.mul_by_fp(scalar),
            self.params.clone(),
        )
    }

    /// Karatsuba multiplication.
    ///
    /// Fp6 = Fp2[w] / (w^3 - xi).
    pub fn mul(&self, other: &Self) -> Self {
        let a = &self.c0;
        let b = &self.c1;
        let c = &self.c2;

        let aa = &other.c0;
        let bb = &other.c1;
        let cc = &other.c2;

        let a_aa = a.mul(aa);
        let b_bb = b.mul(bb);
        let c_cc = c.mul(cc);

        let nr_term = self.mul_by_non_residue(&b.add(c).mul(&bb.add(cc)).sub(&b_bb).sub(&c_cc));

        let r0 = a_aa.add(&nr_term);

        let r1 = a
            .add(b)
            .mul(&aa.add(bb))
            .sub(&a_aa)
            .sub(&b_bb)
            .add(&self.mul_by_non_residue(&c_cc));

        let r2 = a.add(c).mul(&aa.add(cc)).sub(&a_aa).add(&b_bb).sub(&c_cc);

        Self::new(r0, r1, r2, self.params.clone())
    }

    /// CH-SQR2 squaring formula from libff.
    pub fn squared(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;
        let c = &self.c2;

        let s0 = a.squared();

        let ab = a.mul(b);
        let s1 = ab.add(&ab);

        let s2 = a.sub(b).add(c).squared();

        let bc = b.mul(c);
        let s3 = bc.add(&bc);

        let s4 = c.squared();

        let r0 = s0.add(&self.mul_by_non_residue(&s3));

        let r1 = s1.add(&self.mul_by_non_residue(&s4));

        let r2 = s1.add(&s2).add(&s3).sub(&s0).sub(&s4);

        Self::new(r0, r1, r2, self.params.clone())
    }

    pub fn square(&mut self) {
        *self = self.squared();
    }

    /// Algorithm 17:
    /// High-Speed Software Implementation of the Optimal Ate Pairing
    /// over Barreto-Naehrig Curves.
    pub fn inverse(&self) -> Option<Self> {
        let a = &self.c0;
        let b = &self.c1;
        let c = &self.c2;

        let t0 = a.squared();
        let t1 = b.squared();
        let t2 = c.squared();

        let t3 = a.mul(b);
        let t4 = a.mul(c);
        let t5 = b.mul(c);

        let c0 = t0.sub(&self.mul_by_non_residue(&t5));

        let c1 = self.mul_by_non_residue(&t2).sub(&t3);

        // The minus sign is intentional.
        // This is the corrected expression from libff.
        let c2 = t1.sub(&t4);

        let denominator = a
            .mul(&c0)
            .add(&self.mul_by_non_residue(&c.mul(&c1).add(&b.mul(&c2))));

        let denominator_inv = denominator.inverse()?;

        Some(Self::new(
            denominator_inv.mul(&c0),
            denominator_inv.mul(&c1),
            denominator_inv.mul(&c2),
            self.params.clone(),
        ))
    }

    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(inv) => {
                *self = inv;
                true
            }
            None => false,
        }
    }

    /// Frobenius map.
    ///
    /// The coefficient arrays must be initialized for the
    /// exact Fp2 extension and modulus being used.
    pub fn frobenius_map(&self, power: usize) -> Self {
        let k = power % 6;

        Self::new(
            self.c0.frobenius_map(power),
            self.params.frobenius_coeffs_c1[k].mul(&self.c1.frobenius_map(power)),
            self.params.frobenius_coeffs_c2[k].mul(&self.c2.frobenius_map(power)),
            self.params.clone(),
        )
    }

    /// Binary exponentiation.
    ///
    /// exponent is little-endian u64 limbs.
    pub fn pow(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for limb in exponent {
            let mut e = *limb;

            for _ in 0..64 {
                if (e & 1) != 0 {
                    result = result.mul(&base);
                }

                e >>= 1;
                base = base.squared();
            }
        }

        result
    }

    pub fn to_words(&self) -> Vec<u64> {
        let mut words = self.c0.to_words();

        words.extend(self.c1.to_words());
        words.extend(self.c2.to_words());

        words
    }

    /// Each Fp2 element occupies 8 u64 words.
    pub fn from_words(words: &[u64], params: Fp6Params) -> Option<Self> {
        if words.len() != 24 {
            return None;
        }

        let c0 = Fp2::from_words(&words[0..8], params.fp2_params.clone())?;

        let c1 = Fp2::from_words(&words[8..16], params.fp2_params.clone())?;

        let c2 = Fp2::from_words(&words[16..24], params.fp2_params.clone())?;

        Some(Self::new(c0, c1, c2, params))
    }

    /// Tonelli-Shanks square root in Fp6.
    ///
    /// Returns None if the element is not a square.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero(self.params.clone()));
        }

        let modulus = self.params.fp2_params.fp_params.modulus;

        // N = p^6.
        let mut n = [0u64; 24];
        n[0] = 1;

        for _ in 0..6 {
            n = mul_24_by_4(&n, &modulus);
        }

        // N - 1 = q * 2^s.
        let mut q = n;
        sub_one(&mut q);

        let mut s = 0usize;
        while !is_odd(&q) {
            shr1(&mut q);
            s += 1;
        }

        // x = a^((q+1)/2)
        let mut q_plus_one = q;
        add_one(&mut q_plus_one);
        shr1(&mut q_plus_one);

        let mut x = self.pow(&q_plus_one);
        let mut t = self.pow(&q);

        if t.is_one() {
            return Some(x);
        }

        // Find a quadratic non-residue in Fp2, then embed it into Fp6.
        //
        // Since Fp6/Fp2 has odd degree 3, a nonsquare in Fp2
        // remains a nonsquare in Fp6.
        //
        // Euler criterion:
        //     z^((p^6 - 1) / 2) == -1
        let mut legendre_exp = n;
        sub_one(&mut legendre_exp);
        shr1(&mut legendre_exp);

        let mut z: Option<Fp6> = None;

        for candidate in 0u64..4096 {
            let fp_candidate = Fp::from_u64(candidate, self.params.fp2_params.fp_params.clone());

            let fp_one = Fp::one(self.params.fp2_params.fp_params.clone());

            // Candidate = candidate + u.
            let fp2_candidate = Fp2::new(fp_candidate, fp_one, self.params.fp2_params.clone());

            let candidate_z = Self::new(
                fp2_candidate,
                Fp2::zero(self.params.fp2_params.clone()),
                Fp2::zero(self.params.fp2_params.clone()),
                self.params.clone(),
            );

            let test = candidate_z.pow(&legendre_exp);

            if test.equals(&Self::one(self.params.clone()).negate()) {
                z = Some(candidate_z);
                break;
            }
        }

        let z = z?;

        let mut c = z.pow(&q);
        let mut m = s;

        while !t.is_one() {
            let mut i = 1usize;
            let mut t2 = t.squared();

            while i < m && !t2.is_one() {
                t2 = t2.squared();
                i += 1;
            }

            if i == m {
                return None;
            }

            let mut exponent = [0u64; 24];
            exponent[0] = 1u64 << (m - i - 1);

            let b = c.pow(&exponent);
            let b2 = b.squared();

            x = x.mul(&b);
            t = t.mul(&b2);
            c = b2;
            m = i;
        }

        // Defensive verification.
        if x.squared().equals(self) {
            Some(x)
        } else {
            None
        }
    }
}

// -----------------------------------------------------------------------------
// 24-limb integer helpers used by Tonelli-Shanks.
// Limbs are little-endian.
// -----------------------------------------------------------------------------

fn is_odd(a: &[u64; 24]) -> bool {
    (a[0] & 1) == 1
}

fn is_zero_words(a: &[u64; 24]) -> bool {
    a.iter().all(|x| *x == 0)
}

fn shr1(a: &mut [u64; 24]) {
    let mut carry = 0u64;

    for i in (0..24).rev() {
        let next = a[i] & 1;
        a[i] = (a[i] >> 1) | (carry << 63);
        carry = next;
    }
}

fn add_one(a: &mut [u64; 24]) {
    for limb in a.iter_mut() {
        let (value, carry) = limb.overflowing_add(1);
        *limb = value;

        if !carry {
            break;
        }
    }
}

fn sub_one(a: &mut [u64; 24]) {
    for limb in a.iter_mut() {
        let (value, borrow) = limb.overflowing_sub(1);
        *limb = value;

        if !borrow {
            break;
        }
    }
}

/// Multiply a 24-limb integer by a 4-limb integer.
///
/// The result fits in 24 limbs for the p^6 construction
/// used by this implementation.
fn mul_24_by_4(a: &[u64; 24], b: &[u64; 4]) -> [u64; 24] {
    let mut out = [0u64; 24];

    for i in 0..24 {
        let mut carry = 0u128;

        for j in 0..4 {
            let k = i + j;

            if k >= 24 {
                break;
            }

            let value = (a[i] as u128) * (b[j] as u128) + (out[k] as u128) + carry;

            out[k] = value as u64;
            carry = value >> 64;
        }

        let mut k = i + 4;

        while carry != 0 && k < 24 {
            let value = (out[k] as u128) + carry;
            out[k] = value as u64;
            carry = value >> 64;
            k += 1;
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // params() 应返回已经正确初始化的 BN254 Fp6Params。
    // Frobenius 系数必须使用真实扩域参数。
    use crate::algebra::fields::prime_base::fp::Fp;
    use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

    /// Test field:
    ///
    /// Fp = F7
    /// Fp2 = F7[u] / (u^2 + 1)
    /// Fp6 = Fp2[v] / (v^3 - (1 + u))
    ///
    /// -1 is a quadratic non-residue in F7.
    ///
    /// Norm(1 + u) = 2, which is not a cube in F7.
    /// Therefore 1 + u is not a cube in Fp2.
    /// Thus v^3 - (1 + u) is irreducible over Fp2.
    fn params() -> Fp6Params {
        let fp_params = FieldParams::new([7, 0, 0, 0]);

        let one = Fp::one(fp_params.clone());
        let minus_one = Fp::from_u64(6, fp_params.clone());

        // Fp2 = Fp[u] / (u^2 + 1)
        //
        // u^2 = -1
        // Frobenius coefficients: [1, -1]
        let fp2_params = Fp2Params::new(
            fp_params.clone(),
            minus_one,
            [one.clone(), Fp::from_u64(6, fp_params.clone())],
        );

        // xi = 1 + u
        let xi = Fp2::new(one.clone(), one, fp2_params.clone());

        // Frobenius coefficients:
        //
        // c1[k] = xi^((7^k - 1) / 3)
        // c2[k] = xi^(2 * (7^k - 1) / 3)
        let frobenius_coeffs_c1 = std::array::from_fn(|k| {
            let p_k = 7u64.pow(k as u32);
            let exponent = (p_k - 1) / 3;

            xi.pow(&[exponent])
        });

        let frobenius_coeffs_c2 = std::array::from_fn(|k| {
            let p_k = 7u64.pow(k as u32);
            let exponent = 2 * ((p_k - 1) / 3);

            xi.pow(&[exponent])
        });

        Fp6Params {
            fp2_params,
            non_residue: xi,
            frobenius_coeffs_c1,
            frobenius_coeffs_c2,
        }
    }

    fn sample_element(params: Fp6Params) -> Fp6 {
        let fp2_params = params.fp2_params.clone();

        let c0 = Fp2::new(
            Fp::from_u64(2, fp2_params.fp_params.clone()),
            Fp::from_u64(3, fp2_params.fp_params.clone()),
            fp2_params.clone(),
        );

        let c1 = Fp2::new(
            Fp::from_u64(4, fp2_params.fp_params.clone()),
            Fp::from_u64(5, fp2_params.fp_params.clone()),
            fp2_params.clone(),
        );

        let c2 = Fp2::new(
            Fp::from_u64(6, fp2_params.fp_params.clone()),
            Fp::from_u64(1, fp2_params.fp_params.clone()),
            fp2_params,
        );

        Fp6::new(c0, c1, c2, params)
    }

    fn sample_nonzero_element(params: Fp6Params) -> Fp6 {
        let element = sample_element(params);

        assert!(!element.is_zero());

        element
    }
    #[test]
    fn test_add_sub() {
        let params = params();
        let a = sample_element(params.clone());
        let b = sample_element(params.clone());

        assert!(a.add(&b).sub(&b).equals(&a));
    }

    #[test]
    fn test_mul_square() {
        let params = params();
        let a = sample_element(params);

        assert!(a.mul(&a).equals(&a.squared()));
    }

    #[test]
    fn test_inverse() {
        let params = params();
        let a = sample_nonzero_element(params.clone());

        let inv = a.inverse().expect("nonzero element");
        let one = Fp6::one(params);

        assert!(a.mul(&inv).equals(&one));
    }

    #[test]
    fn test_frobenius_period() {
        let params = params();
        let a = sample_element(params);

        assert!(a.frobenius_map(6).equals(&a));
    }

    #[test]
    fn test_serialization() {
        let params = params();
        let a = sample_element(params.clone());

        let words = a.to_words();
        assert_eq!(words.len(), 24);

        let decoded = Fp6::from_words(&words, params).expect("valid Fp6 words");

        assert!(decoded.equals(&a));
    }

    #[test]
    fn test_sqrt() {
        let params = params();
        let a = sample_element(params);

        let square = a.squared();
        let root = square.sqrt().expect("square must have a root");

        assert!(root.squared().equals(&square));
    }
    #[test]
    fn test_bn254_fq6_frobenius_cycle() {
        use crate::algebra::curves::alt_bn128::alt_bn128_init::bn254_fp12_params;

        let params = bn254_fp12_params();
        let fp6_params = params.fp6_params;

        let fp2_params = fp6_params.fp2_params.clone();

        let x = Fp6::new(
            Fp2::from_u64(3, fp2_params.clone()),
            Fp2::from_u64(5, fp2_params.clone()),
            Fp2::from_u64(7, fp2_params),
            fp6_params,
        );

        assert!(
            x.frobenius_map(6).equals(&x),
            "BN254 Fp6 Frobenius^6 must be identity"
        );
    }
}
