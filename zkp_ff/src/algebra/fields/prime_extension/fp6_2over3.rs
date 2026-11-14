use super::super::prime_base::fp::Fp;
use super::fp2::{Fp2, Fp2Params};
use super::fp3::{Fp3, Fp3Params};

/// Fp6 = Fp3[v] / (v² - u)
///
/// Fp3 = Fp[u] / (u³ - eta)
///
/// v² = u
/// u³ = eta
///
/// Element = c0 + c1 * v
#[derive(Clone, Debug)]
pub struct Fp6Params {
    pub fp3_params: Fp3Params,

    /// Cyclotomic square 中使用的 Fp2 参数。
    ///
    /// 必须与 fp3_params 使用相同的底域和非剩余元。
    pub fp2_params: Fp2Params,

    /// Frobenius coefficients for c1.
    pub frobenius_coeffs_c1: [Fp3; 6],
}

#[derive(Clone, Debug)]
pub struct Fp6 {
    pub c0: Fp3,
    pub c1: Fp3,
    pub params: Fp6Params,
}

impl Fp6Params {
    pub fn new(
        fp3_params: Fp3Params,
        fp2_params: Fp2Params,
        frobenius_coeffs_c1: [Fp3; 6],
    ) -> Self {
        Self {
            fp3_params,
            fp2_params,
            frobenius_coeffs_c1,
        }
    }
}

impl Fp6 {
    pub fn new(c0: Fp3, c1: Fp3, params: Fp6Params) -> Self {
        Self { c0, c1, params }
    }

    pub fn zero(params: Fp6Params) -> Self {
        let fp3_params = params.fp3_params.clone();

        Self {
            c0: Fp3::zero(fp3_params.clone()),
            c1: Fp3::zero(fp3_params),
            params,
        }
    }

    pub fn one(params: Fp6Params) -> Self {
        let fp3_params = params.fp3_params.clone();

        Self {
            c0: Fp3::one(fp3_params.clone()),
            c1: Fp3::zero(fp3_params),
            params,
        }
    }

    pub fn from_u64(value: u64, params: Fp6Params) -> Self {
        let fp3_params = params.fp3_params.clone();

        Self {
            c0: Fp3::from_u64(value, fp3_params.clone()),
            c1: Fp3::zero(fp3_params),
            params,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    pub fn equals(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }

    /// 对应 libff:
    ///
    /// Fp6_2over3_model::mul_by_non_residue()
    ///
    /// u * (c0 + c1*u + c2*u²)
    ///   = eta*c2 + c0*u + c1*u²
    ///
    /// 其中 u³ = eta。
    ///
    /// 注意这里的 u 是 Fp3 的生成元，
    /// 不是 Fp6 外层的 v。
    pub fn mul_by_non_residue(&self, elem: &Fp3) -> Fp3 {
        let eta = &self.params.fp3_params.non_residue;

        Fp3::new(
            elem.c2.mul(eta),
            elem.c0.clone(),
            elem.c1.clone(),
            self.params.fp3_params.clone(),
        )
    }

    pub fn add(&self, other: &Self) -> Self {
        Self::new(
            self.c0.add(&other.c0),
            self.c1.add(&other.c1),
            self.params.clone(),
        )
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self::new(
            self.c0.sub(&other.c0),
            self.c1.sub(&other.c1),
            self.params.clone(),
        )
    }

    pub fn negate(&self) -> Self {
        Self::new(self.c0.negate(), self.c1.negate(), self.params.clone())
    }

    /// Fp 标量乘法。
    pub fn mul_by_fp(&self, scalar: &Fp) -> Self {
        Self::new(
            self.c0.mul_by_fp(scalar),
            self.c1.mul_by_fp(scalar),
            self.params.clone(),
        )
    }

    /// Karatsuba multiplication.
    ///
    /// (a + bv)(A + Bv)
    ///
    /// c0 = aA + u*bB
    /// c1 = (a+b)(A+B) - aA - bB
    ///
    /// v² = u
    pub fn mul(&self, other: &Self) -> Self {
        let a = &self.c0;
        let b = &self.c1;
        let a2 = &other.c0;
        let b2 = &other.c1;

        let aa = a.mul(a2);
        let bb = b.mul(b2);

        let beta_bb = self.mul_by_non_residue(&bb);

        let c0 = aa.add(&beta_bb);

        let c1 = a.add(b).mul(&a2.add(b2)).sub(&aa).sub(&bb);

        Self::new(c0, c1, self.params.clone())
    }

    /// 对应 libff 的 mul_by_2345。
    ///
    /// 要求:
    ///
    /// other.c0.c0 == 0
    /// other.c0.c1 == 0
    ///
    /// 即 other.c0 只有 c2 分量。
    pub fn mul_by_2345(&self, other: &Self) -> Self {
        assert!(
            other.c0.c0.is_zero(),
            "mul_by_2345: other.c0.c0 must be zero"
        );

        assert!(
            other.c0.c1.is_zero(),
            "mul_by_2345: other.c0.c1 must be zero"
        );

        let a = &self.c0;
        let b = &self.c1;
        let a2 = &other.c0;
        let b2 = &other.c1;

        let eta = &self.params.fp3_params.non_residue;
        let a2_c2 = &a2.c2;

        // 对应 C++:
        //
        // Fp3(
        //   a.c1 * A.c2 * non_residue,
        //   a.c2 * A.c2 * non_residue,
        //   a.c0 * A.c2
        // )
        let aa = Fp3::new(
            a.c1.mul(a2_c2).mul(eta),
            a.c2.mul(a2_c2).mul(eta),
            a.c0.mul(a2_c2),
            self.params.fp3_params.clone(),
        );

        let bb = b.mul(b2);
        let beta_bb = self.mul_by_non_residue(&bb);

        let c0 = aa.add(&beta_bb);

        let c1 = a.add(b).mul(&a2.add(b2)).sub(&aa).sub(&bb);

        Self::new(c0, c1, self.params.clone())
    }

    /// libff 的复数形式平方算法。
    ///
    /// ab = a*b
    ///
    /// c0 = (a+b)(a+u*b) - ab - u*ab
    /// c1 = 2ab
    pub fn squared(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;

        let ab = a.mul(b);

        let nr_b = self.mul_by_non_residue(b);
        let nr_ab = self.mul_by_non_residue(&ab);

        let c0 = a.add(b).mul(&a.add(&nr_b)).sub(&ab).sub(&nr_ab);

        let c1 = ab.add(&ab);

        Self::new(c0, c1, self.params.clone())
    }

    pub fn square(&self) -> Self {
        self.squared()
    }

    /// (a + bv)^-1 = (a - bv)/(a² - u*b²)
    pub fn inverse(&self) -> Option<Self> {
        let b2 = self.c1.squared();

        let denominator = self.c0.squared().sub(&self.mul_by_non_residue(&b2));

        let denominator_inv = denominator.inverse()?;

        let c0 = self.c0.mul(&denominator_inv);
        let c1 = self.c1.mul(&denominator_inv).negate();

        Some(Self::new(c0, c1, self.params.clone()))
    }

    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(value) => {
                *self = value;
                true
            }
            None => false,
        }
    }

    /// Frobenius map.
    ///
    /// c0 <- Frobenius(c0)
    ///
    /// c1 <- coeff[power % 6] * Frobenius(c1)
    pub fn frobenius_map(&self, power: usize) -> Self {
        let index = power % 6;

        let c0 = self.c0.frobenius_map(power);

        let c1 = self.params.frobenius_coeffs_c1[index].mul(&self.c1.frobenius_map(power));

        Self::new(c0, c1, self.params.clone())
    }

    /// Unitary inverse.
    ///
    /// 仅对 unitary subgroup 中的元素成立。
    pub fn unitary_inverse(&self) -> Self {
        Self::new(self.c0.clone(), self.c1.negate(), self.params.clone())
    }

    /// 普通二进制指数运算。
    ///
    /// exponent 使用 little-endian u64 limbs。
    pub fn pow(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for &limb in exponent {
            let mut e = limb;

            for _ in 0..64 {
                if e & 1 == 1 {
                    result = result.mul(&base);
                }

                e >>= 1;
                base = base.squared();
            }
        }

        result
    }

    /// Cyclotomic square.
    ///
    /// 只适用于 cyclotomic subgroup 中的元素。
    ///
    /// 对应 libff 的 Fp6_2over3_model::cyclotomic_squared。
    pub fn cyclotomic_squared(&self) -> Self {
        let a = Fp2::new(
            self.c0.c0.clone(),
            self.c1.c1.clone(),
            self.params.fp2_params.clone(),
        );

        let b = Fp2::new(
            self.c1.c0.clone(),
            self.c0.c2.clone(),
            self.params.fp2_params.clone(),
        );

        let c = Fp2::new(
            self.c0.c1.clone(),
            self.c1.c2.clone(),
            self.params.fp2_params.clone(),
        );

        let asq = a.square_complex();
        let bsq = b.square_complex();
        let csq = c.square_complex();

        // A_a = 3*asq.c0 - 2*a.c0
        let aa = asq.c0.sub(&a.c0).add(&asq.c0).add(&asq.c0);

        // A_b = 3*asq.c1 + 2*a.c1
        let ab = asq.c1.add(&a.c1).add(&a.c1).add(&asq.c1);

        // B_tmp = eta * csq.c1
        let eta = &self.params.fp3_params.non_residue;
        let b_tmp = csq.c1.mul(eta);

        // B_a = 3*B_tmp + 2*b.c0
        let ba = b_tmp.add(&b.c0).add(&b.c0).add(&b_tmp);

        // B_b = 3*csq.c0 - 2*b.c1
        let bb = csq.c0.sub(&b.c1).add(&csq.c0).add(&csq.c0);

        // C_a = 3*bsq.c0 - 2*c.c0
        let ca = bsq.c0.sub(&c.c0).add(&bsq.c0).add(&bsq.c0);

        // C_b = 3*bsq.c1 + 2*c.c1
        let cb = bsq.c1.add(&c.c1).add(&c.c1).add(&bsq.c1);

        let c0 = Fp3::new(aa, ca, bb, self.params.fp3_params.clone());

        let c1 = Fp3::new(ba, ab, cb, self.params.fp3_params.clone());

        Self::new(c0, c1, self.params.clone())
    }

    /// Cyclotomic exponentiation.
    ///
    /// 使用 cyclotomic square 优化。
    /// exponent 为非负 little-endian u64 limbs。
    pub fn cyclotomic_exp(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for &limb in exponent {
            let mut e = limb;

            for _ in 0..64 {
                if e & 1 == 1 {
                    result = result.mul(&base);
                }

                e >>= 1;
                base = base.cyclotomic_squared();
            }
        }

        result
    }

    /// Tonelli-Shanks square root.
    ///
    /// Fp6 的域大小为 p^6。
    ///
    /// 对 BN254 这类 256-bit 底域，
    /// p^6 - 1 最多需要 1536 bits。
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero(self.params.clone()));
        }

        let modulus = self.params.fp3_params.fp_params.modulus;

        // q = p^6
        let mut q = vec![1u64];

        for _ in 0..6 {
            q = mul_words(&q, &modulus);
        }

        // q_minus_1 = p^6 - 1
        let mut q_minus_1 = q;
        sub_one(&mut q_minus_1);

        // q - 1 = d * 2^s
        let mut d = q_minus_1.clone();
        let mut s = 0usize;

        while !is_odd(&d) {
            shr1(&mut d);
            s += 1;
        }

        // Euler criterion:
        // a^((q-1)/2) == 1
        let mut legendre_exp = q_minus_1.clone();
        shr1(&mut legendre_exp);

        let one = Self::one(self.params.clone());

        let legendre = self.pow(&legendre_exp);

        if !legendre.equals(&one) {
            return None;
        }

        // Search a quadratic non-residue.
        //
        // z = 2, 3, 4, ...
        let mut z = None;

        for candidate in 2u64..4096 {
            // Do not search only in the base field Fp: every nonzero
            // element of Fp embeds as a square in the even-degree
            // extension Fp6. Search in Fp6 \ Fp instead, using k + v.
            let value = Self::new(
                Fp3::from_u64(candidate, self.params.fp3_params.clone()),
                Fp3::one(self.params.fp3_params.clone()),
                self.params.clone(),
            );

            if !value.pow(&legendre_exp).equals(&one) {
                z = Some(value);
                break;
            }
        }

        let z = z?;

        let mut c = z.pow(&d);

        // x = a^((d+1)/2)
        let mut d_plus_one = d.clone();
        add_one(&mut d_plus_one);
        shr1(&mut d_plus_one);

        let mut x = self.pow(&d_plus_one);

        // t = a^d
        let mut t = self.pow(&d);

        let mut m = s;

        while !t.equals(&one) {
            let mut i = 0usize;
            let mut t2 = t.clone();

            loop {
                if i >= m {
                    return None;
                }

                t2 = t2.squared();
                i += 1;

                if t2.equals(&one) {
                    break;
                }
            }

            // b = c^(2^(m-i-1))
            let mut b = c.clone();

            for _ in 0..(m - i - 1) {
                b = b.squared();
            }

            x = x.mul(&b);

            let b2 = b.squared();

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

    /// c0 || c1
    ///
    /// Fp3 = 12 words
    /// Fp6 = 24 words
    pub fn to_words(&self) -> Vec<u64> {
        let mut words = self.c0.to_words();
        words.extend(self.c1.to_words());
        words
    }

    pub fn from_words(words: &[u64], params: Fp6Params) -> Option<Self> {
        if words.len() != 24 {
            return None;
        }

        let c0 = Fp3::from_words(&words[0..12], params.fp3_params.clone())?;

        let c1 = Fp3::from_words(&words[12..24], params.fp3_params.clone())?;

        Some(Self::new(c0, c1, params))
    }
}

// -----------------------------------------------------------------------------
// Big integer helpers for Tonelli-Shanks.
//
// Little-endian u64 limbs.
// -----------------------------------------------------------------------------

fn normalize_words(words: &mut Vec<u64>) {
    while words.len() > 1 && words.last() == Some(&0) {
        words.pop();
    }
}

fn is_zero_words(words: &[u64]) -> bool {
    words.iter().all(|x| *x == 0)
}

fn is_odd(words: &[u64]) -> bool {
    words.first().copied().unwrap_or(0) & 1 == 1
}

fn shr1(words: &mut Vec<u64>) {
    let mut carry = 0u64;

    for i in (0..words.len()).rev() {
        let next_carry = words[i] & 1;

        words[i] = (words[i] >> 1) | (carry << 63);
        carry = next_carry;
    }

    normalize_words(words);
}

fn add_one(words: &mut Vec<u64>) {
    let mut carry = true;

    for word in words.iter_mut() {
        if !carry {
            break;
        }

        let (value, overflow) = word.overflowing_add(1);
        *word = value;
        carry = overflow;
    }

    if carry {
        words.push(1);
    }
}

fn sub_one(words: &mut Vec<u64>) {
    let mut borrow = true;

    for word in words.iter_mut() {
        if !borrow {
            break;
        }

        let (value, underflow) = word.overflowing_sub(1);
        *word = value;
        borrow = underflow;
    }

    debug_assert!(!borrow);
    normalize_words(words);
}

/// Multiplication of arbitrary little-endian limb arrays.
fn mul_words(a: &[u64], b: &[u64]) -> Vec<u64> {
    if is_zero_words(a) || is_zero_words(b) {
        return vec![0];
    }

    let mut out = vec![0u64; a.len() + b.len()];

    for i in 0..a.len() {
        let mut carry = 0u128;

        for j in 0..b.len() {
            let uv = (a[i] as u128) * (b[j] as u128) + (out[i + j] as u128) + carry;

            out[i + j] = uv as u64;
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

    normalize_words(&mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

    /// Small test field:
    ///
    /// Fp = F7
    /// Fp3 = F7[u] / (u^3 - 3)
    /// Fp6 = Fp3[v] / (v^2 - u)
    ///
    /// 3 is both a non-cube and a non-square in F7:
    /// - u^3 - 3 is irreducible over F7;
    /// - Norm_Fp3/Fp(u) = 3 is a non-square, so u is a non-square
    ///   in Fp3 and v^2 - u is irreducible over Fp3.
    fn test_fp6_params() -> Fp6Params {
        let fp_params = FieldParams::new([7, 0, 0, 0]);

        let one = Fp::one(fp_params.clone());
        let two = Fp::from_u64(2, fp_params.clone());
        let three = Fp::from_u64(3, fp_params.clone());
        let four = Fp::from_u64(4, fp_params.clone());
        let five = Fp::from_u64(5, fp_params.clone());
        let six = Fp::from_u64(6, fp_params.clone());

        // Fp3 = Fp[u] / (u^3 - 3).
        //
        // u^7 = 3^2 * u = 2u
        // u^14 = 3^4 * u^2 = 4u^2
        let fp3_params = Fp3Params::new(
            fp_params.clone(),
            three.clone(),
            [one.clone(), two.clone(), four.clone()],
            [one.clone(), four.clone(), two.clone()],
        );

        // Fp2 is only needed by the cyclotomic-square implementation.
        // Use the valid quadratic extension Fp[w] / (w^2 - 3).
        // Since 7 == 3 (mod 4), w^7 / w = -1.
        let fp2_params = Fp2Params::new(fp_params.clone(), three, [one.clone(), six]);

        // For v^2 = u, the outer Frobenius coefficient is
        // u^((7^i - 1) / 2). With u^3 = 3, these are the following
        // base-field constants for i = 0..5:
        //
        // [1, 3, 3^8, 3^57, 3^400, 3^2801] mod 7
        // = [1, 3, 2, 6, 4, 5].
        let frobenius_coeffs_c1 = [
            Fp3::from_u64(1, fp3_params.clone()),
            Fp3::from_u64(3, fp3_params.clone()),
            Fp3::from_u64(2, fp3_params.clone()),
            Fp3::from_u64(6, fp3_params.clone()),
            Fp3::from_u64(4, fp3_params.clone()),
            Fp3::from_u64(5, fp3_params.clone()),
        ];

        Fp6Params::new(fp3_params, fp2_params, frobenius_coeffs_c1)
    }

    fn test_fp6_element() -> Fp6 {
        let params = test_fp6_params();
        let fp = params.fp3_params.fp_params.clone();
        let fp3 = params.fp3_params.clone();

        let c0 = Fp3::new(
            Fp::from_u64(1, fp.clone()),
            Fp::from_u64(2, fp.clone()),
            Fp::from_u64(3, fp.clone()),
            fp3.clone(),
        );

        let c1 = Fp3::new(
            Fp::from_u64(4, fp.clone()),
            Fp::from_u64(5, fp.clone()),
            Fp::from_u64(6, fp),
            fp3,
        );

        Fp6::new(c0, c1, params)
    }

    #[test]
    fn test_fp6_zero_and_one() {
        let params = test_fp6_params();

        assert!(Fp6::zero(params.clone()).is_zero());
        assert!(Fp6::one(params).is_one());
    }

    #[test]
    fn test_fp6_add_sub() {
        let a = test_fp6_element();
        let b = test_fp6_element().negate();

        assert!(a.add(&b).is_zero());
        assert!(a.add(&b).sub(&b).equals(&a));
    }

    #[test]
    fn test_fp6_square_equals_mul() {
        let a = test_fp6_element();
        assert!(a.squared().equals(&a.mul(&a)));
        assert!(a.square().equals(&a.mul(&a)));
    }

    #[test]
    fn test_fp6_inverse() {
        let a = test_fp6_element();

        let inv = a.inverse().expect("nonzero element must be invertible");
        assert!(a.mul(&inv).is_one());
    }

    #[test]
    fn test_fp6_zero_inverse() {
        let zero = Fp6::zero(test_fp6_params());

        assert!(zero.inverse().is_none());

        let mut value = zero;
        assert!(!value.invert());
        assert!(value.is_zero());
    }

    #[test]
    fn test_fp6_frobenius_period() {
        let a = test_fp6_element();

        assert!(a.frobenius_map(6).equals(&a));

        let mut value = a.clone();
        for _ in 0..6 {
            value = value.frobenius_map(1);
        }
        assert!(value.equals(&a));
    }

    #[test]
    fn test_fp6_frobenius_additive_and_multiplicative() {
        let a = test_fp6_element();
        let b = test_fp6_element().add(&Fp6::one(test_fp6_params()));

        assert!(
            a.add(&b)
                .frobenius_map(1)
                .equals(&a.frobenius_map(1).add(&b.frobenius_map(1)))
        );

        assert!(
            a.mul(&b)
                .frobenius_map(1)
                .equals(&a.frobenius_map(1).mul(&b.frobenius_map(1)))
        );
    }

    #[test]
    fn test_fp6_serialization() {
        let a = test_fp6_element();

        let words = a.to_words();
        assert_eq!(words.len(), 24);

        let b = Fp6::from_words(&words, a.params.clone()).unwrap();

        assert!(a.equals(&b));
    }

    #[test]
    fn test_fp6_invalid_serialization_length() {
        let params = test_fp6_params();

        assert!(Fp6::from_words(&[1u64; 23], params).is_none());
    }

    #[test]
    fn test_fp6_sqrt_roundtrip() {
        let a = test_fp6_element();
        let square = a.squared();

        let root = square.sqrt().expect("a square must have a square root");

        assert!(root.squared().equals(&square));
    }

    #[test]
    fn test_fp6_sqrt_zero() {
        let zero = Fp6::zero(test_fp6_params());

        let root = zero.sqrt().unwrap();
        assert!(root.is_zero());
    }
}
