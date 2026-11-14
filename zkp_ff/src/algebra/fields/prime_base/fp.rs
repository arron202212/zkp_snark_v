use crate::algebra::fields::prime_base::fp_montgomery::{
    FieldParams, Limbs, add_mod, cmp, montgomery_mul, sub_mod,
};

#[derive(Clone, Debug)]
pub struct Fp {
    /// Montgomery representation: x * R mod p.
    pub mont_repr: Limbs,

    pub params: FieldParams,
}

impl PartialEq for Fp {
    fn eq(&self, other: &Self) -> bool {
        self.params.modulus == other.params.modulus && self.mont_repr == other.mont_repr
    }
}

impl Eq for Fp {}

impl Fp {
    pub fn new(value: Limbs, params: FieldParams) -> Self {
        assert!(
            cmp(&value, &params.modulus) == core::cmp::Ordering::Less,
            "Value must be smaller than modulus"
        );

        let mont_repr = montgomery_mul(value, params.r2, &params);

        Self { mont_repr, params }
    }

    pub fn from_u64(value: u64, params: FieldParams) -> Self {
        Self::new([value, 0, 0, 0], params)
    }

    pub fn zero(params: FieldParams) -> Self {
        Self {
            mont_repr: [0; 4],
            params,
        }
    }

    pub fn one(params: FieldParams) -> Self {
        Self {
            mont_repr: params.r,
            params,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.mont_repr == [0; 4]
    }

    fn check_field(&self, other: &Self) {
        assert_eq!(
            self.params.modulus, other.params.modulus,
            "Different finite fields"
        );
    }

    /// Convert Montgomery representation into canonical integer.
    pub fn to_canonical(&self) -> Limbs {
        montgomery_mul(self.mont_repr, [1, 0, 0, 0], &self.params)
    }

    /// Addition in Fp.
    pub fn add(&self, other: &Self) -> Self {
        self.check_field(other);

        Self {
            mont_repr: add_mod(self.mont_repr, other.mont_repr, self.params.modulus),
            params: self.params.clone(),
        }
    }

    /// Subtraction in Fp.
    pub fn sub(&self, other: &Self) -> Self {
        self.check_field(other);

        Self {
            mont_repr: sub_mod(self.mont_repr, other.mont_repr, self.params.modulus),
            params: self.params.clone(),
        }
    }

    /// Multiplication in Fp.
    pub fn mul(&self, other: &Self) -> Self {
        self.check_field(other);

        Self {
            mont_repr: montgomery_mul(self.mont_repr, other.mont_repr, &self.params),
            params: self.params.clone(),
        }
    }

    /// Squaring.
    pub fn square(&self) -> Self {
        self.mul(self)
    }

    /// Negation: -x mod p.
    pub fn negate(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }

        Self {
            mont_repr: sub_mod([0; 4], self.mont_repr, self.params.modulus),
            params: self.params.clone(),
        }
    }

    /// Modular exponentiation.
    ///
    /// Exponent is a 256-bit integer in little-endian limbs.
    pub fn pow(&self, exponent: Limbs) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for limb in exponent.iter() {
            let mut e = *limb;

            for _ in 0..64 {
                if e & 1 == 1 {
                    result = result.mul(&base);
                }

                base = base.square();
                e >>= 1;
            }
        }

        result
    }

    /// Multiplicative inverse.
    ///
    /// Fermat's little theorem:
    ///
    /// x^-1 = x^(p-2) mod p
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let exponent = sub_small(self.params.modulus, 2);

        Some(self.pow(exponent))
    }

    /// In-place inverse, matching libff's invert() behavior.
    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(value) => {
                *self = value;
                true
            }
            None => false,
        }
    }

    /// Tonelli-Shanks square root.
    ///
    /// Returns Some(root) when a square root exists.
    ///
    /// The returned root is not guaranteed to be the smaller root.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(self.clone());
        }

        let p = self.params.modulus;

        let p_minus_one = sub_small(p, 1);

        // Legendre exponent = (p - 1) / 2.
        let legendre_exp = shr(p_minus_one, 1);

        let legendre = self.pow(legendre_exp);

        if legendre != Self::one(self.params.clone()) {
            return None;
        }

        // Factor p - 1 = q * 2^s, q odd.
        let mut q = p_minus_one;
        let mut s = 0usize;

        while q[0] & 1 == 0 {
            q = shr(q, 1);
            s += 1;
        }

        // p = 3 mod 4 shortcut.
        if s == 1 {
            let exponent = shr(add_small(p, 1), 2);

            return Some(self.pow(exponent));
        }

        // Find quadratic non-residue z.
        let mut z_value = 2u64;

        let z = loop {
            let candidate = Self::from_u64(z_value, self.params.clone());

            let symbol = candidate.pow(legendre_exp);

            if symbol != Self::one(self.params.clone()) {
                break candidate;
            }

            z_value += 1;
        };

        let mut c = z.pow(q);
        let mut x = self.pow(shr(add_small(q, 1), 1));
        let mut t = self.pow(q);
        let mut m = s;

        while t != Self::one(self.params.clone()) {
            let mut i = 0usize;
            let mut t2 = t.clone();

            while t2 != Self::one(self.params.clone()) {
                t2 = t2.square();
                i += 1;

                if i >= m {
                    return None;
                }
            }

            let mut b = c.clone();

            for _ in 0..(m - i - 1) {
                b = b.square();
            }

            x = x.mul(&b);
            let b2 = b.square();

            t = t.mul(&b2);
            c = b2;
            m = i;
        }

        Some(x)
    }

    /// Frobenius map.
    ///
    /// For the prime field Fp, Frobenius is the identity.
    pub fn frobenius_map(&self, _power: usize) -> Self {
        self.clone()
    }

    /// Convert to four 64-bit canonical words.
    pub fn to_words(&self) -> Vec<u64> {
        self.to_canonical().to_vec()
    }

    /// Read four canonical words.
    ///
    /// Returns false if the input is >= modulus.
    pub fn from_words(words: &[u64], params: FieldParams) -> Option<Self> {
        if words.len() != 4 {
            return None;
        }

        let value: Limbs = [words[0], words[1], words[2], words[3]];

        if cmp(&value, &params.modulus) != core::cmp::Ordering::Less {
            return None;
        }

        Some(Self::new(value, params))
    }
}

/// Subtract a small integer from a 256-bit integer.
fn sub_small(mut a: Limbs, value: u64) -> Limbs {
    let (first, mut borrow) = a[0].overflowing_sub(value);
    a[0] = first;

    for limb in a.iter_mut().skip(1) {
        if !borrow {
            break;
        }

        let (result, next) = limb.overflowing_sub(1);
        *limb = result;
        borrow = next;
    }

    assert!(!borrow);
    a
}

/// Add a small integer.
fn add_small(mut a: Limbs, value: u64) -> Limbs {
    let (first, mut carry) = a[0].overflowing_add(value);
    a[0] = first;

    for limb in a.iter_mut().skip(1) {
        if !carry {
            break;
        }

        let (result, next) = limb.overflowing_add(1);
        *limb = result;
        carry = next;
    }

    assert!(!carry);
    a
}

/// Divide a 256-bit integer by two.
fn shr(mut a: Limbs, bits: usize) -> Limbs {
    assert!(bits < 64);

    if bits == 0 {
        return a;
    }

    for i in 0..3 {
        a[i] = (a[i] >> bits) | (a[i + 1] << (64 - bits));
    }

    a[3] >>= bits;
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fr_params() -> FieldParams {
        // BN254 scalar field modulus.
        FieldParams::new([
            0x43e1_f593_f000_0001,
            0x2833_e848_79b9_7091,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ])
    }

    fn canonical(value: &Fp) -> u64 {
        let words = value.to_canonical();
        assert_eq!(words[1..], [0, 0, 0]);
        words[0]
    }

    #[test]
    fn test_zero_one() {
        let params = fr_params();

        let zero = Fp::zero(params.clone());
        let one = Fp::one(params.clone());

        assert!(zero.is_zero());
        assert_eq!(canonical(&one), 1);
    }

    #[test]
    fn test_add() {
        let params = fr_params();

        let a = Fp::from_u64(7, params.clone());
        let b = Fp::from_u64(9, params);

        assert_eq!(canonical(&a.add(&b)), 16);
    }

    #[test]
    fn test_sub() {
        let params = fr_params();

        let a = Fp::from_u64(9, params.clone());
        let b = Fp::from_u64(7, params);

        assert_eq!(canonical(&a.sub(&b)), 2);
    }

    #[test]
    fn test_sub_underflow() {
        let params = fr_params();

        let a = Fp::from_u64(3, params.clone());
        let b = Fp::from_u64(7, params.clone());

        let result = a.sub(&b);

        let expected = sub_mod([3, 0, 0, 0], [7, 0, 0, 0], params.modulus);

        assert_eq!(result.to_canonical(), expected);
    }

    #[test]
    fn test_mul() {
        let params = fr_params();

        let a = Fp::from_u64(7, params.clone());
        let b = Fp::from_u64(9, params);

        assert_eq!(canonical(&a.mul(&b)), 63);
    }

    #[test]
    fn test_square() {
        let params = fr_params();

        let a = Fp::from_u64(12, params);

        assert_eq!(canonical(&a.square()), 144);
    }

    #[test]
    fn test_negate() {
        let params = fr_params();

        let a = Fp::from_u64(5, params.clone());
        let result = a.add(&a.negate());

        assert!(result.is_zero());
    }

    #[test]
    fn test_inverse() {
        let params = fr_params();

        let a = Fp::from_u64(7, params.clone());

        let inverse = a.inverse().unwrap();

        let result = a.mul(&inverse);

        assert_eq!(result, Fp::one(params));
    }

    #[test]
    fn test_zero_inverse() {
        let params = fr_params();

        let zero = Fp::zero(params);

        assert!(zero.inverse().is_none());
    }

    #[test]
    fn test_pow() {
        let params = fr_params();

        let a = Fp::from_u64(3, params);

        let result = a.pow([4, 0, 0, 0]);

        assert_eq!(canonical(&result), 81);
    }

    #[test]
    fn test_sqrt() {
        let params = fr_params();

        let a = Fp::from_u64(9, params);

        let root = a.sqrt().unwrap();

        assert_eq!(root.square(), a);
    }

    #[test]
    fn test_sqrt_zero() {
        let params = fr_params();

        let zero = Fp::zero(params);

        assert_eq!(zero.sqrt().unwrap(), zero);
    }

    #[test]
    fn test_serialization() {
        let params = fr_params();

        let a = Fp::from_u64(123456, params.clone());

        let words = a.to_words();

        assert_eq!(words.len(), 4);

        let b = Fp::from_words(&words, params).unwrap();

        assert_eq!(a, b);
    }

    #[test]
    fn test_invalid_serialization() {
        let params = fr_params();

        let words = params.modulus;

        assert!(Fp::from_words(&words, params).is_none());
    }

    #[test]
    fn test_frobenius() {
        let params = fr_params();

        let a = Fp::from_u64(123, params);

        assert_eq!(a, a.frobenius_map(1));
        assert_eq!(a, a.frobenius_map(2));
    }
}
