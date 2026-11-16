use core::cmp::Ordering;

use crate::algebra::fields::prime_base::fp_montgomery::{
    FieldParams, LimbsN, add_mod, add_one, cmp, from_montgomery, montgomery_mul, shr1, sub_mod,
    sub_one, to_montgomery,
};

/// Prime field element.
///
/// 默认 N=4，兼容原有 BN254 API。
///
/// Montgomery representation:
///```text
///     mont_repr = x * R mod p
///```
/// 其中 R = 2^(64*N)。
#[derive(Clone, Debug)]
pub struct Fp<const N: usize = 4> {
    pub mont_repr: LimbsN<N>,
    pub params: FieldParams<N>,
}

impl<const N: usize> PartialEq for Fp<N> {
    fn eq(&self, other: &Self) -> bool {
        self.params.modulus == other.params.modulus && self.mont_repr == other.mont_repr
    }
}

impl<const N: usize> Eq for Fp<N> {}

impl<const N: usize> Fp<N> {
    // ---------------------------------------------------------------------
    // Constructors
    // ---------------------------------------------------------------------

    /// Construct from a canonical integer.
    ///
    /// Requires value < modulus.
    pub fn new(value: LimbsN<N>, params: FieldParams<N>) -> Self {
        assert!(N > 0, "limb count must be greater than zero");

        assert!(
            cmp(&value, &params.modulus) == Ordering::Less,
            "Value must be smaller than modulus"
        );

        let mont_repr = to_montgomery(&value, &params);

        Self { mont_repr, params }
    }

    /// Construct from u64.
    pub fn from_u64(value: u64, params: FieldParams<N>) -> Self {
        assert!(N > 0);

        let mut words = [0u64; N];
        words[0] = value;

        Self::new(words, params)
    }

    /// Additive identity.
    pub fn zero(params: FieldParams<N>) -> Self {
        Self {
            mont_repr: [0u64; N],
            params,
        }
    }

    /// Multiplicative identity.
    pub fn one(params: FieldParams<N>) -> Self {
        Self {
            mont_repr: params.r,
            params,
        }
    }

    // ---------------------------------------------------------------------
    // Basic operations
    // ---------------------------------------------------------------------

    pub fn is_zero(&self) -> bool {
        self.mont_repr.iter().all(|&x| x == 0)
    }

    fn check_field(&self, other: &Self) {
        assert_eq!(
            self.params.modulus, other.params.modulus,
            "Different finite fields"
        );
    }

    /// Convert Montgomery representation to canonical integer.
    pub fn to_canonical(&self) -> LimbsN<N> {
        from_montgomery(&self.mont_repr, &self.params)
    }

    /// Addition.
    pub fn add(&self, other: &Self) -> Self {
        self.check_field(other);

        Self {
            mont_repr: add_mod(self.mont_repr, other.mont_repr, self.params.modulus),
            params: self.params.clone(),
        }
    }

    /// Subtraction.
    pub fn sub(&self, other: &Self) -> Self {
        self.check_field(other);

        Self {
            mont_repr: sub_mod(self.mont_repr, other.mont_repr, self.params.modulus),
            params: self.params.clone(),
        }
    }

    /// Multiplication.
    pub fn mul(&self, other: &Self) -> Self {
        self.check_field(other);

        Self {
            mont_repr: montgomery_mul(self.mont_repr, other.mont_repr, &self.params),
            params: self.params.clone(),
        }
    }

    /// Squaring.
    pub fn square(&self) -> Self {
        Self {
            mont_repr: montgomery_mul(self.mont_repr, self.mont_repr, &self.params),
            params: self.params.clone(),
        }
    }

    /// Negation.
    pub fn negate(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }

        Self {
            mont_repr: sub_mod([0u64; N], self.mont_repr, self.params.modulus),
            params: self.params.clone(),
        }
    }

    // ---------------------------------------------------------------------
    // Exponentiation
    // ---------------------------------------------------------------------

    /// Modular exponentiation.
    ///
    /// Exponent is an N-limb little-endian integer.
    ///
    /// Scans all 64*N bits.
    pub fn pow(&self, exponent: LimbsN<N>) -> Self {
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

    // ---------------------------------------------------------------------
    // Inversion
    // ---------------------------------------------------------------------

    /// Multiplicative inverse.
    ///
    /// x^(-1) = x^(p-2) mod p
    ///
    /// Assumes modulus is prime.
    pub fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let exponent = sub_one(&self.params.modulus);

        // p - 2
        let exponent = sub_one(&exponent);

        Some(self.pow(exponent))
    }

    /// In-place inverse.
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

    // ---------------------------------------------------------------------
    // Square root
    // ---------------------------------------------------------------------

    /// Tonelli-Shanks square root.
    ///
    /// Returns Some(root) if a square root exists.
    ///
    /// The returned root is not guaranteed to be the smaller root.
    ///
    /// Assumes modulus is an odd prime.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(self.clone());
        }

        let p = self.params.modulus;

        // p - 1
        let p_minus_one = sub_one(&p);

        // (p - 1) / 2
        let legendre_exp = shr1(&p_minus_one);

        let one = Self::one(self.params.clone());

        // Euler's criterion.
        let legendre = self.pow(legendre_exp);

        if legendre != one {
            return None;
        }

        // Factor:
        //
        // p - 1 = q * 2^s
        //
        // where q is odd.
        let mut q = p_minus_one;
        let mut s = 0usize;

        while q[0] & 1 == 0 {
            q = shr1(&q);
            s += 1;
        }

        // p = 3 mod 4.
        if s == 1 {
            // (p + 1) / 4
            let p_plus_one = add_one(&p);
            let exponent = shr1(&shr1(&p_plus_one));

            return Some(self.pow(exponent));
        }

        // Find quadratic non-residue z.
        let mut z_value = 2u64;

        let z = loop {
            let mut words = [0u64; N];
            words[0] = z_value;

            // Avoid constructing a non-canonical field element.
            if cmp(&words, &p) != Ordering::Less {
                return None;
            }

            let candidate = Self::new(words, self.params.clone());

            let symbol = candidate.pow(legendre_exp);

            if symbol != one {
                break candidate;
            }

            z_value = z_value.checked_add(1)?;
        };

        // c = z^q
        let mut c = z.pow(q);

        // x = a^((q+1)/2)
        let q_plus_one = add_one(&q);
        let x_exp = shr1(&q_plus_one);

        let mut x = self.pow(x_exp);

        // t = a^q
        let mut t = self.pow(q);

        let mut m = s;

        while t != one {
            let mut i = 0usize;
            let mut t2 = t.clone();

            while t2 != one {
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

    // ---------------------------------------------------------------------
    // Frobenius
    // ---------------------------------------------------------------------

    /// Frobenius map over a prime field.
    ///
    /// a -> a^(p^power) = a
    pub fn frobenius_map(&self, _power: usize) -> Self {
        self.clone()
    }

    // ---------------------------------------------------------------------
    // Serialization
    // ---------------------------------------------------------------------

    /// Convert to canonical little-endian words.
    pub fn to_words(&self) -> Vec<u64> {
        self.to_canonical().to_vec()
    }

    /// Read canonical little-endian words.
    ///
    /// Returns None if:
    /// - input length != N
    /// - input >= modulus
    pub fn from_words(words: &[u64], params: FieldParams<N>) -> Option<Self> {
        if words.len() != N {
            return None;
        }

        let mut value = [0u64; N];
        value.copy_from_slice(words);

        if cmp(&value, &params.modulus) != Ordering::Less {
            return None;
        }

        Some(Self::new(value, params))
    }
}

// -------------------------------------------------------------------------
// Type aliases
// -------------------------------------------------------------------------

/// BN254 field element.
pub type Bn254Fp = Fp<4>;

/// MNT4-298 base field element.
pub type Mnt4Fq = Fp<5>;

/// MNT4-298 scalar field element.
pub type Mnt4Fr = Fp<5>;
