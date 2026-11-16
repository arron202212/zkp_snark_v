
use num_complex::Complex64;
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};

use std::ops::{
    Add, AddAssign,
    Sub, SubAssign,
    Mul, MulAssign,
    Neg,
};

/// PI constant.
pub const PI: f64 = 3.141592653589793238460264338328;

/// Complex floating-point type.
///
/// Equivalent to libff::Double.
#[derive(Clone, Copy, Debug, Default)]
pub struct Double {
    pub val: Complex64,
}

/// Operation counters.
#[derive(Clone, Copy, Debug, Default)]
pub struct DoubleOperationCounts {
    pub add_cnt: u64,
    pub sub_cnt: u64,
    pub mul_cnt: u64,
    pub inv_cnt: u64,
}

impl Double {
    /// Construct from real and imaginary parts.
    pub fn new(real: f64, imag: f64) -> Self {
        Self {
            val: Complex64::new(real, imag),
        }
    }

    /// Construct from a real number.
    pub fn from_real(real: f64) -> Self {
        Self::new(real, 0.0)
    }

    /// Construct from a complex number.
    pub fn from_complex(val: Complex64) -> Self {
        Self { val }
    }

    /// Real part.
    pub fn real(&self) -> f64 {
        self.val.re
    }

    /// Imaginary part.
    pub fn imag(&self) -> f64 {
        self.val.im
    }

    /// Multiplicative inverse.
    pub fn inverse(&self) -> Self {
        Self::from_complex(Complex64::new(1.0, 0.0) / self.val)
    }

    /// Square.
    pub fn squared(&self) -> Self {
        Self::from_complex(self.val * self.val)
    }

    /// Exponentiation by usize.
    pub fn pow(&self, power: usize) -> Self {
        Self::from_complex(self.val.powu(power as u32))
    }

    /// Exponentiation by arbitrary precision integer.
    pub fn pow_bigint(&self, power: &BigUint) -> Self {
        let mut exp = power.clone();
        let mut base = self.val;
        let mut result = Complex64::new(1.0, 0.0);

        while !exp.is_zero() {
            if (&exp & BigUint::from(1u8)) == BigUint::from(1u8) {
                result *= base;
            }

            base *= base;
            exp >>= 1;
        }

        Self::from_complex(result)
    }

    /// Convert real part to BigUint.
    pub fn as_bigint(&self) -> BigUint {
        if !self.val.re.is_finite() || self.val.re < 0.0 {
            return BigUint::zero();
        }

        BigUint::from(self.val.re as usize)
    }

    /// Convert real part to unsigned long equivalent.
    pub fn as_ulong(&self) -> u64 {
        if !self.val.re.is_finite() || self.val.re < 0.0 {
            return 0;
        }

        self.val.re.round() as u64
    }

    /// One.
    pub fn one() -> Self {
        Self::from_real(1.0)
    }

    /// Zero.
    pub fn zero() -> Self {
        Self::from_real(0.0)
    }

    /// Random element in [0, 1000].
    pub fn random_element() -> Self {
        use rand::Rng;

        let mut rng = rand::thread_rng();
        Self::from_real(rng.gen_range(0..=1000) as f64)
    }

    /// Geometric generator.
    pub fn geometric_generator() -> Self {
        Self::from_real(2.0)
    }

    /// Arithmetic generator.
    pub fn arithmetic_generator() -> Self {
        Self::from_real(1.0)
    }

    /// Multiplicative generator.
    pub fn multiplicative_generator() -> Self {
        Self::from_real(2.0)
    }

    /// Compare real parts.
    pub fn real_lt(&self, other: &Self) -> bool {
        self.val.re < other.val.re
    }

    /// Compare real parts.
    pub fn real_gt(&self, other: &Self) -> bool {
        self.val.re > other.val.re
    }
}

impl Add for Double {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::from_complex(self.val + rhs.val)
    }
}

impl Sub for Double {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self::from_complex(self.val - rhs.val)
    }
}

impl Mul for Double {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::from_complex(self.val * rhs.val)
    }
}

impl Neg for Double {
    type Output = Self;

    fn neg(self) -> Self {
        Self::from_complex(-self.val)
    }
}

impl AddAssign for Double {
    fn add_assign(&mut self, rhs: Self) {
        self.val += rhs.val;
    }
}

impl SubAssign for Double {
    fn sub_assign(&mut self, rhs: Self) {
        self.val -= rhs.val;
    }
}

impl MulAssign for Double {
    fn mul_assign(&mut self, rhs: Self) {
        self.val *= rhs.val;
    }
}

/// Approximate equality, matching libff's tolerance.
impl PartialEq for Double {
    fn eq(&self, other: &Self) -> bool {
        (self.val.re - other.val.re).abs() < 0.000001
            && (self.val.im - other.val.im).abs() < 0.000001
    }
}

impl From<f64> for Double {
    fn from(value: f64) -> Self {
        Self::from_real(value)
    }
}


// #[cfg(feature = "profile-op-counts")]
// pub mod counters {
//     use std::sync::atomic::AtomicU64;

//     pub static ADD_CNT: AtomicU64 = AtomicU64::new(0);
//     pub static SUB_CNT: AtomicU64 = AtomicU64::new(0);
//     pub static MUL_CNT: AtomicU64 = AtomicU64::new(0);
//     pub static INV_CNT: AtomicU64 = AtomicU64::new(0);
// }

// #[cfg(feature = "profile-op-counts")]
// counters::ADD_CNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constructors() {
        let x = Double::new(3.0, 4.0);

        assert_eq!(x.real(), 3.0);
        assert_eq!(x.imag(), 4.0);
    }

    #[test]
    fn test_add() {
        let a = Double::new(1.0, 2.0);
        let b = Double::new(3.0, 4.0);

        assert_eq!(a + b, Double::new(4.0, 6.0));
    }

    #[test]
    fn test_sub() {
        let a = Double::new(5.0, 6.0);
        let b = Double::new(2.0, 3.0);

        assert_eq!(a - b, Double::new(3.0, 3.0));
    }

    #[test]
    fn test_mul() {
        let a = Double::new(1.0, 2.0);
        let b = Double::new(3.0, 4.0);

        // (1+2i)(3+4i) = -5+10i
        assert_eq!(a * b, Double::new(-5.0, 10.0));
    }

    #[test]
    fn test_inverse() {
        let a = Double::new(3.0, 4.0);
        let result = a * a.inverse();

        assert_eq!(result, Double::one());
    }

    #[test]
    fn test_squared() {
        let a = Double::new(3.0, 4.0);

        assert_eq!(
            a.squared(),
            Double::new(-7.0, 24.0)
        );
    }

    #[test]
    fn test_pow() {
        let a = Double::new(1.0, 1.0);

        assert_eq!(
            a.pow(2),
            Double::new(0.0, 2.0)
        );

        assert_eq!(a.pow(0), Double::one());
    }

    #[test]
    fn test_generators() {
        assert_eq!(Double::zero(), Double::new(0.0, 0.0));
        assert_eq!(Double::one(), Double::new(1.0, 0.0));

        assert_eq!(
            Double::geometric_generator(),
            Double::new(2.0, 0.0)
        );

        assert_eq!(
            Double::arithmetic_generator(),
            Double::new(1.0, 0.0)
        );
    }
}