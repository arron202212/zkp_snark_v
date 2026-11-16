/// 多项式运算需要的有限域接口。
///
/// 系数按低次到高次排列：
/// [a0, a1, a2] 表示 a0 + a1*x + a2*x^2。

pub trait PolynomialField: Clone + PartialEq + std::fmt::Debug {
    type Params: Clone + std::fmt::Debug;

    fn zero(params: &Self::Params) -> Self;

    fn one(params: &Self::Params) -> Self;

    fn from_u64(value: u64, params: &Self::Params) -> Self;

    fn one_like(&self) -> Self;
    fn from_u64_like(&self, value: u64) -> Self;
    fn add(&self, rhs: &Self) -> Self;

    fn sub(&self, rhs: &Self) -> Self;

    fn mul(&self, rhs: &Self) -> Self;

    fn inverse(&self) -> Option<Self>;

    fn is_zero(&self) -> bool;

    fn pow_usize(&self, exponent: usize) -> Self;
}
use zkp_ff::algebra::fields::prime_base::fp::Fp;
use zkp_ff::algebra::fields::prime_base::fp_montgomery::FieldParams;

impl PolynomialField for Fp {
    type Params = FieldParams;

    fn zero(params: &Self::Params) -> Self {
        Fp::zero(params.clone())
    }

    fn one(params: &Self::Params) -> Self {
        Fp::one(params.clone())
    }

    fn from_u64(value: u64, params: &Self::Params) -> Self {
        Fp::from_u64(value, params.clone())
    }

    fn add(&self, rhs: &Self) -> Self {
        Fp::add(self, rhs)
    }

    fn sub(&self, rhs: &Self) -> Self {
        Fp::sub(self, rhs)
    }

    fn mul(&self, rhs: &Self) -> Self {
        Fp::mul(self, rhs)
    }
    fn one_like(&self) -> Self {
        Fp::one(self.params.clone())
    }
    fn from_u64_like(&self, value: u64) -> Self {
        Fp::from_u64(value, self.params.clone())
    }
    fn inverse(&self) -> Option<Self> {
        Fp::inverse(self)
    }

    fn is_zero(&self) -> bool {
        Fp::is_zero(self)
    }

    fn pow_usize(&self, exponent: usize) -> Self {
        self.pow([exponent as u64, 0, 0, 0])
    }
}

/// 支持 radix-2 FFT 的有限域。
///
/// `root_of_unity(order, params)` 必须返回精确阶为 `order` 的单位根，
/// 而不只是满足 omega^order = 1 的元素。
/// 支持 radix-2 FFT 的有限域。
///
/// root_of_unity(log_order) 必须返回精确阶为 2^log_order 的单位根。
pub trait Radix2Field: PolynomialField {
    /// 最大支持的 2^s 阶单位根。
    fn two_adicity(&self) -> u32;

    /// 返回指定阶数 2^log_order 的原始单位根。
    fn root_of_unity(&self, log_order: u32) -> Option<Self>;

    /// 乘法群生成元，用于 coset FFT。
    fn multiplicative_generator(&self) -> Self;
}

impl Radix2Field for Fp {
    fn two_adicity(&self) -> u32 {
        let p_minus_one = sub_one(self.params.modulus);
        trailing_zeros_256(p_minus_one)
    }

    fn root_of_unity(&self, log_order: u32) -> Option<Self> {
        let s = self.two_adicity();

        if log_order > s {
            return None;
        }

        if log_order == 0 {
            return Some(self.one_like());
        }

        let p_minus_one = sub_one(self.params.modulus);
        let exponent = shr_256(p_minus_one, log_order);

        let generator = self.multiplicative_generator();
        let omega = generator.pow(exponent);

        // 验证 omega^(2^log_order) = 1。
        let order = 1usize.checked_shl(log_order)?;
        if omega.pow_usize(order) != self.one_like() {
            return None;
        }

        // 验证 omega 的阶不是 2^(log_order - 1)。
        if omega.pow_usize(order / 2) == self.one_like() {
            return None;
        }

        Some(omega)
    }

    fn multiplicative_generator(&self) -> Self {
        const FR_MODULUS: [u64; 4] = [
            0x43e1_f593_f000_0001,
            0x2833_e848_79b9_7091,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ];

        const FQ_MODULUS: [u64; 4] = [
            0x3c208c16d87cfd47,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ];

        if self.params.modulus == FR_MODULUS {
            Self::from_u64_like(self, 5)
        } else if self.params.modulus == FQ_MODULUS {
            Self::from_u64_like(self, 3)
        } else {
            panic!("Radix2Field: unsupported Fp modulus");
        }
    }
}

/// p - 1
fn sub_one(mut value: [u64; 4]) -> [u64; 4] {
    for limb in &mut value {
        let (result, borrow) = limb.overflowing_sub(1);
        *limb = result;

        if !borrow {
            break;
        }
    }

    value
}

/// Count trailing zero bits in a 256-bit little-endian integer.
fn trailing_zeros_256(value: [u64; 4]) -> u32 {
    for (i, limb) in value.iter().enumerate() {
        if *limb != 0 {
            return (i as u32) * 64 + limb.trailing_zeros();
        }
    }

    256
}

/// Right shift a 256-bit little-endian integer.
fn shr_256(value: [u64; 4], bits: u32) -> [u64; 4] {
    if bits >= 256 {
        return [0; 4];
    }

    let words = (bits / 64) as usize;
    let rem = bits % 64;

    let mut result = [0u64; 4];

    for i in 0..4 {
        let src = i + words;

        if src >= 4 {
            break;
        }

        result[i] = value[src] >> rem;

        if rem != 0 && src + 1 < 4 {
            result[i] |= value[src + 1] << (64 - rem);
        }
    }

    result
}
