use core::cmp::Ordering;

/// 兼容旧代码：BN254 默认使用 4 limbs。
pub type Limbs = [u64; 4];

/// 泛型 limb 数组，用于其他曲线。
pub type LimbsN<const N: usize> = [u64; N];

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FieldParams<const N: usize = 4> {
    pub modulus: LimbsN<N>,
    pub n0_inv: u64,
    pub r: LimbsN<N>,
    pub r2: LimbsN<N>,
}

impl<const N: usize> FieldParams<N> {
    /// 创建 Montgomery 参数。
    ///
    /// R  = 2^(64*N)
    /// R2 = R^2 mod modulus
    ///
    /// n0_inv = -modulus[0]^(-1) mod 2^64
    pub fn new(modulus: LimbsN<N>) -> Self {
        assert!(N > 0, "limb count must be greater than zero");
        assert!(!is_zero(&modulus), "modulus must not be zero");
        assert!(modulus[0] & 1 == 1, "modulus must be odd");
        assert!(!is_one(&modulus), "modulus must be greater than one");

        let n0_inv = compute_n0_inv(modulus[0]);

        // R = 2^(64*N) mod modulus
        let mut r = [0u64; N];
        r[0] = 1;

        for _ in 0..(64 * N) {
            r = add_mod(r, r, modulus);
        }

        // R2 = R * 2^(64*N) mod modulus
        let mut r2 = r;

        for _ in 0..(64 * N) {
            r2 = add_mod(r2, r2, modulus);
        }

        Self {
            modulus,
            n0_inv,
            r,
            r2,
        }
    }

    pub fn validate(&self) -> bool {
        N > 0
            && !is_zero(&self.modulus)
            && !is_one(&self.modulus)
            && self.modulus[0] & 1 == 1
            && self.modulus[0].wrapping_mul(self.n0_inv) == u64::MAX
    }
}

/// 计算 -a^(-1) mod 2^64。
///
/// a 必须为奇数。
fn compute_n0_inv(a: u64) -> u64 {
    assert!(a & 1 == 1);

    let mut x = 1u64;

    // Newton iteration，每次迭代有效位数翻倍。
    for _ in 0..6 {
        x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
    }

    x.wrapping_neg()
}

// -----------------------------------------------------------------------------
// Basic limb operations
// -----------------------------------------------------------------------------

pub fn is_zero<const N: usize>(a: &LimbsN<N>) -> bool {
    a.iter().all(|&x| x == 0)
}

pub fn is_one<const N: usize>(a: &LimbsN<N>) -> bool {
    N > 0 && a[0] == 1 && a[1..].iter().all(|&x| x == 0)
}

pub fn cmp<const N: usize>(a: &LimbsN<N>, b: &LimbsN<N>) -> Ordering {
    for i in (0..N).rev() {
        match a[i].cmp(&b[i]) {
            Ordering::Equal => {}
            other => return other,
        }
    }

    Ordering::Equal
}

/// N-limb 加法，返回结果和最高进位。
pub fn add_raw<const N: usize>(a: LimbsN<N>, b: LimbsN<N>) -> (LimbsN<N>, u64) {
    let mut out = [0u64; N];
    let mut carry = 0u128;

    for i in 0..N {
        let sum = a[i] as u128 + b[i] as u128 + carry;

        out[i] = sum as u64;
        carry = sum >> 64;
    }

    (out, carry as u64)
}

/// N-limb 减法，返回结果和最高借位。
pub fn sub_raw<const N: usize>(a: LimbsN<N>, b: LimbsN<N>) -> (LimbsN<N>, u64) {
    let mut out = [0u64; N];
    let mut borrow = 0u128;

    for i in 0..N {
        let ai = a[i] as u128;
        let bi = b[i] as u128 + borrow;

        if ai >= bi {
            out[i] = (ai - bi) as u64;
            borrow = 0;
        } else {
            out[i] = ((1u128 << 64) + ai - bi) as u64;
            borrow = 1;
        }
    }

    (out, borrow as u64)
}

/// 模加法。
///
/// 前提：a < modulus，b < modulus。
pub fn add_mod<const N: usize>(a: LimbsN<N>, b: LimbsN<N>, modulus: LimbsN<N>) -> LimbsN<N> {
    let (mut sum, carry) = add_raw(a, b);

    if carry != 0 || cmp(&sum, &modulus) != Ordering::Less {
        // 若 carry=1，则 sum < modulus，且扩展减法的借位
        // 会抵消最高位的 carry。
        sum = sub_raw(sum, modulus).0;
    }

    sum
}

/// 模减法。
///
/// 前提：a < modulus，b < modulus。
pub fn sub_mod<const N: usize>(a: LimbsN<N>, b: LimbsN<N>, modulus: LimbsN<N>) -> LimbsN<N> {
    let (mut diff, borrow) = sub_raw(a, b);

    if borrow != 0 {
        // 结果按 2^(64*N) 截断，等价于 a - b + modulus。
        diff = add_raw(diff, modulus).0;
    }

    diff
}

// -----------------------------------------------------------------------------
// Montgomery multiplication
// -----------------------------------------------------------------------------

/// Montgomery multiplication。
///
/// 计算：`result = a * b * R^(-1) mod p`。
///
/// 其中：`R = 2^(64*N)`。
///
/// 前提：`a < p` 且 `b < p`。
///
/// 使用 CIOS Montgomery multiplication。
///
/// 注意：这里采用值传递，兼容旧版 fp.rs 的调用方式。
pub fn montgomery_mul<const N: usize>(
    a: LimbsN<N>,
    b: LimbsN<N>,
    params: &FieldParams<N>,
) -> LimbsN<N> {
    assert!(N > 0);

    let modulus = &params.modulus;
    let n0_inv = params.n0_inv;

    debug_assert!(cmp(&a, modulus) == Ordering::Less);
    debug_assert!(cmp(&b, modulus) == Ordering::Less);

    // N + 2 limbs 的临时工作区。
    //
    // 使用 Vec，避免稳定 Rust 对 [u64; N + 2] 的 const
    // generic 表达式限制。
    let mut t = vec![0u64; N + 2];

    for i in 0..N {
        // ---------------------------------------------------------
        // t += a * b[i]
        // ---------------------------------------------------------
        let mut carry = 0u128;

        for j in 0..N {
            let uv = (a[j] as u128) * (b[i] as u128) + (t[j] as u128) + carry;

            t[j] = uv as u64;
            carry = uv >> 64;
        }

        let uv = t[N] as u128 + carry;
        t[N] = uv as u64;

        let uv2 = t[N + 1] as u128 + (uv >> 64);
        t[N + 1] = uv2 as u64;

        debug_assert_eq!(uv2 >> 64, 0);

        // ---------------------------------------------------------
        // m = t[0] * n0_inv mod 2^64
        // ---------------------------------------------------------
        let m = t[0].wrapping_mul(n0_inv);

        // ---------------------------------------------------------
        // t += m * modulus
        // ---------------------------------------------------------
        carry = 0;

        for j in 0..N {
            let uv = (m as u128) * (modulus[j] as u128) + (t[j] as u128) + carry;

            t[j] = uv as u64;
            carry = uv >> 64;
        }

        let uv = t[N] as u128 + carry;
        t[N] = uv as u64;

        let uv2 = t[N + 1] as u128 + (uv >> 64);
        t[N + 1] = uv2 as u64;

        debug_assert_eq!(uv2 >> 64, 0);

        // t[0] 应当为 0，可以除以 2^64。
        debug_assert_eq!(t[0], 0);

        // ---------------------------------------------------------
        // t >>= 64
        // ---------------------------------------------------------
        for j in 0..(N + 1) {
            t[j] = t[j + 1];
        }

        t[N + 1] = 0;
    }

    let mut result = [0u64; N];
    result.copy_from_slice(&t[..N]);

    // CIOS 最终结果可能占 N+1 limbs。
    //
    // 因为 a,b < p，所以 Montgomery 结果小于 2p。
    // 因此最多需要减去一次 modulus。
    debug_assert!(t[N] <= 1);

    if t[N] != 0 || cmp(&result, modulus) != Ordering::Less {
        let (reduced, borrow) = sub_raw(result, *modulus);

        // 若高位 t[N]=1，低 N limbs 必须借位；
        // 若高位 t[N]=0，则低 N limbs 不应借位。
        debug_assert_eq!(borrow, t[N]);

        result = reduced;
    }

    result
}

/// Montgomery 平方。
pub fn montgomery_square<const N: usize>(a: &LimbsN<N>, params: &FieldParams<N>) -> LimbsN<N> {
    montgomery_mul(*a, *a, params)
}

// -----------------------------------------------------------------------------
// Conversion helpers
// -----------------------------------------------------------------------------

/// 将普通整数 x 转换为 Montgomery 表示。
///
/// `mont(x) = x * R mod p`。
///
/// 通过 montgomery_mul(x, R²) 实现。
///
/// 要求 x < p。
pub fn to_montgomery<const N: usize>(x: &LimbsN<N>, params: &FieldParams<N>) -> LimbsN<N> {
    montgomery_mul(*x, params.r2, params)
}

/// 将 Montgomery 表示转换回普通整数。
///
/// `x = mont(x) * R^(-1) mod p`。
pub fn from_montgomery<const N: usize>(x: &LimbsN<N>, params: &FieldParams<N>) -> LimbsN<N> {
    let mut one = [0u64; N];
    one[0] = 1;

    montgomery_mul(*x, one, params)
}

// -----------------------------------------------------------------------------
// Small integer helpers
// -----------------------------------------------------------------------------

/// 普通整数加一，忽略溢出。
pub fn add_one<const N: usize>(a: &LimbsN<N>) -> LimbsN<N> {
    let mut out = *a;
    let mut carry = 1u64;

    for limb in out.iter_mut() {
        let (value, overflow) = limb.overflowing_add(carry);
        *limb = value;

        if !overflow {
            break;
        }

        carry = 1;
    }

    out
}

/// 普通整数减一，忽略下溢。
pub fn sub_one<const N: usize>(a: &LimbsN<N>) -> LimbsN<N> {
    let mut out = *a;
    let mut borrow = 1u64;

    for limb in out.iter_mut() {
        let (value, underflow) = limb.overflowing_sub(borrow);
        *limb = value;

        if !underflow {
            break;
        }

        borrow = 1;
    }

    out
}

/// 右移一位。
pub fn shr1<const N: usize>(a: &LimbsN<N>) -> LimbsN<N> {
    let mut out = [0u64; N];
    let mut carry = 0u64;

    for i in (0..N).rev() {
        out[i] = (a[i] >> 1) | carry;
        carry = a[i] << 63;
    }

    out
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn bn254_params() -> FieldParams<4> {
        // BN254 Fr modulus
        let modulus = [
            0x43e1_f593_f000_0001,
            0x2833_e848_79b9_7091,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ];

        FieldParams::new(modulus)
    }

    fn small_params() -> FieldParams<1> {
        // p = 17
        let modulus = [17u64];

        FieldParams::new(modulus)
    }

    fn small_params_5_limbs() -> FieldParams<5> {
        // 使用 5 limbs 表示 p = 17，
        // 用于验证泛型 Montgomery 运算。
        let modulus = [17u64, 0, 0, 0, 0];

        FieldParams::new(modulus)
    }

    #[test]
    fn test_add_sub_raw() {
        let a = [u64::MAX, 0, 0, 0];
        let b = [1, 0, 0, 0];

        let (sum, carry) = add_raw(a, b);

        assert_eq!(sum, [0, 1, 0, 0]);
        assert_eq!(carry, 0);

        let (diff, borrow) = sub_raw(sum, b);

        assert_eq!(diff, a);
        assert_eq!(borrow, 0);
    }

    #[test]
    fn test_add_raw_overflow() {
        let a = [u64::MAX, u64::MAX];
        let b = [1, 0];

        let (sum, carry) = add_raw(a, b);

        assert_eq!(sum, [0, 0]);
        assert_eq!(carry, 1);
    }

    #[test]
    fn test_sub_raw_underflow() {
        let a = [0u64, 0];
        let b = [1u64, 0];

        let (diff, borrow) = sub_raw(a, b);

        assert_eq!(diff, [u64::MAX, u64::MAX]);
        assert_eq!(borrow, 1);
    }

    #[test]
    fn test_bn254_params() {
        let params = bn254_params();

        assert!(params.validate());

        assert_eq!(params.modulus[0].wrapping_mul(params.n0_inv), u64::MAX);
    }

    #[test]
    fn test_small_field_roundtrip() {
        let params = small_params();

        let x = [5u64];

        let xm = to_montgomery(&x, &params);
        let recovered = from_montgomery(&xm, &params);

        assert_eq!(recovered, x);
    }

    #[test]
    fn test_small_field_mul() {
        let params = small_params();

        let a = to_montgomery(&[7], &params);
        let b = to_montgomery(&[9], &params);

        let c = montgomery_mul(a, b, &params);
        let result = from_montgomery(&c, &params);

        // 7 * 9 mod 17 = 12
        assert_eq!(result, [12]);
    }

    #[test]
    fn test_bn254_roundtrip() {
        let params = bn254_params();

        let x = [123456789u64, 0, 0, 0];

        let xm = to_montgomery(&x, &params);
        let recovered = from_montgomery(&xm, &params);

        assert_eq!(recovered, x);
    }

    #[test]
    fn test_bn254_mul() {
        let params = bn254_params();

        let x = [123u64, 0, 0, 0];
        let y = [456u64, 0, 0, 0];

        let xm = to_montgomery(&x, &params);
        let ym = to_montgomery(&y, &params);

        let zm = montgomery_mul(xm, ym, &params);
        let z = from_montgomery(&zm, &params);

        assert_eq!(z, [56088, 0, 0, 0]);
    }

    #[test]
    fn test_5_limb_roundtrip() {
        let params = small_params_5_limbs();

        let x = [5u64, 0, 0, 0, 0];

        let xm = to_montgomery(&x, &params);
        let recovered = from_montgomery(&xm, &params);

        assert_eq!(recovered, x);
    }

    #[test]
    fn test_5_limb_mul() {
        let params = small_params_5_limbs();

        let x = [7u64, 0, 0, 0, 0];
        let y = [9u64, 0, 0, 0, 0];

        let xm = to_montgomery(&x, &params);
        let ym = to_montgomery(&y, &params);

        let zm = montgomery_mul(xm, ym, &params);
        let z = from_montgomery(&zm, &params);

        assert_eq!(z, [12, 0, 0, 0, 0]);
    }
}
