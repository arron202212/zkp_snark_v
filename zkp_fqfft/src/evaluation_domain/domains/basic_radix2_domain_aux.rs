use crate::polynomial_arithmetic::field::PolynomialField;

/// FFT 错误。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FftError {
    EmptyInput,
    SizeNotPowerOfTwo,
    InvalidFieldSize,
    ZeroInverse,
    SizeOverflow,
}

/// 计算 log2(n)。
fn log2_exact(n: usize) -> Result<u32, FftError> {
    if n == 0 {
        return Err(FftError::EmptyInput);
    }

    if !n.is_power_of_two() {
        return Err(FftError::SizeNotPowerOfTwo);
    }

    Ok(n.trailing_zeros())
}

/// 计算 bit-reversal。
fn bit_reverse(value: usize, bits: u32) -> usize {
    value.reverse_bits() >> (usize::BITS - bits)
}

/// 计算 radix-2 FFT。
///
/// 输入：
///```text           
///     a[i] = 多项式的第 i 个系数
///
/// 输出：
///
///     a[i] = f(omega^i)
///```
/// 注意：这里不进行 1/n 缩放。
pub fn basic_radix2_fft<F: PolynomialField>(a: &mut [F], omega: &F) -> Result<(), FftError> {
    let n = a.len();
    let log_n = log2_exact(n)?;

    // Bit-reversal permutation.
    for k in 0..n {
        let rk = bit_reverse(k, log_n);

        if k < rk {
            a.swap(k, rk);
        }
    }

    // Iterative Cooley-Tukey FFT.
    let mut m = 1usize;

    for _s in 1..=log_n {
        let two_m = m.checked_mul(2).ok_or(FftError::SizeOverflow)?;

        // w_m = omega^(n / (2*m))
        let exponent = n / two_m;
        let w_m = omega.pow_usize(exponent);

        let mut k = 0usize;

        while k < n {
            let mut w = omega.one_like();

            for j in 0..m {
                let t = w.mul(&a[k + j + m]);

                let left = a[k + j].clone();

                a[k + j + m] = left.sub(&t);
                a[k + j] = left.add(&t);

                w = w.mul(&w_m);
            }

            k += two_m;
        }

        m = two_m;
    }

    Ok(())
}

/// 多线程 FFT 接口。
///
/// 当前先使用串行算法，保持与 C++ 串行实现相同的数学结果。
/// 后续可以在不改变调用接口的情况下替换为 Rayon 实现。
pub fn parallel_basic_radix2_fft<F: PolynomialField>(
    a: &mut [F],
    omega: &F,
) -> Result<(), FftError> {
    basic_radix2_fft(a, omega)
}

/// 将多项式系数转换到 coset。
///
/// a[i] = a[i] * g^i
///
/// 对应 C++：
/// _multiply_by_coset(a, g)
pub fn multiply_by_coset<F: PolynomialField>(a: &mut [F], g: &F) {
    let mut u = g.clone();

    for value in a.iter_mut().skip(1) {
        *value = value.mul(&u);
        u = u.mul(g);
    }
}

/// 计算全部 Lagrange 基函数在 t 处的值。
///```text
/// S = {1, omega, omega^2, ..., omega^(m-1)}
///
/// 返回：
///
///     [L_0(t), L_1(t), ..., L_(m-1)(t)]
///```
/// 要求 omega 是 m 阶原始单位根。
pub fn basic_radix2_evaluate_all_lagrange_polynomials<F: PolynomialField>(
    m: usize,
    t: &F,
    omega: &F,
    params: &F::Params,
) -> Result<Vec<F>, FftError> {
    if m == 0 {
        return Err(FftError::EmptyInput);
    }

    if m == 1 {
        return Ok(vec![F::one(params)]);
    }

    if !m.is_power_of_two() {
        return Err(FftError::SizeNotPowerOfTwo);
    }

    let m_u64 = u64::try_from(m).map_err(|_| FftError::SizeOverflow)?;

    let one = F::one(params);

    let mut result = vec![F::zero(params); m];

    // Z(t) = t^m - 1
    let z = t.pow_usize(m).sub(&one);

    // If t is one of the roots of unity:
    //
    // L_i(t) = 1 if t = omega^i
    // L_i(t) = 0 otherwise
    if z.is_zero() {
        let mut omega_i = one.clone();

        for value in result.iter_mut() {
            if omega_i == *t {
                *value = one;
                return Ok(result);
            }

            omega_i = omega_i.mul(omega);
        }

        // If omega is a valid m-th root, this branch should not occur.
        return Err(FftError::InvalidFieldSize);
    }

    // l = Z(t) / m
    let m_element = F::from_u64(m_u64, params);

    let m_inverse = m_element.inverse().ok_or(FftError::ZeroInverse)?;

    let mut l = z.mul(&m_inverse);
    let mut r = one;

    // C++:
    //
    // u[i] = l * (t - r).inverse()
    // l *= omega
    // r *= omega
    //
    for value in result.iter_mut() {
        let denominator = t.sub(&r);

        let denominator_inverse = denominator.inverse().ok_or(FftError::ZeroInverse)?;

        *value = l.mul(&denominator_inverse);

        l = l.mul(omega);
        r = r.mul(omega);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    use zkp_ff::algebra::fields::prime_base::fp::Fp;
    use zkp_ff::algebra::fields::prime_base::fp_montgomery::FieldParams;

    fn fr_params() -> FieldParams {
        FieldParams::new([
            0x43e1_f593_f000_0001,
            0x2833_e848_79b9_7091,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ])
    }

    fn canonical(value: &Fp) -> [u64; 4] {
        value.to_canonical()
    }

    #[test]
    fn test_bit_reverse() {
        assert_eq!(bit_reverse(0, 3), 0);
        assert_eq!(bit_reverse(1, 3), 4);
        assert_eq!(bit_reverse(2, 3), 2);
        assert_eq!(bit_reverse(3, 3), 6);
        assert_eq!(bit_reverse(4, 3), 1);
    }

    #[test]
    fn test_coset_multiplication() {
        let params = fr_params();

        let g = Fp::from_u64(3, params.clone());

        let mut values = vec![
            Fp::from_u64(2, params.clone()),
            Fp::from_u64(3, params.clone()),
            Fp::from_u64(4, params.clone()),
        ];

        multiply_by_coset(&mut values, &g);

        assert_eq!(canonical(&values[0]), [2, 0, 0, 0]);

        assert_eq!(canonical(&values[1]), [9, 0, 0, 0]);

        assert_eq!(canonical(&values[2]), [36, 0, 0, 0]);
    }

    #[test]
    fn test_fft_size_validation() {
        let params = fr_params();

        let omega = Fp::one(params.clone());

        let mut values = vec![
            Fp::from_u64(1, params.clone()),
            Fp::from_u64(2, params.clone()),
            Fp::from_u64(3, params.clone()),
        ];

        assert_eq!(
            basic_radix2_fft(&mut values, &omega),
            Err(FftError::SizeNotPowerOfTwo)
        );
    }

    #[test]
    fn test_lagrange_m_equals_one() {
        let params = fr_params();

        let t = Fp::from_u64(123, params.clone());
        let omega = Fp::one(params.clone());

        let result =
            basic_radix2_evaluate_all_lagrange_polynomials(1, &t, &omega, &params).unwrap();

        assert_eq!(result.len(), 1);
        assert_eq!(result[0], Fp::one(params));
    }
}
