use std::fmt;

use super::field::PolynomialField;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PolynomialError {
    DivisionByZeroPolynomial,
    NonInvertibleLeadingCoefficient,
    InvalidTransposeSize {
        a_len: usize,
        n: usize,
        c_len: usize,
    },
}

impl fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DivisionByZeroPolynomial => {
                write!(f, "division by zero polynomial")
            }
            Self::NonInvertibleLeadingCoefficient => {
                write!(f, "leading coefficient is not invertible")
            }
            Self::InvalidTransposeSize { a_len, n, c_len } => {
                write!(
                    f,
                    "invalid transpose multiplication sizes: a_len={}, n={}, c_len={}",
                    a_len, n, c_len
                )
            }
        }
    }
}

impl std::error::Error for PolynomialError {}

/// 判断多项式是否为零。
///
/// 空向量也表示零多项式。
pub fn is_zero<F: PolynomialField>(a: &[F]) -> bool {
    a.iter().all(|x| x.is_zero())
}

/// 删除最高次项的连续零系数。
///
/// 零多项式统一表示为空 Vec。
pub fn condense<F: PolynomialField>(a: &mut Vec<F>) {
    while a.last().is_some_and(|x| x.is_zero()) {
        a.pop();
    }
}

/// 反转多项式系数，并调整为指定长度。
///
/// 与 C++ _reverse 对应：
/// 1. reverse
/// 2. resize(n)
///
/// 如果 n 大于原长度，在末尾补零。
pub fn reverse<F: PolynomialField>(a: &mut Vec<F>, n: usize, params: &F::Params) {
    a.reverse();
    a.resize_with(n, || F::zero(params));
}

/// 多项式加法。
pub fn polynomial_addition<F: PolynomialField>(a: &[F], b: &[F], params: &F::Params) -> Vec<F> {
    if is_zero(a) {
        let mut result = b.to_vec();
        condense(&mut result);
        return result;
    }

    if is_zero(b) {
        let mut result = a.to_vec();
        condense(&mut result);
        return result;
    }

    let size = a.len().max(b.len());
    let mut result = Vec::with_capacity(size);

    for i in 0..size {
        let value = match (a.get(i), b.get(i)) {
            (Some(x), Some(y)) => x.add(y),
            (Some(x), None) => x.clone(),
            (None, Some(y)) => y.clone(),
            (None, None) => F::zero(params),
        };

        result.push(value);
    }

    condense(&mut result);
    result
}

/// 多项式减法。
pub fn polynomial_subtraction<F: PolynomialField>(a: &[F], b: &[F], params: &F::Params) -> Vec<F> {
    if is_zero(b) {
        let mut result = a.to_vec();
        condense(&mut result);
        return result;
    }

    if is_zero(a) {
        let mut result: Vec<F> = b.iter().map(|x| F::zero(params).sub(x)).collect();

        condense(&mut result);
        return result;
    }

    let size = a.len().max(b.len());
    let mut result = Vec::with_capacity(size);

    for i in 0..size {
        let value = match (a.get(i), b.get(i)) {
            (Some(x), Some(y)) => x.sub(y),
            (Some(x), None) => x.clone(),
            (None, Some(y)) => F::zero(params).sub(y),
            (None, None) => F::zero(params),
        };

        result.push(value);
    }

    condense(&mut result);
    result
}

/// 多项式乘法：朴素 O(n*m) 算法。
///
/// 结果系数顺序仍然是低次到高次。
///
/// 这是正确性基线。后续接入 radix-2 FFT 后，
/// polynomial_multiplication 可以切换到 FFT 实现。
pub fn polynomial_multiplication<F: PolynomialField>(
    a: &[F],
    b: &[F],
    params: &F::Params,
) -> Vec<F> {
    if is_zero(a) || is_zero(b) {
        return Vec::new();
    }

    let result_len = a.len() + b.len() - 1;

    let mut result: Vec<F> = (0..result_len).map(|_| F::zero(params)).collect();

    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            let product = ai.mul(bj);
            result[i + j] = result[i + j].add(&product);
        }
    }

    condense(&mut result);
    result
}

/// 与 libfqfft 的转置多项式乘法对应。
///
/// 输入：
/// - n：要求的输出参数
/// - a：长度 m 的多项式
/// - c：另一个多项式
///
/// 输出：
/// product(reverse(a, m) * c)[m-1 .. m+n]
///
/// 返回 n+1 个系数。
pub fn polynomial_multiplication_transpose<F: PolynomialField>(
    n: usize,
    a: &[F],
    c: &[F],
    params: &F::Params,
) -> Result<Vec<F>, PolynomialError> {
    let m = a.len();

    if m == 0 {
        return Ok(Vec::new());
    }

    // C++ 原条件：
    // c.size() - 1 > m + n
    //
    // 对空 c 的情况单独处理，避免 size_t 下溢。
    if c.is_empty() || c.len() - 1 > m.saturating_add(n) {
        return Err(PolynomialError::InvalidTransposeSize {
            a_len: m,
            n,
            c_len: c.len(),
        });
    }

    let mut reversed_a = a.to_vec();
    reverse(&mut reversed_a, m, params);

    let product = polynomial_multiplication(&reversed_a, c, params);

    let start = m - 1;
    let end = m.saturating_add(n);

    let mut result = Vec::with_capacity(n.saturating_add(1));

    for i in start..end {
        result.push(product.get(i).cloned().unwrap_or_else(|| F::zero(params)));
    }

    Ok(result)
}

/// 多项式除法。
///
/// 返回 (quotient, remainder)，满足：
///```text
///     a = quotient * b + remainder
///```
/// 且 remainder 的次数小于 b 的次数。
///
/// q 和 r 都采用低次项在前的系数顺序。
pub fn polynomial_division<F: PolynomialField>(
    a: &[F],
    b: &[F],
    params: &F::Params,
) -> Result<(Vec<F>, Vec<F>), PolynomialError> {
    let mut dividend = a.to_vec();
    let mut divisor = b.to_vec();

    condense(&mut dividend);
    condense(&mut divisor);

    if divisor.is_empty() {
        return Err(PolynomialError::DivisionByZeroPolynomial);
    }

    if dividend.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    if dividend.len() < divisor.len() {
        return Ok((Vec::new(), dividend));
    }

    let leading_inverse = divisor
        .last()
        .expect("divisor is nonempty")
        .inverse()
        .ok_or(PolynomialError::NonInvertibleLeadingCoefficient)?;

    let quotient_len = dividend.len() - divisor.len() + 1;

    let mut quotient: Vec<F> = (0..quotient_len).map(|_| F::zero(params)).collect();

    let mut remainder = dividend;

    while !remainder.is_empty() && remainder.len() >= divisor.len() {
        let shift = remainder.len() - divisor.len();

        let lead = remainder
            .last()
            .expect("remainder is nonempty")
            .mul(&leading_inverse);

        quotient[shift] = quotient[shift].add(&lead);

        // remainder -= divisor * lead * x^shift
        for (i, coefficient) in divisor.iter().enumerate() {
            let product = coefficient.mul(&lead);
            let index = i + shift;

            remainder[index] = remainder[index].sub(&product);
        }

        condense(&mut remainder);
    }

    condense(&mut quotient);
    condense(&mut remainder);

    Ok((quotient, remainder))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct TestField(u64);

    const MODULUS: u64 = 97;

    impl PolynomialField for TestField {
        type Params = u64;

        fn zero(_: &u64) -> Self {
            Self(0)
        }

        fn one(_: &u64) -> Self {
            Self(1)
        }

        fn from_u64(value: u64, _: &u64) -> Self {
            Self(value % MODULUS)
        }

        fn add(&self, rhs: &Self) -> Self {
            Self((self.0 + rhs.0) % MODULUS)
        }

        fn sub(&self, rhs: &Self) -> Self {
            Self((self.0 + MODULUS - rhs.0) % MODULUS)
        }

        fn mul(&self, rhs: &Self) -> Self {
            Self((self.0 * rhs.0) % MODULUS)
        }
        fn one_like(&self) -> Self {
            Self(1)
        }
        fn from_u64_like(&self, value: u64) -> Self {
            TestField(value % 97)
        }
        fn inverse(&self) -> Option<Self> {
            if self.0 == 0 {
                return None;
            }

            let mut t: i64 = 0;
            let mut new_t: i64 = 1;
            let mut r: i64 = MODULUS as i64;
            let mut new_r: i64 = self.0 as i64;

            while new_r != 0 {
                let q = r / new_r;

                (t, new_t) = (new_t, t - q * new_t);
                (r, new_r) = (new_r, r - q * new_r);
            }

            if r != 1 {
                return None;
            }

            Some(Self(t.rem_euclid(MODULUS as i64) as u64))
        }

        fn is_zero(&self) -> bool {
            self.0 == 0
        }

        fn pow_usize(&self, exponent: usize) -> Self {
            let mut base = self.clone();
            let mut exp = exponent;
            let mut result = Self(1);

            while exp > 0 {
                if exp & 1 == 1 {
                    result = result.mul(&base);
                }

                base = base.mul(&base);
                exp >>= 1;
            }

            result
        }
    }

    fn f(x: u64) -> TestField {
        TestField(x % MODULUS)
    }

    fn values(a: &[TestField]) -> Vec<u64> {
        a.iter().map(|x| x.0).collect()
    }

    #[test]
    fn test_is_zero_and_condense() {
        assert!(is_zero::<TestField>(&[]));
        assert!(is_zero(&[f(0), f(0)]));

        let mut p = vec![f(1), f(2), f(0), f(0)];
        condense(&mut p);

        assert_eq!(values(&p), vec![1, 2]);
    }

    #[test]
    fn test_addition() {
        let a = vec![f(1), f(2)];
        let b = vec![f(3), f(4), f(5)];

        let result = polynomial_addition(&a, &b, &97);

        assert_eq!(values(&result), vec![4, 6, 5]);
    }

    #[test]
    fn test_subtraction() {
        let a = vec![f(1), f(2)];
        let b = vec![f(3), f(4), f(5)];

        let result = polynomial_subtraction(&a, &b, &97);

        assert_eq!(values(&result), vec![95, 95, 92]);
    }

    #[test]
    fn test_multiplication() {
        // (1 + 2x) * (3 + 4x)
        // = 3 + 10x + 8x^2

        let a = vec![f(1), f(2)];
        let b = vec![f(3), f(4)];

        let result = polynomial_multiplication(&a, &b, &97);

        assert_eq!(values(&result), vec![3, 10, 8]);
    }

    #[test]
    fn test_transpose_multiplication() {
        // a = [1, 2], m = 2
        // reverse(a, 2) = [2, 1]
        // c = [3, 4, 5]
        //
        // product = [6, 11, 14, 5]
        //
        // n = 1
        // result = product[1..3] = [11, 14]

        let a = vec![f(1), f(2)];
        let c = vec![f(3), f(4), f(5)];

        let result = polynomial_multiplication_transpose(1, &a, &c, &97).unwrap();

        assert_eq!(values(&result), vec![11, 14]);
    }

    #[test]
    fn test_division_exact() {
        // (1 + 3x + 2x²) / (1 + x)
        // = 1 + 2x

        let a = vec![f(1), f(3), f(2)];
        let b = vec![f(1), f(1)];

        let (q, r) = polynomial_division(&a, &b, &97).unwrap();

        assert_eq!(values(&q), vec![1, 2]);
        assert!(r.is_empty());
    }

    #[test]
    fn test_division_with_remainder() {
        // (1 + x²) / (1 + x)
        // = (x - 1) * (x + 1) + 2
        //
        // Mod 97: x - 1 = 96 + x

        let a = vec![f(1), f(0), f(1)];
        let b = vec![f(1), f(1)];

        let (q, r) = polynomial_division(&a, &b, &97).unwrap();

        assert_eq!(values(&q), vec![96, 1]);
        assert_eq!(values(&r), vec![2]);
    }

    #[test]
    fn test_division_by_zero() {
        let a = vec![f(1), f(2)];
        let b: Vec<TestField> = vec![];

        assert_eq!(
            polynomial_division(&a, &b, &97),
            Err(PolynomialError::DivisionByZeroPolynomial),
        );
    }

    #[test]
    fn test_zero_polynomial_multiplication() {
        let a = vec![f(1), f(2)];
        let b = vec![f(0), f(0)];

        let result = polynomial_multiplication(&a, &b, &97);

        assert!(result.is_empty());
    }
}
