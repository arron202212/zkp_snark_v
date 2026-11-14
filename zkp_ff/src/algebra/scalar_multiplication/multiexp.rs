use std::cmp::Ordering;

use crate::algebra::field_utils::bigint::BigInt;

/// MSM 标量接口。
///
/// 约定：
/// - bit(0) 是最低有效位
/// - bit_len() 返回有效位数
/// - sub_assign(rhs) 要求 self >= rhs
pub trait MsmScalar: Clone {
    fn bit_len(&self) -> usize;

    fn bit(&self, index: usize) -> bool;

    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    fn cmp_scalar(&self, rhs: &Self) -> Ordering;

    /// self -= rhs。
    ///
    /// 调用前必须保证 self >= rhs。
    fn sub_assign(&mut self, rhs: &Self);
}

/// BigInt 实现 MSM 标量接口。
impl<const N: usize> MsmScalar for BigInt<N> {
    fn bit_len(&self) -> usize {
        self.num_bits()
    }

    fn bit(&self, index: usize) -> bool {
        self.test_bit(index)
    }

    fn is_zero(&self) -> bool {
        BigInt::is_zero(self)
    }

    fn is_one(&self) -> bool {
        self.num_bits() == 1 && self.data[0] == 1
    }

    fn cmp_scalar(&self, rhs: &Self) -> Ordering {
        BigInt::cmp(self, rhs)
    }

    fn sub_assign(&mut self, rhs: &Self) {
        assert!(
            BigInt::cmp(&*self, rhs) != Ordering::Less,
            "MsmScalar subtraction underflow"
        );

        let mut borrow = false;

        for i in 0..N {
            let (v1, b1) = self.data[i].overflowing_sub(rhs.data[i]);
            let (v2, b2) = v1.overflowing_sub(borrow as u64);

            self.data[i] = v2;
            borrow = b1 || b2;
        }

        debug_assert!(!borrow);
    }
}

/// 群元素接口。
pub trait MsmGroup<S: MsmScalar>: Clone {
    fn zero() -> Self;

    fn add(&self, rhs: &Self) -> Self;

    fn double(&self) -> Self {
        self.add(self)
    }

    fn mul_scalar(&self, scalar: &S) -> Self;

    fn mul_wnaf(&self, scalar: &S, _scalar_bits: usize) -> Self {
        self.mul_scalar(scalar)
    }

    fn is_zero(&self) -> bool;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiExpMethod {
    Naive,
    NaivePlain,
    Bdlo12,
    BosCoster,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MultiExpError {
    LengthMismatch { bases: usize, scalars: usize },
    InvalidChunkCount,
    WindowTooLarge,
    InvalidWindow,
}

fn check_lengths<G, S>(bases: &[G], scalars: &[S]) -> Result<(), MultiExpError> {
    if bases.len() != scalars.len() {
        return Err(MultiExpError::LengthMismatch {
            bases: bases.len(),
            scalars: scalars.len(),
        });
    }

    Ok(())
}

/// Naive：每一项使用 WNAF 标量乘法，然后累加。
pub fn multi_exp_naive<G, S>(bases: &[G], scalars: &[S]) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    check_lengths(bases, scalars)?;

    let mut result = G::zero();

    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        if scalar.is_zero() {
            continue;
        }

        result = result.add(&base.mul_wnaf(scalar, scalar.bit_len()));
    }

    Ok(result)
}

/// NaivePlain：普通标量乘法，然后累加。
pub fn multi_exp_naive_plain<G, S>(bases: &[G], scalars: &[S]) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    check_lengths(bases, scalars)?;

    let mut result = G::zero();

    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        if scalar.is_zero() {
            continue;
        }

        result = result.add(&base.mul_scalar(scalar));
    }

    Ok(result)
}

/// 根据基点数量选择 BDLO12/Pippenger 窗口大小。
fn bdlo12_window_size(length: usize) -> usize {
    let log2_length = if length <= 1 {
        0
    } else {
        usize::BITS as usize - 1 - length.leading_zeros() as usize
    };

    log2_length - (log2_length / 3).saturating_sub(2)
}

/// BDLO12/Pippenger 多标量乘法。
pub fn multi_exp_bdlo12<G, S>(
    bases: &[G],
    scalars: &[S],
    max_window: usize,
) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    check_lengths(bases, scalars)?;

    if bases.is_empty() {
        return Ok(G::zero());
    }

    let num_bits = scalars.iter().map(|s| s.bit_len()).max().unwrap_or(0);

    if num_bits == 0 {
        return Ok(G::zero());
    }

    let c = bdlo12_window_size(bases.len());

    if c == 0 || c >= usize::BITS as usize || c > max_window {
        return Err(MultiExpError::WindowTooLarge);
    }

    let bucket_count = 1usize
        .checked_shl(c as u32)
        .ok_or(MultiExpError::WindowTooLarge)?;

    let num_groups = (num_bits + c - 1) / c;

    let mut result: Option<G> = None;

    // 从最高窗口向最低窗口处理。
    for group in (0..num_groups).rev() {
        if let Some(current) = result.take() {
            let mut shifted = current;

            for _ in 0..c {
                shifted = shifted.double();
            }

            result = Some(shifted);
        }

        let mut buckets: Vec<Option<G>> = (0..bucket_count).map(|_| None).collect();

        // 将基点放入对应的桶。
        for (base, scalar) in bases.iter().zip(scalars.iter()) {
            let mut id = 0usize;

            for j in 0..c {
                if scalar.bit(group * c + j) {
                    id |= 1usize << j;
                }
            }

            if id == 0 {
                continue;
            }

            buckets[id] = Some(match buckets[id].take() {
                Some(old) => old.add(base),
                None => base.clone(),
            });
        }

        // 从高编号桶向低编号桶计算 running sum。
        let mut running_sum: Option<G> = None;

        for i in (1..bucket_count).rev() {
            if let Some(bucket) = buckets[i].take() {
                running_sum = Some(match running_sum.take() {
                    Some(old) => old.add(&bucket),
                    None => bucket,
                });
            }

            if let Some(running) = running_sum.as_ref() {
                result = Some(match result.take() {
                    Some(old) => old.add(running),
                    None => running.clone(),
                });
            }
        }
    }

    Ok(result.unwrap_or_else(G::zero))
}

/// Bos-Coster 多标量乘法。
///
/// a*A + b*B = (a-b)*A + b*(A+B)
pub fn multi_exp_bos_coster<G, S>(bases: &[G], scalars: &[S]) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    check_lengths(bases, scalars)?;

    if bases.is_empty() {
        return Ok(G::zero());
    }

    if bases.len() == 1 {
        if scalars[0].is_zero() {
            return Ok(G::zero());
        }

        return Ok(bases[0].mul_scalar(&scalars[0]));
    }

    let mut g = bases.to_vec();
    let mut exponents = scalars.to_vec();

    let mut result = G::zero();

    loop {
        let mut max_idx: Option<usize> = None;
        let mut second_idx: Option<usize> = None;

        // 找到最大和次大指数。
        for i in 0..exponents.len() {
            if exponents[i].is_zero() {
                continue;
            }

            match max_idx {
                None => {
                    max_idx = Some(i);
                }
                Some(m) => {
                    if exponents[i].cmp_scalar(&exponents[m]) == Ordering::Greater {
                        second_idx = max_idx;
                        max_idx = Some(i);
                    } else {
                        match second_idx {
                            None => {
                                second_idx = Some(i);
                            }
                            Some(s) => {
                                if exponents[i].cmp_scalar(&exponents[s]) == Ordering::Greater {
                                    second_idx = Some(i);
                                }
                            }
                        }
                    }
                }
            }
        }

        let a_idx = match max_idx {
            Some(i) => i,
            None => break,
        };

        let b_idx = match second_idx {
            Some(i) => i,
            None => {
                result =
                    result.add(&g[a_idx].mul_wnaf(&exponents[a_idx], exponents[a_idx].bit_len()));

                break;
            }
        };

        let a_bits = exponents[a_idx].bit_len();
        let b_bits = exponents[b_idx].bit_len();

        if b_bits > a_bits {
            return Err(MultiExpError::LengthMismatch {
                bases: a_bits,
                scalars: b_bits,
            });
        }

        // 保留原有的直接计算启发式。
        let bit_diff = a_bits - b_bits;
        let limit = bit_diff.min(20);
        let direct_threshold = 1usize << limit;

        if b_bits < direct_threshold {
            result = result.add(&g[a_idx].mul_wnaf(&exponents[a_idx], a_bits));

            // 使用临时变量避免可变和不可变借用冲突。
            let old = exponents[a_idx].clone();
            exponents[a_idx].sub_assign(&old);
        } else {
            // (a-b)*A + b*(A+B)
            let b_exp = exponents[b_idx].clone();

            exponents[a_idx].sub_assign(&b_exp);
            g[b_idx] = g[b_idx].add(&g[a_idx]);
        }
    }

    Ok(result)
}

/// 统一 MSM 入口。
pub fn multi_exp<G, S>(
    bases: &[G],
    scalars: &[S],
    method: MultiExpMethod,
    chunks: usize,
) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    check_lengths(bases, scalars)?;

    if chunks == 0 {
        return Err(MultiExpError::InvalidChunkCount);
    }

    let total = bases.len();

    if total == 0 {
        return Ok(G::zero());
    }

    if total < chunks || chunks == 1 {
        return multi_exp_inner(bases, scalars, method);
    }

    let one = total / chunks;
    let mut final_result = G::zero();

    for i in 0..chunks {
        let start = i * one;

        let end = if i == chunks - 1 {
            total
        } else {
            (i + 1) * one
        };

        let partial = multi_exp_inner(&bases[start..end], &scalars[start..end], method)?;

        final_result = final_result.add(&partial);
    }

    Ok(final_result)
}

fn multi_exp_inner<G, S>(
    bases: &[G],
    scalars: &[S],
    method: MultiExpMethod,
) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    match method {
        MultiExpMethod::Naive => multi_exp_naive(bases, scalars),
        MultiExpMethod::NaivePlain => multi_exp_naive_plain(bases, scalars),
        MultiExpMethod::Bdlo12 => multi_exp_bdlo12(bases, scalars, 16),
        MultiExpMethod::BosCoster => multi_exp_bos_coster(bases, scalars),
    }
}

/// mixed addition 版本。
///
/// 标量为 0：跳过。
/// 标量为 1：直接累加基点。
/// 其他标量：交给 MSM 算法。
pub fn multi_exp_with_mixed_addition<G, S>(
    bases: &[G],
    scalars: &[S],
    method: MultiExpMethod,
    chunks: usize,
) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    check_lengths(bases, scalars)?;

    let mut acc = G::zero();

    let mut remaining_bases = Vec::new();
    let mut remaining_scalars = Vec::new();

    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        if scalar.is_zero() {
            continue;
        }

        if scalar.is_one() {
            acc = acc.add(base);
        } else {
            remaining_bases.push(base.clone());
            remaining_scalars.push(scalar.clone());
        }
    }

    let rest = multi_exp(&remaining_bases, &remaining_scalars, method, chunks)?;

    Ok(acc.add(&rest))
}

/// 固定基点窗口表。
///
/// table[outer][inner] = inner * (2^(outer*window) * base)
pub type WindowTable<G> = Vec<Vec<G>>;

/// 构建固定窗口表。
pub fn get_window_table<G, S>(
    scalar_size: usize,
    window: usize,
    base: &G,
) -> Result<WindowTable<G>, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    let _ = std::marker::PhantomData::<S>;

    if window == 0 || window >= usize::BITS as usize {
        return Err(MultiExpError::InvalidWindow);
    }

    if scalar_size == 0 {
        return Ok(Vec::new());
    }

    let in_window = 1usize
        .checked_shl(window as u32)
        .ok_or(MultiExpError::WindowTooLarge)?;

    let outer_count = (scalar_size + window - 1) / window;

    let last_bits = scalar_size - (outer_count - 1) * window;
    let last_count = 1usize << last_bits;

    let mut table = Vec::with_capacity(outer_count);
    let mut g_outer = base.clone();

    for outer in 0..outer_count {
        let mut row = Vec::with_capacity(in_window);
        let mut g_inner = G::zero();

        let current_count = if outer == outer_count - 1 {
            last_count
        } else {
            in_window
        };

        for _ in 0..current_count {
            row.push(g_inner.clone());
            g_inner = g_inner.add(&g_outer);
        }

        table.push(row);

        for _ in 0..window {
            g_outer = g_outer.double();
        }
    }

    Ok(table)
}

/// 固定窗口标量乘法。
pub fn windowed_exp<G, S>(
    scalar_size: usize,
    window: usize,
    table: &WindowTable<G>,
    exponent: &S,
) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    if window == 0 || window >= usize::BITS as usize {
        return Err(MultiExpError::InvalidWindow);
    }

    if scalar_size == 0 {
        return Ok(G::zero());
    }

    let outer_count = (scalar_size + window - 1) / window;

    if table.len() < outer_count {
        return Err(MultiExpError::InvalidWindow);
    }

    let mut result = G::zero();

    for outer in 0..outer_count {
        let mut inner = 0usize;

        for i in 0..window {
            if exponent.bit(outer * window + i) {
                inner |= 1usize << i;
            }
        }

        if inner >= table[outer].len() {
            return Err(MultiExpError::InvalidWindow);
        }

        result = result.add(&table[outer][inner]);
    }

    Ok(result)
}

/// 批量固定窗口标量乘法。
pub fn batch_exp<G, S>(
    scalar_size: usize,
    window: usize,
    table: &WindowTable<G>,
    exponents: &[S],
) -> Result<Vec<G>, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    exponents
        .iter()
        .map(|e| windowed_exp(scalar_size, window, table, e))
        .collect()
}

/// 批量固定窗口标量乘法，并乘以系数。
pub fn batch_exp_with_coeff<G, S, F>(
    scalar_size: usize,
    window: usize,
    table: &WindowTable<G>,
    coeff: &S,
    exponents: &[S],
    multiply: F,
) -> Result<Vec<G>, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
    F: Fn(&S, &S) -> S,
{
    exponents
        .iter()
        .map(|e| {
            let scaled = multiply(coeff, e);
            windowed_exp(scalar_size, window, table, &scaled)
        })
        .collect()
}

/// 内积：sum(a[i] * b[i])。
pub fn inner_product<G, S>(a: &[G], b: &[S]) -> Result<G, MultiExpError>
where
    S: MsmScalar,
    G: MsmGroup<S>,
{
    multi_exp(a, b, MultiExpMethod::NaivePlain, 1)
}
