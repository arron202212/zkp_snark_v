use crate::algebra::field_utils::bigint::BigInt;

/// Find the Windowed Non-Adjacent Form (wNAF) of a scalar.
///
/// The returned digits satisfy:
///
/// scalar = sum(digit[i] * 2^i)
///
/// For window size w:
/// - nonzero digits are odd;
/// - |digit| <= 2^w - 1;
/// - any two nonzero digits are separated by at least w zero digits.
///
/// Digits are stored in little-endian order.
pub fn find_wnaf<const N: usize>(
    window_size: usize,
    scalar: &BigInt<N>,
) -> Result<Vec<i64>, &'static str> {
    if window_size == 0 || window_size >= 63 {
        return Err("window_size must be in 1..63");
    }

    let length = BigInt::<N>::max_bits();

    let mut result = vec![0i64; length + 1];
    let mut c = *scalar;

    let modulus = 1u64 << (window_size + 1);
    let half = 1u64 << window_size;

    let mut j = 0usize;

    while !c.is_zero() {
        let mut u = 0i64;

        if (c.data[0] & 1) == 1 {
            let residue = c.data[0] & (modulus - 1);

            if residue > half {
                u = residue as i64 - modulus as i64;
            } else {
                u = residue as i64;
            }

            if u > 0 {
                sub_small(&mut c, u as u64)?;
            } else if u < 0 {
                add_small(&mut c, (-u) as u64)?;
            }
        }

        result[j] = u;
        j += 1;

        // c = floor(c / 2)
        shift_right_one(&mut c);
    }

    result.truncate(j);

    Ok(result)
}

/// Subtract a small integer from a fixed-width bigint.
fn sub_small<const N: usize>(value: &mut BigInt<N>, rhs: u64) -> Result<(), &'static str> {
    if N == 0 {
        return Err("cannot subtract from zero-width bigint");
    }

    let (first, mut borrow) = value.data[0].overflowing_sub(rhs);
    value.data[0] = first;

    for i in 1..N {
        if !borrow {
            break;
        }

        let (next, next_borrow) = value.data[i].overflowing_sub(1);
        value.data[i] = next;
        borrow = next_borrow;
    }

    if borrow {
        return Err("bigint subtraction underflow");
    }

    Ok(())
}

/// Add a small integer to a fixed-width bigint.
fn add_small<const N: usize>(value: &mut BigInt<N>, rhs: u64) -> Result<(), &'static str> {
    if N == 0 {
        return Err("cannot add to zero-width bigint");
    }

    let (first, mut carry) = value.data[0].overflowing_add(rhs);
    value.data[0] = first;

    for i in 1..N {
        if !carry {
            break;
        }

        let (next, next_carry) = value.data[i].overflowing_add(1);
        value.data[i] = next;
        carry = next_carry;
    }

    if carry {
        return Err("bigint addition overflow");
    }

    Ok(())
}

/// Shift a bigint right by one bit.
fn shift_right_one<const N: usize>(value: &mut BigInt<N>) {
    let mut carry = 0u64;

    for i in (0..N).rev() {
        let next_carry = value.data[i] & 1;

        value.data[i] = (value.data[i] >> 1) | (carry << 63);

        carry = next_carry;
    }
}

/// Group element interface required by wNAF scalar multiplication.
///
/// Implement this for G1, G2, or other additive groups.
pub trait WnafGroup: Clone {
    fn zero() -> Self;

    fn add(&self, rhs: &Self) -> Self;

    fn sub(&self, rhs: &Self) -> Self;

    fn double(&self) -> Self;
}

/// Fixed-window wNAF exponentiation.
///
/// This is equivalent to libff::fixed_window_wnaf_exp().
pub fn fixed_window_wnaf_exp<T, const N: usize>(
    window_size: usize,
    base: &T,
    scalar: &BigInt<N>,
) -> Result<T, &'static str>
where
    T: WnafGroup,
{
    if window_size == 0 || window_size >= usize::BITS as usize {
        return Err("invalid window size");
    }

    let naf = find_wnaf(window_size, scalar)?;

    let table_size = 1usize
        .checked_shl((window_size - 1) as u32)
        .ok_or("window size too large")?;

    // table[i] = (2*i + 1) * base
    let mut table = Vec::with_capacity(table_size);

    let mut tmp = base.clone();
    let dbl = base.double();

    for _ in 0..table_size {
        table.push(tmp.clone());
        tmp = tmp.add(&dbl);
    }

    let mut result = T::zero();
    let mut found_nonzero = false;

    for i in (0..naf.len()).rev() {
        if found_nonzero {
            result = result.double();
        }

        let digit = naf[i];

        if digit != 0 {
            found_nonzero = true;

            let index = (digit.unsigned_abs() / 2) as usize;

            let point = table
                .get(index)
                .ok_or("wNAF digit exceeds precomputed table")?;

            if digit > 0 {
                result = result.add(point);
            } else {
                result = result.sub(point);
            }
        }
    }

    Ok(result)
}

/// Select the wNAF window size using a group-specific lookup table.
///
/// The table must be sorted in ascending order, matching the
/// interpretation used by libff.
///
/// If no entry matches, use ordinary double-and-add.
pub fn opt_window_wnaf_exp<T, const N: usize>(
    base: &T,
    scalar: &BigInt<N>,
    scalar_bits: usize,
    window_table: &[usize],
) -> Result<T, &'static str>
where
    T: WnafGroup,
{
    let mut best = 0usize;

    for i in (0..window_table.len()).rev() {
        if scalar_bits >= window_table[i] {
            best = i + 1;
            break;
        }
    }

    if best > 0 {
        fixed_window_wnaf_exp(best, base, scalar)
    } else {
        Ok(double_and_add(base, scalar))
    }
}

/// Ordinary double-and-add scalar multiplication.
pub fn double_and_add<T, const N: usize>(base: &T, scalar: &BigInt<N>) -> T
where
    T: WnafGroup,
{
    let mut result = T::zero();

    for i in (0..scalar.num_bits()).rev() {
        result = result.double();

        if scalar.test_bit(i) {
            result = result.add(base);
        }
    }

    result
}
