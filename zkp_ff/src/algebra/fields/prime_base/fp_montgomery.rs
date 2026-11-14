/// Four-limb Montgomery arithmetic.
///
/// Limbs are stored in little-endian order:
/// value = limbs[0] + limbs[1] * 2^64 + ...
///
/// R = 2^256.

pub type Limbs = [u64; 4];

const LIMBS: usize = 4;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FieldParams {
    pub modulus: Limbs,
    /// -modulus[0]^{-1} mod 2^64
    pub n0_inv: u64,
    /// R mod p
    pub r: Limbs,
    /// R^2 mod p
    pub r2: Limbs,
}

impl FieldParams {
    pub fn new(modulus: Limbs) -> Self {
        assert!(modulus != [0; 4], "Modulus must be non-zero");
        assert!(modulus != [1, 0, 0, 0], "Modulus must be greater than one");
        assert!(modulus[0] & 1 == 1, "Modulus must be odd");

        // Four limbs provide 256-bit storage.
        // The actual modulus may be 254-bit, as in BN254.
        let n0_inv = inverse_mod_2_64(modulus[0]).wrapping_neg();

        // R = 2^256 mod p
        let mut r = [1, 0, 0, 0];

        for _ in 0..256 {
            r = add_mod(r, r, modulus);
        }

        // R^2 = 2^512 mod p
        let mut r2 = r;

        for _ in 0..256 {
            r2 = add_mod(r2, r2, modulus);
        }

        Self {
            modulus,
            n0_inv,
            r,
            r2,
        }
    }
}

/// Newton iteration for inverse modulo 2^64.
fn inverse_mod_2_64(x: u64) -> u64 {
    assert!(x & 1 == 1);

    let mut y = x;

    for _ in 0..6 {
        y = y.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(y)));
    }

    y
}

/// Compare two 256-bit integers.
pub fn cmp(a: &Limbs, b: &Limbs) -> core::cmp::Ordering {
    for i in (0..LIMBS).rev() {
        match a[i].cmp(&b[i]) {
            core::cmp::Ordering::Equal => {}
            other => return other,
        }
    }

    core::cmp::Ordering::Equal
}

/// Add two 256-bit integers.
/// Returns (result, carry).
pub fn add_raw(a: Limbs, b: Limbs) -> (Limbs, u64) {
    let mut out = [0u64; LIMBS];
    let mut carry = 0u128;

    for i in 0..LIMBS {
        let sum = a[i] as u128 + b[i] as u128 + carry;
        out[i] = sum as u64;
        carry = sum >> 64;
    }

    (out, carry as u64)
}

/// Subtract two 256-bit integers.
/// Returns (result, borrow).
pub fn sub_raw(a: Limbs, b: Limbs) -> (Limbs, u64) {
    let mut out = [0u64; LIMBS];
    let mut borrow = 0u128;

    for i in 0..LIMBS {
        let ai = a[i] as u128;
        let bi = b[i] as u128 + borrow;

        out[i] = ai.wrapping_sub(bi) as u64;
        borrow = (ai < bi) as u128;
    }

    (out, borrow as u64)
}

/// Modular addition.
/// Requires a < p and b < p.
pub fn add_mod(a: Limbs, b: Limbs, p: Limbs) -> Limbs {
    let (mut out, carry) = add_raw(a, b);

    if carry != 0 || cmp(&out, &p) != core::cmp::Ordering::Less {
        let (reduced, borrow) = sub_raw(out, p);

        // If carry = 0, subtraction must not borrow.
        // If carry = 1, subtraction must borrow.
        debug_assert_eq!(borrow, carry);

        out = reduced;
    }

    out
}

/// Modular subtraction.
/// Requires a < p and b < p.
pub fn sub_mod(a: Limbs, b: Limbs, p: Limbs) -> Limbs {
    let (mut out, borrow) = sub_raw(a, b);

    if borrow != 0 {
        let (reduced, carry) = add_raw(out, p);

        debug_assert_eq!(carry, 1);

        out = reduced;
    }

    out
}

/// CIOS Montgomery multiplication.
///
/// Computes:
///     a * b * R^-1 mod p
///
/// R = 2^256.
///
/// Requirements:
/// - p is odd
/// - a < p
/// - b < p
pub fn montgomery_mul(a: Limbs, b: Limbs, params: &FieldParams) -> Limbs {
    let p = &params.modulus;
    let n0_inv = params.n0_inv;

    assert!(cmp(&a, p) == core::cmp::Ordering::Less);
    assert!(cmp(&b, p) == core::cmp::Ordering::Less);

    // Extra limbs preserve intermediate carries.
    let mut t = [0u64; 6];

    for i in 0..LIMBS {
        // t += a * b[i]
        let mut carry = 0u128;

        for j in 0..LIMBS {
            let uv = (a[j] as u128) * (b[i] as u128) + (t[j] as u128) + carry;

            t[j] = uv as u64;
            carry = uv >> 64;
        }

        let (v, overflow) = t[4].overflowing_add(carry as u64);
        t[4] = v;

        if overflow {
            let (v, overflow2) = t[5].overflowing_add(1);
            t[5] = v;
            assert!(!overflow2, "CIOS carry overflow");
        }

        // m = t[0] * (-p[0]^-1) mod 2^64
        let m = t[0].wrapping_mul(n0_inv);

        // t += m * p
        carry = 0;

        for j in 0..LIMBS {
            let uv = (m as u128) * (p[j] as u128) + (t[j] as u128) + carry;

            t[j] = uv as u64;
            carry = uv >> 64;
        }

        let (v, overflow) = t[4].overflowing_add(carry as u64);
        t[4] = v;

        if overflow {
            let (v, overflow2) = t[5].overflowing_add(1);
            t[5] = v;
            assert!(!overflow2, "CIOS carry overflow");
        }

        debug_assert_eq!(t[0], 0);

        // Divide by 2^64.
        for j in 0..5 {
            t[j] = t[j + 1];
        }

        t[5] = 0;
    }

    let mut result = [t[0], t[1], t[2], t[3]];

    // Final result is less than 2p.
    if t[4] != 0 {
        debug_assert_eq!(t[4], 1);

        let (reduced, borrow) = sub_raw(result, *p);

        debug_assert_eq!(borrow, 1);

        result = reduced;
    } else if cmp(&result, p) != core::cmp::Ordering::Less {
        result = sub_raw(result, *p).0;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const FQ: Limbs = [
        0x3c208c16d87cfd47,
        0x97816a916871ca8d,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ];

    const FR: Limbs = [
        0x43e1f593f0000001,
        0x2833e84879b97091,
        0xb85045b68181585d,
        0x30644e72e131a029,
    ];

    #[test]
    fn test_bn254_modulus_width() {
        let fq = FieldParams::new(FQ);
        let fr = FieldParams::new(FR);

        assert_eq!(FQ[3] >> 63, 0);
        assert_eq!(FR[3] >> 63, 0);

        assert_eq!(fq.modulus, FQ);
        assert_eq!(fr.modulus, FR);
    }

    #[test]
    fn test_n0_inv() {
        let fq = FieldParams::new(FQ);
        let fr = FieldParams::new(FR);

        assert_eq!(FQ[0].wrapping_mul(fq.n0_inv), u64::MAX);

        assert_eq!(FR[0].wrapping_mul(fr.n0_inv), u64::MAX);
    }

    #[test]
    fn test_montgomery_parameters() {
        let fq = FieldParams::new(FQ);
        let fr = FieldParams::new(FR);

        assert_ne!(fq.r, [0; 4]);
        assert_ne!(fq.r2, [0; 4]);

        assert_ne!(fr.r, [0; 4]);
        assert_ne!(fr.r2, [0; 4]);
    }

    #[test]
    fn test_add_sub_mod() {
        let p = FieldParams::new(FQ).modulus;

        let a = [10, 0, 0, 0];
        let b = [20, 0, 0, 0];

        assert_eq!(add_mod(a, b, p), [30, 0, 0, 0]);
        assert_eq!(sub_mod(b, a, p), [10, 0, 0, 0]);
        assert_eq!(sub_mod(a, b, p), sub_raw(p, [10, 0, 0, 0]).0);
    }
}
