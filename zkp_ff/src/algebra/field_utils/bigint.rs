use core::cmp::Ordering;
use core::fmt;
use core::str::FromStr;

/// Fixed-width unsigned big integer.
///
/// Limbs use little-endian order:
/// data[0] is the least significant 64-bit limb.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BigInt<const N: usize> {
    pub data: [u64; N],
}

pub type BigInt256 = BigInt<4>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BigIntError {
    InvalidDecimal,
    Overflow,
    InvalidHex,
}

impl<const N: usize> Default for BigInt<N> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<const N: usize> BigInt<N> {
    /// Construct zero.
    pub const fn zero() -> Self {
        Self { data: [0; N] }
    }

    /// Construct one.
    pub fn one() -> Self {
        let mut result = Self::zero();

        if N > 0 {
            result.data[0] = 1;
        }

        result
    }

    /// Construct from u64.
    pub fn from_u64(value: u64) -> Self {
        let mut result = Self::zero();

        if N > 0 {
            result.data[0] = value;
        }

        result
    }

    /// Construct from little-endian limbs.
    pub fn from_limbs(data: [u64; N]) -> Self {
        Self { data }
    }

    /// Construct from little-endian words.
    ///
    /// Missing words are zero-filled.
    /// Extra words are rejected.
    pub fn from_words(words: &[u64]) -> Option<Self> {
        if words.len() > N {
            return None;
        }

        let mut result = Self::zero();

        result.data[..words.len()].copy_from_slice(words);

        Some(result)
    }

    /// Clear all limbs.
    pub fn clear(&mut self) {
        self.data.fill(0);
    }

    /// Check whether the integer is zero.
    pub fn is_zero(&self) -> bool {
        self.data.iter().all(|&x| x == 0)
    }

    /// Check whether the integer is even.
    pub fn is_even(&self) -> bool {
        N == 0 || (self.data[0] & 1) == 0
    }

    /// Return the maximum number of representable bits.
    pub const fn max_bits() -> usize {
        N * 64
    }

    /// Return the actual number of significant bits.
    ///
    /// Equivalent to libff::bigint::num_bits().
    pub fn num_bits(&self) -> usize {
        for i in (0..N).rev() {
            let x = self.data[i];

            if x != 0 {
                return i * 64 + (64 - x.leading_zeros() as usize);
            }
        }

        0
    }

    /// Return the lowest limb.
    ///
    /// Equivalent to libff::bigint::as_ulong().
    pub fn as_ulong(&self) -> u64 {
        if N == 0 { 0 } else { self.data[0] }
    }

    /// Test a bit.
    ///
    /// Bits outside the fixed-width representation are false.
    pub fn test_bit(&self, bitno: usize) -> bool {
        if bitno >= Self::max_bits() {
            return false;
        }

        let limb = bitno / 64;
        let bit = bitno % 64;

        ((self.data[limb] >> bit) & 1) != 0
    }

    /// Set a bit.
    pub fn set_bit(&mut self, bitno: usize) {
        if bitno >= Self::max_bits() {
            return;
        }

        let limb = bitno / 64;
        let bit = bitno % 64;

        self.data[limb] |= 1u64 << bit;
    }

    /// Clear a bit.
    pub fn clear_bit(&mut self, bitno: usize) {
        if bitno >= Self::max_bits() {
            return;
        }

        let limb = bitno / 64;
        let bit = bitno % 64;

        self.data[limb] &= !(1u64 << bit);
    }

    /// Return little-endian limbs.
    pub fn to_words(&self) -> Vec<u64> {
        self.data.to_vec()
    }

    /// Return the number of significant bytes.
    pub fn num_bytes(&self) -> usize {
        (self.num_bits() + 7) / 8
    }

    /// Compare two integers.
    ///
    /// Equivalent to GMP mpn_cmp().
    pub fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..N).rev() {
            match self.data[i].cmp(&other.data[i]) {
                Ordering::Equal => {}
                result => return result,
            }
        }

        Ordering::Equal
    }

    /// Checked addition.
    pub fn overflowing_add(&self, rhs: &Self) -> (Self, bool) {
        let mut result = Self::zero();
        let mut carry = false;

        for i in 0..N {
            let (a, c1) = self.data[i].overflowing_add(rhs.data[i]);
            let (b, c2) = a.overflowing_add(carry as u64);

            result.data[i] = b;
            carry = c1 || c2;
        }

        (result, carry)
    }

    /// Checked multiplication by a small integer.
    ///
    /// Returns true if the result overflowed.
    pub fn overflowing_mul_u64(&mut self, rhs: u64) -> bool {
        let mut carry = 0u128;

        for i in 0..N {
            let value = self.data[i] as u128 * rhs as u128 + carry;

            self.data[i] = value as u64;
            carry = value >> 64;
        }

        carry != 0
    }

    /// Add a small integer.
    ///
    /// Returns true if overflow occurred.
    pub fn overflowing_add_u64(&mut self, rhs: u64) -> bool {
        if N == 0 {
            return rhs != 0;
        }

        let (value, mut carry) = self.data[0].overflowing_add(rhs);
        self.data[0] = value;

        for i in 1..N {
            if !carry {
                break;
            }

            let (value, next_carry) = self.data[i].overflowing_add(1);
            self.data[i] = value;
            carry = next_carry;
        }

        carry
    }

    /// Divide by a small integer.
    ///
    /// Returns the remainder.
    fn div_rem_u64(&mut self, divisor: u64) -> u64 {
        assert!(divisor != 0);

        let mut remainder = 0u128;

        for i in (0..N).rev() {
            let value = (remainder << 64) | self.data[i] as u128;

            self.data[i] = (value / divisor as u128) as u64;
            remainder = value % divisor as u128;
        }

        remainder as u64
    }

    /// Parse an unsigned decimal integer.
    ///
    /// Equivalent to libff::bigint(const char*).
    pub fn from_decimal_str(s: &str) -> Result<Self, BigIntError> {
        if s.is_empty() {
            return Err(BigIntError::InvalidDecimal);
        }

        let mut result = Self::zero();

        for byte in s.bytes() {
            if !byte.is_ascii_digit() {
                return Err(BigIntError::InvalidDecimal);
            }

            let digit = (byte - b'0') as u64;

            let overflow_mul = result.overflowing_mul_u64(10);
            let overflow_add = result.overflowing_add_u64(digit);

            if overflow_mul || overflow_add {
                return Err(BigIntError::Overflow);
            }
        }

        Ok(result)
    }

    /// Convert to decimal string.
    ///
    /// Equivalent to libff::bigint::to_mpz() followed by output.
    pub fn to_decimal_string(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }

        let mut value = *self;
        let mut digits = Vec::new();

        while !value.is_zero() {
            let remainder = value.div_rem_u64(10);
            digits.push((b'0' + remainder as u8) as char);
        }

        digits.iter().rev().collect()
    }

    /// Convert to hexadecimal string.
    pub fn to_hex_string(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }

        let mut output = String::new();
        let mut started = false;

        for i in (0..N).rev() {
            if !started {
                if self.data[i] == 0 {
                    continue;
                }

                output.push_str(&format!("{:x}", self.data[i]));
                started = true;
            } else {
                output.push_str(&format!("{:016x}", self.data[i]));
            }
        }

        output
    }

    /// Serialize as fixed-width little-endian bytes.
    ///
    /// Equivalent to BINARY_OUTPUT mode in libff.
    pub fn to_le_bytes(&self) -> Vec<u8> {
        let mut result = Vec::with_capacity(N * 8);

        for limb in self.data {
            result.extend_from_slice(&limb.to_le_bytes());
        }

        result
    }

    /// Deserialize fixed-width little-endian bytes.
    pub fn from_le_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != N * 8 {
            return None;
        }

        let mut result = Self::zero();

        for i in 0..N {
            let start = i * 8;

            result.data[i] = u64::from_le_bytes(bytes[start..start + 8].try_into().ok()?);
        }

        Some(result)
    }

    /// Randomize using a caller-provided random-word generator.
    ///
    /// The caller is responsible for providing a suitable random source.
    ///
    /// This mirrors libff::bigint::randomize() without coupling the
    /// bigint implementation to a particular RNG crate.
    pub fn randomize_with<F>(&mut self, mut next_u64: F)
    where
        F: FnMut() -> u64,
    {
        for limb in self.data.iter_mut() {
            *limb = next_u64();
        }
    }

    /// Convert to a fixed-width hexadecimal representation.
    pub fn print_hex(&self) {
        println!("{}", self.to_hex_string());
    }

    /// Print decimal representation.
    pub fn print(&self) {
        println!("{}", self.to_decimal_string());
    }
}

impl<const N: usize> From<u64> for BigInt<N> {
    fn from(value: u64) -> Self {
        Self::from_u64(value)
    }
}

impl<const N: usize> FromStr for BigInt<N> {
    type Err = BigIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_decimal_str(s)
    }
}

impl<const N: usize> Ord for BigInt<N> {
    fn cmp(&self, other: &Self) -> Ordering {
        BigInt::cmp(self, other)
    }
}

impl<const N: usize> PartialOrd for BigInt<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> fmt::Display for BigInt<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_decimal_string())
    }
}

impl<const N: usize> fmt::LowerHex for BigInt<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_hex_string())
    }
}

impl<const N: usize> fmt::UpperHex for BigInt<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_hex_string().to_uppercase())
    }
}

/// Write a bigint in binary form.
pub fn write_binary<const N: usize>(
    value: &BigInt<N>,
    output: &mut impl std::io::Write,
) -> std::io::Result<()> {
    output.write_all(&value.to_le_bytes())
}

/// Read a bigint in binary form.
pub fn read_binary<const N: usize>(input: &mut impl std::io::Read) -> std::io::Result<BigInt<N>> {
    let mut bytes = vec![0u8; N * 8];

    input.read_exact(&mut bytes)?;

    BigInt::from_le_bytes(&bytes).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "invalid bigint byte length",
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    type U256 = BigInt<4>;

    #[test]
    fn test_zero_and_one() {
        assert!(U256::zero().is_zero());
        assert_eq!(U256::one().data[0], 1);
        assert!(!U256::one().is_zero());
    }

    #[test]
    fn test_from_u64() {
        let x = U256::from_u64(123456);

        assert_eq!(x.data[0], 123456);
        assert_eq!(x.data[1], 0);
        assert_eq!(x.num_bits(), 17);
    }

    #[test]
    fn test_limb_order() {
        let x = U256::from_limbs([
            0x0123456789abcdef,
            0xfedcba9876543210,
            0x1111222233334444,
            0xaaaabbbbccccdddd,
        ]);

        assert_eq!(x.data[0], 0x0123456789abcdef);
        assert_eq!(x.data[3], 0xaaaabbbbccccdddd);
    }

    #[test]
    fn test_num_bits() {
        assert_eq!(U256::zero().num_bits(), 0);
        assert_eq!(U256::from_u64(1).num_bits(), 1);
        assert_eq!(U256::from_u64(2).num_bits(), 2);
        assert_eq!(U256::from_u64(255).num_bits(), 8);
        assert_eq!(U256::from_u64(256).num_bits(), 9);

        let x = U256::from_limbs([0, 0, 0, 1 << 63]);

        assert_eq!(x.num_bits(), 256);
    }

    #[test]
    fn test_test_bit() {
        let x = U256::from_limbs([1, 1, 0, 1 << 63]);

        assert!(x.test_bit(0));
        assert!(x.test_bit(64));
        assert!(x.test_bit(255));

        assert!(!x.test_bit(1));
        assert!(!x.test_bit(128));
        assert!(!x.test_bit(256));
    }

    #[test]
    fn test_set_and_clear_bit() {
        let mut x = U256::zero();

        x.set_bit(100);
        assert!(x.test_bit(100));

        x.clear_bit(100);
        assert!(!x.test_bit(100));
    }

    #[test]
    fn test_even() {
        assert!(U256::from_u64(0).is_even());
        assert!(U256::from_u64(2).is_even());
        assert!(!U256::from_u64(1).is_even());
        assert!(!U256::from_u64(3).is_even());
    }

    #[test]
    fn test_comparison() {
        let a = U256::from_limbs([u64::MAX, 0, 0, 0]);
        let b = U256::from_limbs([0, 1, 0, 0]);

        assert!(a < b);
        assert!(b > a);
        assert_eq!(a, a);
    }

    #[test]
    fn test_decimal_roundtrip() {
        let values = [
            "0",
            "1",
            "9",
            "10",
            "123456789",
            "18446744073709551615",
            "340282366920938463463374607431768211455",
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        ];

        for value in values {
            let x = U256::from_decimal_str(value).unwrap();

            assert_eq!(x.to_decimal_string(), value);
        }
    }

    #[test]
    fn test_decimal_overflow() {
        let too_large =
            "115792089237316195423570985008687907853269984665640564039457584007913129639936";

        assert_eq!(
            U256::from_decimal_str(too_large),
            Err(BigIntError::Overflow)
        );
    }

    #[test]
    fn test_invalid_decimal() {
        assert_eq!(
            U256::from_decimal_str("12abc"),
            Err(BigIntError::InvalidDecimal)
        );

        assert_eq!(U256::from_decimal_str(""), Err(BigIntError::InvalidDecimal));
    }

    #[test]
    fn test_hex() {
        let x = U256::from_limbs([0x0123456789abcdef, 0, 0, 0]);

        assert_eq!(x.to_hex_string(), "123456789abcdef");
    }

    #[test]
    fn test_binary_roundtrip() {
        let x = U256::from_limbs([
            0x0123456789abcdef,
            0xfedcba9876543210,
            0x1111222233334444,
            0xaaaabbbbccccdddd,
        ]);

        let bytes = x.to_le_bytes();
        let restored = U256::from_le_bytes(&bytes).unwrap();

        assert_eq!(x, restored);
    }

    #[test]
    fn test_binary_stream() {
        let x = U256::from_limbs([1, 2, 3, 4]);

        let mut buffer = Vec::new();

        write_binary(&x, &mut buffer).unwrap();

        let restored: U256 = read_binary(&mut buffer.as_slice()).unwrap();

        assert_eq!(x, restored);
    }

    #[test]
    fn test_overflowing_add() {
        let a = U256::from_limbs([u64::MAX, 0, 0, 0]);
        let b = U256::from_u64(1);

        let (result, overflow) = a.overflowing_add(&b);

        assert!(!overflow);
        assert_eq!(result.data, [0, 1, 0, 0]);
    }

    #[test]
    fn test_overflowing_add_max() {
        let a = U256::from_limbs([u64::MAX; 4]);
        let b = U256::from_u64(1);

        let (result, overflow) = a.overflowing_add(&b);

        assert!(overflow);
        assert!(result.is_zero());
    }

    #[test]
    fn test_randomize() {
        let mut x = U256::zero();

        let mut i = 1u64;

        x.randomize_with(|| {
            let result = i;
            i += 1;
            result
        });

        assert_eq!(x.data, [1, 2, 3, 4]);
    }
}
