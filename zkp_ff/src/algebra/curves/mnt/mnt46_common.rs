use std::sync::OnceLock;

use crate::algebra::field_utils::bigint::BigInt;

/// MNT4-6 曲线公共参数。
///
/// 模数位数。
pub const MNT46_A_BITCOUNT: usize = 298;
pub const MNT46_B_BITCOUNT: usize = 298;

/// 每个 limb 的位数。
pub const GMP_NUMB_BITS: usize = 64;

/// 模数所需 limb 数。
pub const MNT46_A_LIMBS: usize = (MNT46_A_BITCOUNT + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

pub const MNT46_B_LIMBS: usize = (MNT46_B_BITCOUNT + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;

/// MNT4-6 的 A 模数。
pub static MNT46_MODULUS_A: OnceLock<BigInt<MNT46_A_LIMBS>> = OnceLock::new();

/// MNT4-6 的 B 模数。
pub static MNT46_MODULUS_B: OnceLock<BigInt<MNT46_B_LIMBS>> = OnceLock::new();

/// 初始化 A 模数。
pub fn init_mnt46_modulus_a(modulus: BigInt<MNT46_A_LIMBS>) -> Result<(), BigInt<MNT46_A_LIMBS>> {
    MNT46_MODULUS_A.set(modulus)
}

/// 初始化 B 模数。
pub fn init_mnt46_modulus_b(modulus: BigInt<MNT46_B_LIMBS>) -> Result<(), BigInt<MNT46_B_LIMBS>> {
    MNT46_MODULUS_B.set(modulus)
}

/// 获取 A 模数。
pub fn mnt46_modulus_a() -> &'static BigInt<MNT46_A_LIMBS> {
    MNT46_MODULUS_A
        .get()
        .expect("MNT46 modulus A has not been initialized")
}

/// 获取 B 模数。
pub fn mnt46_modulus_b() -> &'static BigInt<MNT46_B_LIMBS> {
    MNT46_MODULUS_B
        .get()
        .expect("MNT46 modulus B has not been initialized")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mnt46_limb_count() {
        assert_eq!(MNT46_A_BITCOUNT, 298);
        assert_eq!(MNT46_B_BITCOUNT, 298);

        assert_eq!(MNT46_A_LIMBS, 5);
        assert_eq!(MNT46_B_LIMBS, 5);
    }
}
