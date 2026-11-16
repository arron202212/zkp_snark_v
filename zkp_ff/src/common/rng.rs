use sha2::{Digest, Sha512};

/// 为 SHA512_rng 提供有限域所需的信息。
///
/// N 是模数的 limb 数，每个 limb 为 64 bit。
pub trait Sha512RngField<const N: usize>: Sized {
    /// 小端 limb 顺序的模数。
    const MODULUS: [u64; N];

    /// 从规范形式（小于模数）的 limb 构造有限域元素。
    fn from_canonical_limbs(limbs: [u64; N]) -> Self;
}

/// 比较两个小端 limb 表示的无符号整数。
///
/// 返回：
/// - Less：a < b
/// - Equal：a == b
/// - Greater：a > b
fn cmp_limbs<const N: usize>(
    a: &[u64; N],
    b: &[u64; N],
) -> std::cmp::Ordering {
    for i in (0..N).rev() {
        match a[i].cmp(&b[i]) {
            std::cmp::Ordering::Equal => continue,
            other => return other,
        }
    }

    std::cmp::Ordering::Equal
}

/// 计算模数的 bit 长度。
fn modulus_bit_length<const N: usize>(modulus: &[u64; N]) -> usize {
    for i in (0..N).rev() {
        if modulus[i] != 0 {
            return i * 64 + (64 - modulus[i].leading_zeros() as usize);
        }
    }

    0
}

/// 使用 SHA-512 生成确定性的有限域元素。
///
/// 相同的 idx 和相同的有限域参数会得到相同结果。
pub fn sha512_rng<F, const N: usize>(idx: u64) -> F
where
    F: Sha512RngField<N>,
{
    assert!(N > 0, "field must have at least one limb");
    assert!(N <= 8, "SHA-512 only provides 8 limbs of 64 bits");

    let modulus = F::MODULUS;
    let bit_len = modulus_bit_length(&modulus);

    assert!(bit_len > 0, "field modulus must be non-zero");
    assert!(bit_len <= 512, "field modulus exceeds SHA-512 output");

    let mut iter: u64 = 0;

    loop {
        // SHA512(idx || iter)
        //
        // 两个整数均使用 8 字节小端序。
        let mut hasher = Sha512::new();

        hasher.update(idx.to_le_bytes());
        hasher.update(iter.to_le_bytes());

        let digest = hasher.finalize();

        // 将摘要的前 N 个 64-bit 小端 limb 转换成整数。
        let mut limbs = [0u64; N];

        for i in 0..N {
            let offset = i * 8;

            limbs[i] = u64::from_le_bytes(
                digest[offset..offset + 8]
                    .try_into()
                    .expect("SHA-512 digest chunk must contain 8 bytes"),
            );
        }

        // 清除模数 bit 长度以上的所有位。
        //
        // 对应 C++：
        // while (FieldT::mod.test_bit(bitno) == false) {
        //     rval.data[part] &= ~(one << bit);
        //     bitno--;
        // }
        for bit in bit_len..(N * 64) {
            limbs[bit / 64] &= !(1u64 << (bit % 64));
        }

        // Rejection sampling：
        // 只有 candidate < modulus 时才接受。
        if cmp_limbs(&limbs, &modulus) == std::cmp::Ordering::Less {
            return F::from_canonical_limbs(limbs);
        }

        // C++ uint64_t 溢出时按无符号整数回绕。
        iter = iter.wrapping_add(1);
    }
}