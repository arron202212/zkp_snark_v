use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
use crate::algebra::fields::prime_extension::fp2::Fp2Params;
use crate::algebra::fields::prime_extension::fp4::Fp4Params;

/// MNT4-298 scalar field.
pub type Mnt4Fr = Fp<5>;

/// MNT4-298 base field.
pub type Mnt4Fq = Fp<5>;

/// MNT4-298 quadratic extension.
pub type Mnt4Fq2 = crate::algebra::fields::prime_extension::fp2::Fp2<5>;

/// MNT4-298 quartic extension.
pub type Mnt4Fq4 = crate::algebra::fields::prime_extension::fp4::Fp4<5>;

/// MNT4 target group.
pub type Mnt4GT = Mnt4Fq4;

pub type Mnt4Fq4Params = Fp4Params<5>;

/// MNT4-298 scalar field modulus.
///
/// Decimal:
/// 475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137
pub const MNT4_MODULUS_R: [u64; 5] = [
    0xbb43_34a4_0000_0001,
    0xfb49_4c07_925d_6ad3,
    0xcaee_c963_5cf4_4194,
    0xa266_249d_a7b0_548e,
    0x0000_03bc_f7bc_d473,
];

/// MNT4-298 base field modulus.
///
/// Decimal:
/// 475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081
pub const MNT4_MODULUS_Q: [u64; 5] = [
    0xc90c_d65a_7166_0001,
    0x41a9_e35e_5120_0e12,
    0xcaee_c963_5d13_30ea,
    0xa266_249d_a7b0_548e,
    0x0000_03bc_f7bc_d473,
];

/// Initialize MNT4 scalar field Fr.
pub fn mnt4_fr_params() -> FieldParams<5> {
    FieldParams::new(MNT4_MODULUS_R)
}

/// Initialize MNT4 base field Fq.
pub fn mnt4_fq_params() -> FieldParams<5> {
    FieldParams::new(MNT4_MODULUS_Q)
}

/// Convert canonical Fq limbs into a Montgomery Fq element.
fn fq_from_words(words: [u64; 5], params: &FieldParams<5>) -> Mnt4Fq {
    Mnt4Fq::new(words, params.clone())
}

/// Fq2 = Fq[u] / (u² - 17).
///
/// Frobenius:
/// (a + bu)^q = a - bu
pub fn mnt4_fq2_params() -> Fp2Params<5> {
    let fq_params = mnt4_fq_params();

    let non_residue = Mnt4Fq::from_u64(17, fq_params.clone());

    let one = Mnt4Fq::one(fq_params.clone());
    let minus_one = one.negate();

    Fp2Params::new(fq_params, non_residue, [one, minus_one])
}

/// Fq4 Frobenius coefficients from libff mnt4_fields.cpp.
///
/// These are canonical integers, not Montgomery limbs.
const MNT4_FQ4_FROBENIUS_C1: [[u64; 5]; 4] = [
    [
        0x0000_0000_0000_0001,
        0x0000_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x0000_0000_0000_0000,
    ],
    [
        0x94dd_5d7d_ef69_80c4,
        0x8cd9_fae5_c1f7_bdcf,
        0x8d53_4beb_17da_f751,
        0x9916_dfdc_c2fd_1f96,
        0x0000_000f_7377_9fe0,
    ],
    [
        0xc90c_d65a_7166_0000,
        0x41a9_e35e_5120_0e12,
        0xcaee_c963_5d13_30ea,
        0xa266_249d_a7b0_548e,
        0x0000_03bc_f7bc_d473,
    ],
    [
        0x342f_78dc_81fc_7f3d,
        0xb4cf_e878_8f28_5043,
        0x3d9b_7d78_4538_3998,
        0x094f_44c0_e4b3_34f8,
        0x0000_03ad_8445_3493,
    ],
];

/// Build the four Fq4 Frobenius coefficients.
fn mnt4_fq4_frobenius_coeffs(fq_params: &FieldParams<5>) -> [Mnt4Fq; 4] {
    [
        fq_from_words(MNT4_FQ4_FROBENIUS_C1[0], fq_params),
        fq_from_words(MNT4_FQ4_FROBENIUS_C1[1], fq_params),
        fq_from_words(MNT4_FQ4_FROBENIUS_C1[2], fq_params),
        fq_from_words(MNT4_FQ4_FROBENIUS_C1[3], fq_params),
    ]
}

/// Fq4 = Fq2[v] / (v² - u).
///
/// Fq2 = Fq[u] / (u² - 17).
///
/// Therefore:
/// v² = u, u² = 17.
pub fn mnt4_fq4_params() -> Mnt4Fq4Params {
    let fq_params = mnt4_fq_params();
    let fq2_params = mnt4_fq2_params();

    // Corresponds to:
    // mnt4_Fq4::non_residue = mnt4_Fq("17");
    let non_residue = Mnt4Fq::from_u64(17, fq_params.clone());

    let frobenius_coeffs_c1 = mnt4_fq4_frobenius_coeffs(&fq_params);

    Mnt4Fq4Params::new(fq2_params, non_residue, frobenius_coeffs_c1)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_mnt4_field_moduli() {
        let fr = mnt4_fr_params();
        let fq = mnt4_fq_params();

        assert_eq!(fr.modulus[4], 0x0000_03bc_f7bc_d473);
        assert_eq!(fq.modulus[4], 0x0000_03bc_f7bc_d473);

        assert_eq!(fr.modulus, MNT4_MODULUS_R);

        assert_eq!(fq.modulus, MNT4_MODULUS_Q);
    }
}
