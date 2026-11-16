use std::sync::Arc;

use super::mnt4_fields::mnt4_fq_params;
use super::mnt4_g1::Mnt4G1Params;
use super::mnt4_g2::{Mnt4G2Params, mnt4_g2_params};

// ============================================================
// Pairing parameters
// ============================================================

/// MNT4 Ate loop count.
///
/// libff:
/// 689871209842287392837045615510547309923794944
///
/// Little-endian u64 limbs.
pub const MNT4_ATE_LOOP_COUNT: [u64; 5] = [
    0x0dc9_a1b6_7166_0000,
    0x4660_9756_bec2_a33f,
    0x0000_0000_001e_ef55,
    0,
    0,
];

/// libff:
/// mnt4_ate_is_loop_count_neg = false
pub const MNT4_ATE_IS_LOOP_COUNT_NEG: bool = false;

/// MNT4 final exponent.
///
/// C++ type:
/// bigint<4 * mnt4_q_limbs>
///
/// MNT4 q_limbs = 5, therefore 20 limbs.
///
/// The decimal value has 14 significant limbs.
/// The remaining 6 limbs are zero.
pub const MNT4_FINAL_EXPONENT: [u64; 20] = [
    0xe7e6_9541_c598_0000,
    0x065e_41b0_1242_2619,
    0x36c7_a580_1f55_ac8e,
    0xa570_04bb_9a1e_7a60,
    0xb351_384d_8485_e987,
    0x11ce_e431_be83_48de,
    0x4e61_f58b_68ea_590b,
    0x8507_4066_8186_8763,
    0xba73_85b8_d20f_ea4e,
    0x95e8_7dcb_6c97_234c,
    0x42e0_0292_64c0_324a,
    0x35ac_ca5a_07b2_394f,
    0xefe2_16b3_7afb_6d30,
    0x343c_7ac3_174c_87a1,
    0,
    0,
    0,
    0,
    0,
    0,
];

/// Absolute value of the last chunk w0.
///
/// libff:
/// 689871209842287392837045615510547309923794945
pub const MNT4_FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: [u64; 5] = [
    0x0dc9_a1b6_7166_0001,
    0x4660_9756_bec2_a33f,
    0x0000_0000_001e_ef55,
    0,
    0,
];

/// libff:
/// mnt4_final_exponent_last_chunk_is_w0_neg = false
pub const MNT4_FINAL_EXPONENT_LAST_CHUNK_IS_W0_NEG: bool = false;

/// Last chunk w1 = 1.
pub const MNT4_FINAL_EXPONENT_LAST_CHUNK_W1: [u64; 5] = [1, 0, 0, 0, 0];

// ============================================================
// Pairing parameter structure
// ============================================================

#[derive(Clone, Debug)]
pub struct Mnt4PairingParams {
    pub ate_loop_count: [u64; 5],

    pub ate_is_loop_count_neg: bool,

    pub final_exponent: [u64; 20],

    pub final_exponent_last_chunk_abs_of_w0: [u64; 5],

    pub final_exponent_last_chunk_is_w0_neg: bool,

    pub final_exponent_last_chunk_w1: [u64; 5],
}

impl Mnt4PairingParams {
    pub fn new() -> Self {
        Self {
            ate_loop_count: MNT4_ATE_LOOP_COUNT,

            ate_is_loop_count_neg: MNT4_ATE_IS_LOOP_COUNT_NEG,

            final_exponent: MNT4_FINAL_EXPONENT,

            final_exponent_last_chunk_abs_of_w0: MNT4_FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0,

            final_exponent_last_chunk_is_w0_neg: MNT4_FINAL_EXPONENT_LAST_CHUNK_IS_W0_NEG,

            final_exponent_last_chunk_w1: MNT4_FINAL_EXPONENT_LAST_CHUNK_W1,
        }
    }
}

impl Default for Mnt4PairingParams {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================
// Unified MNT4 parameters
// ============================================================

/// Initialized MNT4 curve and pairing parameters.
///
/// This corresponds to the parameters initialized by:
///
/// libff::init_mnt4_params()
pub struct Mnt4InitializedParams {
    /// G1 curve parameters.
    pub g1: Arc<Mnt4G1Params>,

    /// G2 curve and twist parameters.
    pub g2: Mnt4G2Params,

    /// Pairing parameters.
    pub pairing: Mnt4PairingParams,
}

// ============================================================
// Initialization
// ============================================================

/// Initialize all MNT4 parameters.
///
/// Equivalent to libff:
///
/// void init_mnt4_params()
///
/// The field parameters are initialized through mnt4_fields,
/// while G1 and G2 parameters are delegated to their own modules.
pub fn init_mnt4_params() -> Mnt4InitializedParams {
    // --------------------------------------------------------
    // G1
    // --------------------------------------------------------

    let fq_params = mnt4_fq_params();

    let g1 = Arc::new(Mnt4G1Params::new(fq_params));

    // --------------------------------------------------------
    // G2
    // --------------------------------------------------------

    // mnt4_g2_params() initializes:
    //
    // - twist
    // - twist_coeff_a
    // - twist_coeff_b
    // - twist_mul_by_a_c0 / c1
    // - twist_mul_by_b_c0 / c1
    // - twist_mul_by_q_x / y
    // - G2 cofactor
    //
    // These values should be defined in mnt4_g2.rs.
    let g2 = mnt4_g2_params();

    // --------------------------------------------------------
    // Pairing
    // --------------------------------------------------------

    let pairing = Mnt4PairingParams::new();

    Mnt4InitializedParams { g1, g2, pairing }
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mnt4_ate_loop_count() {
        let params = Mnt4PairingParams::new();

        assert_eq!(
            params.ate_loop_count,
            [
                0x0dc9_a1b6_7166_0000,
                0x4660_9756_bec2_a33f,
                0x0000_0000_001e_ef55,
                0,
                0,
            ]
        );

        assert!(!params.ate_is_loop_count_neg);
    }

    #[test]
    fn test_mnt4_final_exponent() {
        let params = Mnt4PairingParams::new();

        assert_eq!(params.final_exponent.len(), 20);

        assert_eq!(params.final_exponent, MNT4_FINAL_EXPONENT);

        assert_eq!(params.final_exponent[14..], [0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_mnt4_final_exponent_last_chunk() {
        let params = Mnt4PairingParams::new();

        assert_eq!(
            params.final_exponent_last_chunk_abs_of_w0,
            [
                0x0dc9_a1b6_7166_0001,
                0x4660_9756_bec2_a33f,
                0x0000_0000_001e_ef55,
                0,
                0,
            ]
        );

        assert!(!params.final_exponent_last_chunk_is_w0_neg);

        assert_eq!(params.final_exponent_last_chunk_w1, [1, 0, 0, 0, 0]);
    }

    #[test]
    fn test_mnt4_init_params() {
        let params = init_mnt4_params();

        assert_eq!(params.pairing.ate_loop_count, MNT4_ATE_LOOP_COUNT);

        assert_eq!(params.pairing.final_exponent, MNT4_FINAL_EXPONENT);

        assert_eq!(params.pairing.final_exponent_last_chunk_w1, [1, 0, 0, 0, 0]);
    }
}
