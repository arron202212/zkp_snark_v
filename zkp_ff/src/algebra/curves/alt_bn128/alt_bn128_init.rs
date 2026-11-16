use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};

use crate::algebra::fields::prime_extension::fp6_3over2::{Fp6, Fp6Params};

use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12Params;

pub fn bn254_fq_params() -> FieldParams<4> {
    FieldParams::new([
        0x3c20_8c16_d87c_fd47,
        0x9781_6a91_6871_ca8d,
        0xb850_45b6_8181_585d,
        0x3064_4e72_e131_a029,
    ])
}

pub fn bn254_fq2_params() -> Fp2Params<4> {
    let fp_params = bn254_fq_params();

    let one = Fp::<4>::one(fp_params.clone());
    let minus_one = one.negate();

    Fp2Params::new(fp_params, minus_one.clone(), [one, minus_one])
}
pub fn bn254_fp12_params() -> Fp12Params {
    let fp_params = bn254_fq_params();
    let fp2_params = bn254_fq2_params();

    // Fq2 non-residue: u^2 = -1
    // Fq6 non-residue: xi = 9 + u
    let xi = fq2_from_decimal("9", "1", &fp2_params);

    // -------------------------------------------------------------
    // Fq6 Frobenius coefficients
    // -------------------------------------------------------------

    let fq6_c1 = [
        ("1", "0"),
        (
            "21575463638280843010398324269430826099269044274347216827212613867836435027261",
            "10307601595873709700152284273816112264069230130616436755625194854815875713954",
        ),
        (
            "21888242871839275220042445260109153167277707414472061641714758635765020556616",
            "0",
        ),
        (
            "3772000881919853776433695186713858239009073593817195771773381919316419345261",
            "2236595495967245188281701248203181795121068902605861227855261137820944008926",
        ),
        (
            "2203960485148121921418603742825762020974279258880205651966",
            "0",
        ),
        (
            "18429021223477853657660792034369865839114504446431234726392080002137598044644",
            "9344045779998320333812420223237981029506012124075525679208581902008406485703",
        ),
    ];

    let fq6_c2 = [
        ("1", "0"),
        (
            "2581911344467009335267311115468803099551665605076196740867805258568234346338",
            "19937756971775647987995932169929341994314652964949448313374472400716661030",
        ),
        (
            "2203960485148121921418603742825762020974279258880205651966",
            "0",
        ),
        (
            "5324479202449903542726783395506214481928257762400643279780343368557297135718",
            "16208900380737693084919495127334387981393726419856888799917914180988844123039",
        ),
        (
            "21888242871839275220042445260109153167277707414472061641714758635765020556616",
            "0",
        ),
        (
            "13981852324922362344252311234282257507216387789820983642040889267519694726527",
            "7629828391165209371577384193250820201684255241773809077146787135900891633097",
        ),
    ];

    let fq6_params = Fp6Params {
        fp2_params: fp2_params.clone(),
        non_residue: xi,
        frobenius_coeffs_c1: fq6_c1.map(|(a, b)| fq2_from_decimal(a, b, &fp2_params)),
        frobenius_coeffs_c2: fq6_c2.map(|(a, b)| fq2_from_decimal(a, b, &fp2_params)),
    };

    // -------------------------------------------------------------
    // Fq12 Frobenius coefficients
    //
    // C++ stores each coefficient as Fq2.
    // Rust's Fp12Params stores each coefficient as Fp6.
    // Embed Fq2 into Fp6.c0.
    // -------------------------------------------------------------

    let fq12_coeffs = [
        ("1", "0"),
        (
            "8376118865763821496583973867626364092589906065868298776909617916018768340080",
            "16469823323077808223889137241176536799009286646108169935659301613961712198316",
        ),
        (
            "21888242871839275220042445260109153167277707414472061641714758635765020556617",
            "0",
        ),
        (
            "11697423496358154304825782922584725312912383441159505038794027105778954184319",
            "303847389135065887422783454877609941456349188919719272345083954437860409601",
        ),
        (
            "21888242871839275220042445260109153167277707414472061641714758635765020556616",
            "0",
        ),
        (
            "3321304630594332808241809054958361220322477375291206261884409189760185844239",
            "5722266937896532885780051958958348231143373700109372999374820235121374419868",
        ),
        (
            "21888242871839275222246405745257275088696311157297823662689037894645226208582",
            "0",
        ),
        (
            "13512124006075453725662431877630910996106405091429524885779419978626457868503",
            "5418419548761466998357268504080738289687024511189653727029736280683514010267",
        ),
        (
            "2203960485148121921418603742825762020974279258880205651966",
            "0",
        ),
        (
            "10190819375481120917420622822672549775783927716138318623895010788866272024264",
            "21584395482704209334823622290379665147239961968378104390343953940207365798982",
        ),
        (
            "2203960485148121921418603742825762020974279258880205651967",
            "0",
        ),
        (
            "18566938241244942414004596690298913868373833782006617400804628704885040364344",
            "16165975933942742336466353786298926857552937457188450663314217659523851788715",
        ),
    ];

    let zero2 = Fp2::zero(fp2_params.clone());

    let fq12_frobenius_coeffs = fq12_coeffs.map(|(a, b)| {
        let coeff = fq2_from_decimal(a, b, &fp2_params);

        Fp6::new(coeff, zero2.clone(), zero2.clone(), fq6_params.clone())
    });

    Fp12Params {
        fp6_params: fq6_params,
        frobenius_coeffs_c1: fq12_frobenius_coeffs,
    }
}
/// BN254 / alt_bn128 Ate loop count.
pub const ALT_BN128_ATE_LOOP_COUNT: [u64; 4] = [0x9d79_7039_be76_3ba8, 0x0000_0000_0000_0001, 0, 0];

/// C++: alt_bn128_ate_is_loop_count_neg = false
pub const ALT_BN128_ATE_IS_LOOP_COUNT_NEG: bool = false;

/// BN parameter z.
pub const ALT_BN128_FINAL_EXPONENT_Z: [u64; 4] = [0x44e9_92b4_4a69_09f1, 0, 0, 0];

/// C++: alt_bn128_final_exponent_is_z_neg = false
pub const ALT_BN128_FINAL_EXPONENT_IS_Z_NEG: bool = false;

/// G1 cofactor = 1.
pub const ALT_BN128_G1_COFACTOR: [u64; 1] = [1];

/// G2 cofactor.
pub const ALT_BN128_G2_COFACTOR: [u64; 4] = [
    0x345f_2299_c0f9_fa8d,
    0x06ce_ecda_572a_2489,
    0xb850_45b6_8181_585e,
    0x3064_4e72_e131_a029,
];

fn fq(value: [u64; 4], params: &FieldParams<4>) -> Fp<4> {
    Fp::new(value, params.clone())
}

fn fq2(c0: [u64; 4], c1: [u64; 4], params: &Fp2Params<4>) -> Fp2<4> {
    Fp2::new(
        fq(c0, &params.fp_params),
        fq(c1, &params.fp_params),
        params.clone(),
    )
}

/// 将十进制字符串转换为小端 4-limb 整数。
///
/// 输入必须是非负十进制整数，且小于 2^256。
fn decimal_to_limbs4(s: &str) -> [u64; 4] {
    let mut limbs = [0u64; 4];

    assert!(!s.is_empty(), "empty integer");

    for byte in s.bytes() {
        assert!(byte.is_ascii_digit(), "invalid decimal integer");

        let digit = (byte - b'0') as u128;
        let mut carry = digit;

        for limb in &mut limbs {
            let value = (*limb as u128) * 10 + carry;
            *limb = value as u64;
            carry = value >> 64;
        }

        assert_eq!(carry, 0, "integer exceeds 256 bits");
    }

    limbs
}

fn fq_from_decimal(s: &str, params: &FieldParams<4>) -> Fp<4> {
    Fp::new(decimal_to_limbs4(s), params.clone())
}

fn fq2_from_decimal(c0: &str, c1: &str, params: &Fp2Params<4>) -> Fp2<4> {
    Fp2::new(
        fq_from_decimal(c0, &params.fp_params),
        fq_from_decimal(c1, &params.fp_params),
        params.clone(),
    )
}

#[cfg(test)]
mod bn254_fp12_params_tests {
    use super::*;
    use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12;
    #[test]
    fn test_bn254_fp6_frobenius_period() {
        let params = bn254_fp12_params();
        let fp6_params = params.fp6_params.clone();

        let a = Fp6::new(
            Fp2::from_u64(7, fp6_params.fp2_params.clone()),
            Fp2::from_u64(11, fp6_params.fp2_params.clone()),
            Fp2::from_u64(13, fp6_params.fp2_params.clone()),
            fp6_params,
        );

        assert!(
            a.frobenius_map(6).equals(&a),
            "Fp6 Frobenius^6 must be identity"
        );
    }

    #[test]
    fn test_bn254_fp12_frobenius_period() {
        let params = bn254_fp12_params();

        let a = Fp12::new(
            Fp6::one(params.fp6_params.clone()),
            Fp6::new(
                Fp2::from_u64(7, params.fp6_params.fp2_params.clone()),
                Fp2::from_u64(11, params.fp6_params.fp2_params.clone()),
                Fp2::from_u64(13, params.fp6_params.fp2_params.clone()),
                params.fp6_params.clone(),
            ),
            params,
        );

        assert!(
            a.frobenius_map(12).equals(&a),
            "Fp12 Frobenius^12 must be identity"
        );
    }

    #[cfg(test)]
    mod bn254_fp12_params_tests {
        use super::*;

        #[test]
        fn test_bn254_fp12_params_init() {
            let params = bn254_fp12_params();

            assert_eq!(params.frobenius_coeffs_c1.len(), 12);

            assert_eq!(params.fp6_params.frobenius_coeffs_c1.len(), 6);

            assert_eq!(params.fp6_params.frobenius_coeffs_c2.len(), 6);

            assert_eq!(
                params.fp6_params.non_residue.c0,
                Fp::from_u64(9, params.fp6_params.fp2_params.fp_params.clone())
            );

            assert_eq!(
                params.fp6_params.non_residue.c1,
                Fp::one(params.fp6_params.fp2_params.fp_params.clone())
            );
        }
    }
}
