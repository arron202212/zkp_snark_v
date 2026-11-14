// src/curves/alt_bn128/params.rs
use crate::algebra::curves::alt_bn128::alt_bn128_g1::G1;
use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

pub const BN254_FQ_MODULUS: [u64; 4] = [
    0x3c20_8c16_d87c_fd47,
    0x9781_6a91_6871_ca8d,
    0xb850_45b6_8181_585d,
    0x3064_4e72_e131_a029,
];

pub const BN254_FR_MODULUS: [u64; 4] = [
    0x43e1_f593_f000_0001,
    0x2833_e848_79b9_7091,
    0xb850_45b6_8181_585d,
    0x3064_4e72_e131_a029,
];

pub const BN254_G1_COFACTOR: [u64; 1] = [1];

pub const BN254_G1_GENERATOR_X: [u64; 1] = [1];
pub const BN254_G1_GENERATOR_Y: [u64; 1] = [2];

// impl G1Params {
//     pub fn bn254() -> Self {
//         let fq = FieldParams::new(BN254_FQ_MODULUS);

//         let b = Fp::from_u64(3, fq);

//         Self {
//             fq,
//             b,
//             scalar_modulus: BN254_FR_MODULUS,
//         }
//     }

//     pub fn zero(&self) -> Fp {
//         Fp::zero(self.fq)
//     }

//     pub fn one(&self) -> Fp {
//         Fp::one(self.fq)
//     }

//     pub fn curve_b(&self) -> Fp {
//         self.b.clone()
//     }
// }

// src/curves/alt_bn128/g1.rs

// use crate::algebra::fields::prime_base::fp::{FieldParams, Fp};

#[derive(Clone, Debug)]
pub struct G1Params {
    pub fq: FieldParams,
    pub b: Fp,
    pub scalar_modulus: [u64; 4],
}

impl G1Params {
    pub fn bn254() -> Self {
        let fq = FieldParams::new([
            0x3c20_8c16_d87c_fd47,
            0x9781_6a91_6871_ca8d,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ]);

        let b = Fp::from_u64(3, fq.clone());

        Self {
            fq,
            b,
            scalar_modulus: [
                0x43e1_f593_f000_0001,
                0x2833_e848_79b9_7091,
                0xb850_45b6_8181_585d,
                0x3064_4e72_e131_a029,
            ],
        }
    }
}

// impl G1 {
//     /// libff: G1_zero = (0, 1, 0)
//     pub fn zero(params: G1Params) -> Self {
//         Self {
//             x: Fp::zero(params.fq.clone()),
//             y: Fp::one(params.fq.clone()),
//             z: Fp::zero(params.fq.clone()),
//             params,
//         }
//     }

//     /// libff: G1_one = (1, 2, 1)
//     ///
//     /// 注意：这里的 one 表示群生成元，不是群单位元。
//     pub fn generator(params: G1Params) -> Self {
//         Self {
//             x: Fp::from_u64(1, params.fq.clone()),
//             y: Fp::from_u64(2, params.fq.clone()),
//             z: Fp::one(params.fq.clone()),
//             params,
//         }
//     }

//     /// 对应 libff 的 G1_one()
//     pub fn one(params: G1Params) -> Self {
//         Self::generator(params)
//     }

//     pub fn mul_by_cofactor(&self) -> Self {
//         // BN254 G1 的 cofactor h = 1
//         self.clone()
//     }

//     pub fn field_char(&self) -> [u64; 4] {
//         self.params.fq.modulus
//     }

//     pub fn order(&self) -> [u64; 4] {
//         self.params.scalar_modulus
//     }

//     pub fn size_in_bits() -> usize {
//         255
//     }
// }

pub const G1_WNAF_WINDOW_TABLE: [usize; 4] = [11, 24, 60, 127];

pub const G1_FIXED_BASE_EXP_WINDOW_TABLE: [usize; 22] = [
    1, 5, 11, 32, 55, 162, 360, 815, 2373, 6978, 7122, 0, 57818, 0, 169679, 439759, 936073, 0,
    4666555, 7580404, 0, 34552892,
];

pub const BN254_ATE_LOOP_COUNT1: [u128; 1] = [29_793_968_203_157_093_288];

pub const BN254_FINAL_EXPONENT_Z1: [u128; 1] = [4_965_661_367_192_848_881];

pub const BN254_ATE_IS_LOOP_COUNT_NEG: bool = false;
pub const BN254_FINAL_EXPONENT_IS_Z_NEG: bool = false;

pub const BN254_ATE_LOOP_COUNT: [u64; 4] = [0x9f4f_1714_6b2a_4c48, 0x0000_0001, 0, 0];

pub const BN254_FINAL_EXPONENT_Z: [u64; 4] = [0x44e9_92b4_3a8c_3e11, 0, 0, 0];
