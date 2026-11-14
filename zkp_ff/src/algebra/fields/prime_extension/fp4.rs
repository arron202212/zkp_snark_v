use super::super::prime_base::fp::Fp;
use super::fp2::{Fp2, Fp2Params};

/// Fp4 = Fp2[v] / (v² - u)
///
/// Fp2 = Fp[u] / (u² - non_residue)
///
/// 因此 v² = u，且 u² = non_residue。
#[derive(Clone, Debug)]
pub struct Fp4Params {
    pub fp2_params: Fp2Params,

    /// Fp2 中 u² 对应的底域非剩余元。
    pub non_residue: Fp,

    /// Fp4 的 Frobenius 系数。
    pub frobenius_coeffs_c1: [Fp2; 4],
}

#[derive(Clone, Debug)]
pub struct Fp4 {
    /// 常数项
    pub c0: Fp2,

    /// v 的系数
    pub c1: Fp2,

    pub params: Fp4Params,
}

impl Fp4Params {
    pub fn new(fp2_params: Fp2Params, non_residue: Fp, frobenius_coeffs_c1: [Fp2; 4]) -> Self {
        Self {
            fp2_params,
            non_residue,
            frobenius_coeffs_c1,
        }
    }
}

impl Fp4 {
    pub fn new(c0: Fp2, c1: Fp2, params: Fp4Params) -> Self {
        Self { c0, c1, params }
    }

    pub fn zero(params: Fp4Params) -> Self {
        let fp2_params = params.fp2_params.clone();

        Self {
            c0: Fp2::zero(fp2_params.clone()),
            c1: Fp2::zero(fp2_params),
            params,
        }
    }

    pub fn one(params: Fp4Params) -> Self {
        let fp2_params = params.fp2_params.clone();

        Self {
            c0: Fp2::one(fp2_params.clone()),
            c1: Fp2::zero(fp2_params),
            params,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    pub fn equals(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }

    /// 对应 libff:
    ///
    /// Fp4_model::mul_by_non_residue()
    ///
    /// u * (c0 + c1*u) = non_residue*c1 + c0*u
    ///
    /// 这里的 u 是 Fp2 的扩域生成元。
    pub fn mul_by_non_residue(&self, value: &Fp2) -> Fp2 {
        Fp2::new(
            value.c1.mul(&self.params.non_residue),
            value.c0.clone(),
            self.params.fp2_params.clone(),
        )
    }

    pub fn add(&self, other: &Self) -> Self {
        Self::new(
            self.c0.add(&other.c0),
            self.c1.add(&other.c1),
            self.params.clone(),
        )
    }

    pub fn sub(&self, other: &Self) -> Self {
        Self::new(
            self.c0.sub(&other.c0),
            self.c1.sub(&other.c1),
            self.params.clone(),
        )
    }

    pub fn negate(&self) -> Self {
        Self::new(self.c0.negate(), self.c1.negate(), self.params.clone())
    }

    /// Fp 标量乘法。
    pub fn mul_by_fp(&self, scalar: &Fp) -> Self {
        Self::new(
            self.c0.mul_by_fp(scalar),
            self.c1.mul_by_fp(scalar),
            self.params.clone(),
        )
    }

    /// Fp2 标量乘法。
    pub fn mul_by_fp2(&self, scalar: &Fp2) -> Self {
        Self::new(
            self.c0.mul(scalar),
            self.c1.mul(scalar),
            self.params.clone(),
        )
    }

    /// Karatsuba multiplication.
    ///
    /// (a + bv)(A + Bv)
    ///
    /// = (aA + u*bB)
    ///   + ((a+b)(A+B)-aA-bB)v
    ///
    /// where v² = u.
    pub fn mul(&self, other: &Self) -> Self {
        let a = &self.c0;
        let b = &self.c1;
        let a2 = &other.c0;
        let b2 = &other.c1;

        let aa = a.mul(a2);
        let bb = b.mul(b2);

        let beta_bb = self.mul_by_non_residue(&bb);

        let c0 = aa.add(&beta_bb);

        let c1 = a.add(b).mul(&a2.add(b2)).sub(&aa).sub(&bb);

        Self::new(c0, c1, self.params.clone())
    }

    /// 对应 libff 的 mul_by_023。
    ///
    /// 要求 other.c0.c1 == 0，即 other.c0 是 Fp 的嵌入元素。
    pub fn mul_by_023(&self, other: &Self) -> Self {
        assert!(
            other.c0.c1.is_zero(),
            "mul_by_023 requires other.c0.c1 == 0"
        );

        let a = &self.c0;
        let b = &self.c1;
        let a2 = &other.c0;
        let b2 = &other.c1;

        // a * A，其中 A 只有 Fp 常数项。
        let aa = Fp2::new(
            a.c0.mul(&a2.c0),
            a.c1.mul(&a2.c0),
            self.params.fp2_params.clone(),
        );

        let bb = b.mul(b2);
        let beta_bb = self.mul_by_non_residue(&bb);

        let c0 = aa.add(&beta_bb);

        let c1 = a.add(b).mul(&a2.add(b2)).sub(&aa).sub(&bb);

        Self::new(c0, c1, self.params.clone())
    }

    /// libff 的复数形式平方算法。
    pub fn square(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;

        let ab = a.mul(b);

        let nr_b = self.mul_by_non_residue(b);
        let nr_ab = self.mul_by_non_residue(&ab);

        let c0 = a.add(b).mul(&a.add(&nr_b)).sub(&ab).sub(&nr_ab);

        let c1 = ab.add(&ab);

        Self::new(c0, c1, self.params.clone())
    }

    pub fn squared(&self) -> Self {
        self.square()
    }

    /// Fp4 逆元。
    ///
    /// (a + bv)^-1
    /// = (a - bv) / (a² - u*b²)
    pub fn inverse(&self) -> Option<Self> {
        let b2 = self.c1.square();

        let denominator = self.c0.square().sub(&self.mul_by_non_residue(&b2));

        let denominator_inv = denominator.inverse()?;

        let c0 = self.c0.mul(&denominator_inv);
        let c1 = self.c1.mul(&denominator_inv).negate();

        Some(Self::new(c0, c1, self.params.clone()))
    }

    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(value) => {
                *self = value;
                true
            }
            None => false,
        }
    }

    pub fn unitary_inverse(&self) -> Self {
        Self::new(self.c0.clone(), self.c1.negate(), self.params.clone())
    }

    /// Frobenius map.
    ///
    /// c0 <- Frobenius(c0)
    /// c1 <- coeff[power % 4] * Frobenius(c1)
    pub fn frobenius_map(&self, power: usize) -> Self {
        let index = power % 4;

        let c0 = self.c0.frobenius_map(power);

        let c1 = self.params.frobenius_coeffs_c1[index].mul(&self.c1.frobenius_map(power));

        Self::new(c0, c1, self.params.clone())
    }

    /// 普通二进制幂运算。
    ///
    /// exponent 使用 little-endian u64 limbs。
    pub fn pow(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for &limb in exponent {
            let mut e = limb;

            for _ in 0..64 {
                if e & 1 == 1 {
                    result = result.mul(&base);
                }

                e >>= 1;
                base = base.square();
            }
        }

        result
    }

    /// Cyclotomic square.
    ///
    /// 仅适用于 cyclotomic subgroup 中的元素。
    pub fn cyclotomic_squared(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;

        let aa = b.square();
        let bb = b.add(a);
        let cc = bb.square().sub(&aa);

        let dd = self.mul_by_non_residue(&aa);
        let ee = cc.sub(&dd);
        let ff = dd.add(&dd).add(&Fp2::one(self.params.fp2_params.clone()));
        let gg = ee.sub(&Fp2::one(self.params.fp2_params.clone()));

        Self::new(ff, gg, self.params.clone())
    }

    /// Cyclotomic exponentiation.
    ///
    /// 这里使用普通二进制幂运算，保证语义正确。
    /// 若需要与 libff 一样使用 wNAF，可再替换为 wNAF 实现。
    pub fn cyclotomic_exp(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for &limb in exponent {
            let mut e = limb;

            for _ in 0..64 {
                if e & 1 == 1 {
                    result = result.mul(&base);
                }

                e >>= 1;
                base = base.cyclotomic_squared();
            }
        }

        result
    }

    /// 序列化为 16 个 u64。
    ///
    /// c0 || c1，每个 Fp2 为 8 个 u64。
    pub fn to_words(&self) -> Vec<u64> {
        let mut words = self.c0.to_words();
        words.extend(self.c1.to_words());
        words
    }

    /// 从 16 个 u64 反序列化。
    pub fn from_words(words: &[u64], params: Fp4Params) -> Option<Self> {
        if words.len() != 16 {
            return None;
        }

        let c0 = Fp2::from_words(&words[0..8], params.fp2_params.clone())?;

        let c1 = Fp2::from_words(&words[8..16], params.fp2_params.clone())?;

        Some(Self::new(c0, c1, params))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
    fn test_fp4_params() -> Fp4Params {
        let fp_params = FieldParams::new([7, 0, 0, 0]);

        // Fp2 = F7[u] / (u² - 3)
        // 3 在 F7 中不是平方，因此 Fp2 是扩域。
        let nr = Fp::from_u64(3, fp_params.clone());

        let fp2_params = Fp2Params::new(
            fp_params.clone(),
            nr.clone(),
            [
                Fp::one(fp_params.clone()),
                Fp::from_u64(6, fp_params.clone()),
            ],
        );

        let _zero = Fp2::zero(fp2_params.clone());
        let one = Fp2::one(fp2_params.clone());

        // Frobenius 系数需要根据具体扩域参数计算。
        // 此测试主要验证乘法、平方、求逆和序列化。
        Fp4Params::new(fp2_params, nr, [one.clone(), one.clone(), one.clone(), one])
    }

    #[test]
    fn test_fp4_zero_one() {
        let params = test_fp4_params();

        let zero = Fp4::zero(params.clone());
        let one = Fp4::one(params);

        assert!(zero.is_zero());
        assert!(one.is_one());
        assert!(!zero.is_one());
    }

    #[test]
    fn test_fp4_add_sub() {
        let params = test_fp4_params();

        let a = Fp4::one(params.clone());
        let b = a.add(&a);

        assert!(b.sub(&a).equals(&a));
    }

    #[test]
    fn test_fp4_mul_square() {
        let params = test_fp4_params();

        let a = Fp4::one(params.clone());

        assert!(a.mul(&a).equals(&a.square()));
    }

    #[test]
    fn test_fp4_inverse() {
        let params = test_fp4_params();

        let a = Fp4::one(params.clone());

        let inverse = a.inverse().unwrap();

        assert!(a.mul(&inverse).is_one());
    }

    #[test]
    fn test_fp4_serialization() {
        let params = test_fp4_params();

        let a = Fp4::one(params.clone());

        let words = a.to_words();

        assert_eq!(words.len(), 16);

        let b = Fp4::from_words(&words, params.clone()).unwrap();

        assert!(a.equals(&b));

        assert!(Fp4::from_words(&words[..8], params).is_none());
    }
}
