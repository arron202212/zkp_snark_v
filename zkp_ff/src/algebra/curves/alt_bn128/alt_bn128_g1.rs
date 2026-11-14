use crate::algebra::field_utils::bigint::BigInt;
use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
use crate::algebra::scalar_multiplication::multiexp::MsmGroup;

/// BN254 Fr 标量域模数，小端 u64 limbs。
pub const BN254_FR_MODULUS: [u64; 4] = [
    0x43e1_f593_f000_0001,
    0x2833_e848_79b9_7091,
    0xb850_45b6_8181_585d,
    0x3064_4e72_e131_a029,
];

/// BN254 Fq 基域模数，小端 u64 limbs。
pub const BN254_FQ_MODULUS: [u64; 4] = [
    0x3c20_8c16_d87c_fd47,
    0x9781_6a91_6871_ca8d,
    0xb850_45b6_8181_585d,
    0x3064_4e72_e131_a029,
];

/// BN254 G1 曲线参数。
///
/// 曲线方程：y² = x³ + b
#[derive(Clone, Debug)]
pub struct G1Params {
    pub fq: FieldParams,
    pub b: Fp,
}

impl G1Params {
    pub fn new(fq: FieldParams, b: Fp) -> Self {
        Self { fq, b }
    }

    pub fn zero(&self) -> Fp {
        Fp::zero(self.fq.clone())
    }

    pub fn one(&self) -> Fp {
        Fp::one(self.fq.clone())
    }

    /// BN254 G1 参数。
    ///
    /// Fq = BN254 base field
    /// E: y² = x³ + 3
    pub fn bn254() -> Self {
        let fq = FieldParams::new(BN254_FQ_MODULUS);
        let b = Fp::from_u64(3, fq.clone());

        Self { fq, b }
    }
}

/// 返回 BN254 G1 参数。
pub fn bn254_g1_params() -> G1Params {
    G1Params::bn254()
}

/// BN254 G1 点，使用 Jacobian 坐标。
///
/// affine 坐标：
///  ```text
///     x = X / Z²
///     y = Y / Z³
///
/// 无穷远点：
///
///     (0 : 1 : 0)
/// ```
#[derive(Clone, Debug)]
pub struct G1 {
    pub x: Fp,
    pub y: Fp,
    pub z: Fp,
    pub params: G1Params,
}

impl G1 {
    /// 使用 Jacobian 坐标构造点。
    pub fn new(x: Fp, y: Fp, z: Fp, params: G1Params) -> Self {
        Self { x, y, z, params }
    }

    /// 无穷远点。
    pub fn zero(params: G1Params) -> Self {
        Self {
            x: params.zero(),
            y: params.one(),
            z: params.zero(),
            params,
        }
    }

    /// BN254 G1 生成元。
    ///
    /// G = (1, 2)
    pub fn generator(params: G1Params) -> Self {
        let x = Fp::from_u64(1, params.fq.clone());
        let y = Fp::from_u64(2, params.fq.clone());

        Self::from_affine(x, y, params)
    }

    /// 根据 affine 坐标构造点。
    pub fn from_affine(x: Fp, y: Fp, params: G1Params) -> Self {
        let z = params.one();

        Self { x, y, z, params }
    }

    /// 根据 affine 坐标构造无穷远点或普通点。
    pub fn from_affine_option(point: Option<(Fp, Fp)>, params: G1Params) -> Self {
        match point {
            Some((x, y)) => Self::from_affine(x, y, params),
            None => Self::zero(params),
        }
    }

    /// 判断是否为无穷远点。
    pub fn is_zero(&self) -> bool {
        self.z.is_zero()
    }

    /// 判断是否为 special/affine 坐标。
    pub fn is_special(&self) -> bool {
        self.is_zero() || self.z == self.params.one()
    }

    /// Jacobian -> affine。
    ///
    /// x_affine = X / Z²
    /// y_affine = Y / Z³
    pub fn to_affine_coordinates(&mut self) {
        if self.is_zero() {
            self.x = self.params.zero();
            self.y = self.params.one();
            self.z = self.params.zero();
            return;
        }

        let z_inv = self.z.inverse().expect("non-zero Z must be invertible");

        let z2_inv = z_inv.square();
        let z3_inv = z2_inv.mul(&z_inv);

        self.x = self.x.mul(&z2_inv);
        self.y = self.y.mul(&z3_inv);
        self.z = self.params.one();
    }

    pub fn to_special(&mut self) {
        self.to_affine_coordinates();
    }

    pub fn into_affine(mut self) -> Self {
        self.to_affine_coordinates();
        self
    }

    /// 判断两个 Jacobian 点是否表示相同的 affine 点。
    pub fn equals(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero();
        }

        if other.is_zero() {
            return false;
        }

        let z1_2 = self.z.square();
        let z2_2 = other.z.square();

        let x1 = self.x.mul(&z2_2);
        let x2 = other.x.mul(&z1_2);

        if x1 != x2 {
            return false;
        }

        let z1_3 = self.z.mul(&z1_2);
        let z2_3 = other.z.mul(&z2_2);

        let y1 = self.y.mul(&z2_3);
        let y2 = other.y.mul(&z1_3);

        y1 == y2
    }

    /// 点取负：(X:Y:Z) -> (X:-Y:Z)。
    pub fn negate(&self) -> Self {
        Self {
            x: self.x.clone(),
            y: self.y.negate(),
            z: self.z.clone(),
            params: self.params.clone(),
        }
    }

    /// Jacobian 点加法。
    pub fn add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        let z1z1 = self.z.square();
        let z2z2 = other.z.square();

        let u1 = self.x.mul(&z2z2);
        let u2 = other.x.mul(&z1z1);

        let z1_3 = self.z.mul(&z1z1);
        let z2_3 = other.z.mul(&z2z2);

        let s1 = self.y.mul(&z2_3);
        let s2 = other.y.mul(&z1_3);

        // 相同 affine 横坐标。
        if u1 == u2 {
            // 同一个点：倍点。
            if s1 == s2 {
                return self.double();
            }

            // 互为相反点：无穷远点。
            return Self::zero(self.params.clone());
        }

        let h = u2.sub(&u1);
        let s2_minus_s1 = s2.sub(&s1);

        let two_h = h.add(&h);
        let i = two_h.square();

        let j = h.mul(&i);

        let r = s2_minus_s1.add(&s2_minus_s1);

        let v = u1.mul(&i);

        let x3 = r.square().sub(&j).sub(&v.add(&v));

        let s1j = s1.mul(&j);
        let y3 = r.mul(&v.sub(&x3)).sub(&s1j.add(&s1j));

        let z_sum = self.z.add(&other.z);

        let z3 = z_sum.square().sub(&z1z1).sub(&z2z2).mul(&h);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: self.params.clone(),
        }
    }

    /// Jacobian 倍点。
    ///
    /// 曲线 y² = x³ + b，其中 a = 0。
    pub fn double(&self) -> Self {
        if self.is_zero() || self.y.is_zero() {
            return Self::zero(self.params.clone());
        }

        let a = self.x.square();
        let b = self.y.square();
        let c = b.square();

        let x_plus_b = self.x.add(&b);

        let d0 = x_plus_b.square().sub(&a).sub(&c);

        let d = d0.add(&d0);

        let e = a.add(&a).add(&a);
        let f = e.square();

        let x3 = f.sub(&d.add(&d));

        let eight_c = c.add(&c).add(&c).add(&c).add(&c).add(&c).add(&c).add(&c);

        let y3 = e.mul(&d.sub(&x3)).sub(&eight_c);

        let yz = self.y.mul(&self.z);
        let z3 = yz.add(&yz);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: self.params.clone(),
        }
    }

    /// Jacobian + affine。
    ///
    /// 要求 other.Z == 1，或者 other 是无穷远点。
    pub fn mixed_add(&self, other: &Self) -> Self {
        debug_assert!(other.is_special());

        if self.is_zero() {
            return other.clone();
        }

        if other.is_zero() {
            return self.clone();
        }

        let z1z1 = self.z.square();

        let u1 = self.x.clone();
        let u2 = other.x.mul(&z1z1);

        let z1_3 = self.z.mul(&z1z1);

        let s1 = self.y.clone();
        let s2 = other.y.mul(&z1_3);

        if u1 == u2 {
            if s1 == s2 {
                return self.double();
            }

            return Self::zero(self.params.clone());
        }

        let h = u2.sub(&self.x);
        let hh = h.square();

        let mut i = hh.add(&hh);
        i = i.add(&i);

        let j = h.mul(&i);

        let mut r = s2.sub(&self.y);
        r = r.add(&r);

        let v = self.x.mul(&i);

        let x3 = r.square().sub(&j).sub(&v).sub(&v);

        let y1j = self.y.mul(&j);

        let y3 = r.mul(&v.sub(&x3)).sub(&y1j).sub(&y1j);

        let z3 = self.z.add(&h).square().sub(&z1z1).sub(&hh);

        Self {
            x: x3,
            y: y3,
            z: z3,
            params: self.params.clone(),
        }
    }

    /// 判断点是否满足曲线方程：
    ///
    /// Y² = X³ + b*Z⁶
    pub fn is_well_formed(&self) -> bool {
        if self.is_zero() {
            return true;
        }

        let x2 = self.x.square();
        let y2 = self.y.square();
        let z2 = self.z.square();

        let x3 = self.x.mul(&x2);

        let z3 = self.z.mul(&z2);
        let z6 = z3.square();

        let bz6 = self.params.b.mul(&z6);

        y2 == x3.add(&bz6)
    }

    /// BN254 G1 的 cofactor = 1。
    pub fn mul_by_cofactor(&self) -> Self {
        self.clone()
    }

    /// 检查曲线方程。
    ///
    /// 对 BN254 G1，cofactor = 1。
    pub fn validate(&self) -> bool {
        self.is_well_formed()
    }

    /// 输出 affine 坐标。
    pub fn affine_coordinates(&self) -> Option<(Fp, Fp)> {
        if self.is_zero() {
            return None;
        }

        let z_inv = self.z.inverse()?;
        let z2_inv = z_inv.square();
        let z3_inv = z2_inv.mul(&z_inv);

        Some((self.x.mul(&z2_inv), self.y.mul(&z3_inv)))
    }

    /// 使用 little-endian u64 limbs 进行标量乘法。
    ///
    /// scalar[0] 是最低有效 limb。
    pub fn mul_scalar(&self, scalar: &[u64]) -> Self {
        let mut result = Self::zero(self.params.clone());
        let mut base = self.clone();

        for &limb in scalar {
            let mut value = limb;

            for _ in 0..64 {
                if (value & 1) != 0 {
                    result = result.add(&base);
                }

                value >>= 1;
                base = base.double();
            }
        }

        result
    }

    /// 使用 BigInt<4> 进行标量乘法。
    pub fn scalar_mul(&self, scalar: &BigInt<4>) -> Self {
        if scalar.is_zero() || self.is_zero() {
            return Self::zero(self.params.clone());
        }

        let mut result = Self::zero(self.params.clone());

        for i in (0..scalar.num_bits()).rev() {
            result = result.double();

            if scalar.test_bit(i) {
                result = result.add(self);
            }
        }

        result
    }

    /// 将所有非零点转换为 affine/special 坐标。
    ///
    /// 无穷远点不参与批量求逆。
    pub fn batch_to_special_all_non_zeros(points: &mut [G1]) {
        let mut indices = Vec::new();
        let mut prefix_products = Vec::new();

        let mut acc: Option<Fp> = None;

        for (i, point) in points.iter().enumerate() {
            if point.is_zero() {
                continue;
            }

            let next = match acc.as_ref() {
                Some(product) => product.mul(&point.z),
                None => point.z.clone(),
            };

            prefix_products.push(next.clone());
            indices.push(i);
            acc = Some(next);
        }

        if indices.is_empty() {
            return;
        }

        let mut inverse = acc
            .expect("non-empty points")
            .inverse()
            .expect("product of non-zero field elements must be non-zero");

        for j in (0..indices.len()).rev() {
            let point_index = indices[j];

            let z_inverse = if j == 0 {
                inverse.clone()
            } else {
                prefix_products[j - 1].mul(&inverse)
            };

            inverse = inverse.mul(&points[point_index].z);

            let z2_inverse = z_inverse.square();
            let z3_inverse = z2_inverse.mul(&z_inverse);

            points[point_index].x = points[point_index].x.mul(&z2_inverse);

            points[point_index].y = points[point_index].y.mul(&z3_inverse);

            points[point_index].z = points[point_index].params.one();
        }
    }
}

impl PartialEq for G1 {
    fn eq(&self, other: &Self) -> bool {
        self.equals(other)
    }
}

impl Eq for G1 {}

/// 接入通用 MSM 框架。
impl MsmGroup<BigInt<4>> for G1 {
    fn zero() -> Self {
        G1::zero(bn254_g1_params())
    }

    fn add(&self, rhs: &Self) -> Self {
        G1::add(self, rhs)
    }

    fn double(&self) -> Self {
        G1::double(self)
    }

    fn mul_scalar(&self, scalar: &BigInt<4>) -> Self {
        self.scalar_mul(scalar)
    }

    fn is_zero(&self) -> bool {
        G1::is_zero(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_generator() -> G1 {
        G1::generator(bn254_g1_params())
    }

    #[test]
    fn test_g1_identity() {
        let params = bn254_g1_params();
        let zero = G1::zero(params);

        assert!(zero.is_zero());
        assert!(zero.is_well_formed());
        assert_eq!(zero, zero.clone());
    }

    #[test]
    fn test_g1_generator_is_well_formed() {
        let g = test_generator();

        assert!(!g.is_zero());
        assert!(g.is_well_formed());

        let (x, y) = g.affine_coordinates().unwrap();

        assert_eq!(x, Fp::from_u64(1, g.params.fq.clone()));
        assert_eq!(y, Fp::from_u64(2, g.params.fq.clone()));
    }

    #[test]
    fn test_g1_add_zero() {
        let g = test_generator();
        let zero = G1::zero(g.params.clone());

        assert_eq!(g.add(&zero), g);
        assert_eq!(zero.add(&g), g);
    }

    #[test]
    fn test_g1_double() {
        let g = test_generator();

        let doubled = g.add(&g);

        assert_eq!(doubled, g.double());
        assert!(doubled.is_well_formed());
    }

    #[test]
    fn test_g1_add_inverse() {
        let g = test_generator();

        let negative = g.negate();
        let result = g.add(&negative);

        assert!(result.is_zero());
        assert_eq!(result.x, result.params.zero());
        assert_eq!(result.y, result.params.one());
    }

    #[test]
    fn test_g1_affine_conversion() {
        let g = test_generator();

        let p = g.double();
        let mut affine = p.clone();

        affine.to_affine_coordinates();

        assert!(affine.is_special());
        assert_eq!(affine, p);
    }

    #[test]
    fn test_g1_mixed_add() {
        let g = test_generator();

        let mut affine = g.clone();
        affine.to_special();

        let mixed = g.mixed_add(&affine);
        let regular = g.add(&affine);

        assert_eq!(mixed, regular);
        assert!(mixed.is_well_formed());
    }

    #[test]
    fn test_g1_mixed_add_inverse() {
        let g = test_generator();

        let mut negative = g.negate();
        negative.to_special();

        let result = g.mixed_add(&negative);

        assert!(result.is_zero());
    }

    #[test]
    fn test_batch_to_special() {
        let g = test_generator();

        let p1 = g.double();
        let p2 = g.add(&p1);
        let zero = G1::zero(g.params.clone());

        let mut points = vec![p1.clone(), zero.clone(), p2.clone(), g.clone()];

        G1::batch_to_special_all_non_zeros(&mut points);

        assert!(points[0].is_special());
        assert!(points[1].is_zero());
        assert!(points[2].is_special());
        assert!(points[3].is_special());

        assert_eq!(points[0], p1);
        assert_eq!(points[1], zero);
        assert_eq!(points[2], p2);
        assert_eq!(points[3], g);
    }

    #[test]
    fn test_g1_scalar_mul_zero() {
        let g = test_generator();

        let result = g.mul_scalar(&[0, 0, 0, 0]);

        assert!(result.is_zero());
    }

    #[test]
    fn test_g1_scalar_mul_one() {
        let g = test_generator();

        let result = g.mul_scalar(&[1, 0, 0, 0]);

        assert_eq!(result, g);
    }

    #[test]
    fn test_g1_scalar_mul_two() {
        let g = test_generator();

        let result = g.mul_scalar(&[2, 0, 0, 0]);

        assert_eq!(result, g.double());
    }

    #[test]
    fn test_g1_scalar_mul_group_order() {
        let g = test_generator();

        let result = g.mul_scalar(&BN254_FR_MODULUS);

        assert!(result.is_zero());
    }

    #[test]
    fn test_g1_bigint_scalar_mul() {
        let g = test_generator();

        let zero = BigInt::<4>::zero();
        let one = BigInt::<4>::one();
        let two = BigInt::<4>::from_u64(2);

        assert!(g.scalar_mul(&zero).is_zero());
        assert_eq!(g.scalar_mul(&one), g);
        assert_eq!(g.scalar_mul(&two), g.double());
    }

    #[test]
    fn test_g1_bigint_scalar_mul_group_order() {
        let g = test_generator();
        let scalar = BigInt::<4>::from_limbs(BN254_FR_MODULUS);

        assert!(g.scalar_mul(&scalar).is_zero());
    }

    #[test]
    fn test_g1_msm_group_trait() {
        let g = test_generator();
        let scalar = BigInt::<4>::from_u64(2);

        let result = <G1 as MsmGroup<BigInt<4>>>::mul_scalar(&g, &scalar);

        assert_eq!(result, g.double());
    }
}
