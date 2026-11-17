use crate::algebra::fields::prime_base::fp::Fp;
use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};
use crate::algebra::fields::prime_extension::fp6_3over2::{Fp6, Fp6Params};

/// Fp12 = Fp6[v] / (v^2 - w)
///
/// Fp6 = Fp2[w] / (w^3 - xi)
///
/// Element:
///     c0 + c1 * v
///
/// Fp6 coefficient layout:
///     c0 = a0 + a1*w + a2*w^2
///     c1 = b0 + b1*w + b2*w^2
#[derive(Clone)]
pub struct Fp12Params {
    pub fp6_params: Fp6Params,

    /// Frobenius coefficients for the c1 component.
    ///
    /// coeff[k] * Frobenius(c1, k)
    pub frobenius_coeffs_c1: [Fp6; 12],
}

#[derive(Clone)]
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
    pub params: Fp12Params,
}

impl Fp12 {
    pub fn new(c0: Fp6, c1: Fp6, params: Fp12Params) -> Self {
        Self { c0, c1, params }
    }

    pub fn zero(params: Fp12Params) -> Self {
        Self::new(
            Fp6::zero(params.fp6_params.clone()),
            Fp6::zero(params.fp6_params.clone()),
            params,
        )
    }

    pub fn one(params: Fp12Params) -> Self {
        Self::new(
            Fp6::one(params.fp6_params.clone()),
            Fp6::zero(params.fp6_params.clone()),
            params,
        )
    }

    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    pub fn equals(&self, other: &Self) -> bool {
        self.c0.equals(&other.c0) && self.c1.equals(&other.c1)
    }

    // -------------------------------------------------------------------------
    // Multiplication by the quadratic non-residue.
    //
    // In Fp6:
    //     w^3 = xi
    //
    // Therefore:
    //     w * (x0 + x1*w + x2*w^2)
    //       = xi*x2 + x0*w + x1*w^2
    //
    // This exactly matches libff's:
    //     Fp6(non_residue * elt.c2, elt.c0, elt.c1)
    // -------------------------------------------------------------------------

    pub fn mul_by_non_residue(&self, elt: &Fp6) -> Fp6 {
        let xi = &self.params.fp6_params.non_residue;

        Fp6::new(
            xi.mul(&elt.c2),
            elt.c0.clone(),
            elt.c1.clone(),
            self.params.fp6_params.clone(),
        )
    }

    // -------------------------------------------------------------------------
    // Addition / subtraction / negation
    // -------------------------------------------------------------------------

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

    // -------------------------------------------------------------------------
    // Scalar multiplication
    // -------------------------------------------------------------------------

    /// Fp * Fp12
    pub fn mul_by_fp(&self, scalar: &Fp) -> Self {
        let fp6 = |x: &Fp6| {
            Fp6::new(
                x.c0.mul_by_fp(scalar),
                x.c1.mul_by_fp(scalar),
                x.c2.mul_by_fp(scalar),
                self.params.fp6_params.clone(),
            )
        };

        Self::new(fp6(&self.c0), fp6(&self.c1), self.params.clone())
    }

    /// Fp2 * Fp12
    ///
    /// This corresponds to libff's operator* overload.
    pub fn mul_by_fp2(&self, scalar: &Fp2) -> Self {
        let fp6 = |x: &Fp6| {
            Fp6::new(
                scalar.mul(&x.c0),
                scalar.mul(&x.c1),
                scalar.mul(&x.c2),
                self.params.fp6_params.clone(),
            )
        };

        Self::new(fp6(&self.c0), fp6(&self.c1), self.params.clone())
    }

    /// Fp6 * Fp12
    pub fn mul_by_fp6(&self, scalar: &Fp6) -> Self {
        Self::new(
            scalar.mul(&self.c0),
            scalar.mul(&self.c1),
            self.params.clone(),
        )
    }

    // -------------------------------------------------------------------------
    // Karatsuba multiplication
    //
    // (a + bv)(A + Bv)
    //   = aA + bB*w + ((a+b)(A+B)-aA-bB)v
    //
    // v^2 = w
    // -------------------------------------------------------------------------

    pub fn mul(&self, other: &Self) -> Self {
        let a = &self.c0;
        let b = &self.c1;

        let aa = &other.c0;
        let bb = &other.c1;

        let a_a = a.mul(aa);
        let b_b = b.mul(bb);

        let r0 = a_a.add(&self.mul_by_non_residue(&b_b));

        let r1 = a.add(b).mul(&aa.add(bb)).sub(&a_a).sub(&b_b);

        Self::new(r0, r1, self.params.clone())
    }

    // -------------------------------------------------------------------------
    // Karatsuba squaring
    //
    // This is the direct translation of squared_karatsuba().
    // -------------------------------------------------------------------------

    pub fn squared_karatsuba(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;

        let a_sq = a.squared();
        let b_sq = b.squared();

        let r0 = a_sq.add(&self.mul_by_non_residue(&b_sq));

        let r1 = a.add(b).squared().sub(&a_sq).sub(&b_sq);

        Self::new(r0, r1, self.params.clone())
    }

    // -------------------------------------------------------------------------
    // Complex squaring
    //
    // libff's default squared() calls squared_complex().
    //
    // (a + bv)^2:
    //   c0 = (a+b)(a+w*b) - ab - w*ab
    //   c1 = 2ab
    // -------------------------------------------------------------------------

    pub fn squared_complex(&self) -> Self {
        let a = &self.c0;
        let b = &self.c1;

        let ab = a.mul(b);

        let a_plus_b = a.add(b);
        let a_plus_wb = a.add(&self.mul_by_non_residue(b));

        let r0 = a_plus_b
            .mul(&a_plus_wb)
            .sub(&ab)
            .sub(&self.mul_by_non_residue(&ab));

        let r1 = ab.add(&ab);

        Self::new(r0, r1, self.params.clone())
    }

    /// Same behavior as libff's squared().
    pub fn squared(&self) -> Self {
        self.squared_complex()
    }

    pub fn square(&mut self) {
        *self = self.squared();
    }

    // -------------------------------------------------------------------------
    // Inverse
    //
    // Algorithm 8:
    //
    // (a + bv)^-1 = (a - bv) / (a^2 - w*b^2)
    // -------------------------------------------------------------------------

    pub fn inverse(&self) -> Option<Self> {
        let a = &self.c0;
        let b = &self.c1;

        let t0 = a.squared();
        let t1 = b.squared();

        let t2 = t0.sub(&self.mul_by_non_residue(&t1));

        let t3 = t2.inverse()?;

        let c0 = a.mul(&t3);
        let c1 = b.mul(&t3).negate();

        Some(Self::new(c0, c1, self.params.clone()))
    }

    pub fn invert(&mut self) -> bool {
        match self.inverse() {
            Some(inv) => {
                *self = inv;
                true
            }
            None => false,
        }
    }

    // -------------------------------------------------------------------------
    // Frobenius map
    //
    // c0' = Frobenius(c0, power)
    // c1' = coeff[power % 12] * Frobenius(c1, power)
    // -------------------------------------------------------------------------

    pub fn frobenius_map(&self, power: usize) -> Self {
        let k = power % 12;

        Self::new(
            self.c0.frobenius_map(power),
            self.params.frobenius_coeffs_c1[k].mul(&self.c1.frobenius_map(power)),
            self.params.clone(),
        )
    }

    // -------------------------------------------------------------------------
    // Unitary inverse
    //
    // For elements in the unitary subgroup:
    //     (a + bv)^-1 = a - bv
    // -------------------------------------------------------------------------

    pub fn unitary_inverse(&self) -> Self {
        Self::new(self.c0.clone(), self.c1.negate(), self.params.clone())
    }

    // -------------------------------------------------------------------------
    // Cyclotomic squaring
    //
    // This is the optimized six-Fp2 implementation in libff.
    //
    // Coefficient mapping:
    //
    // z0 = c0.c0
    // z4 = c0.c1
    // z3 = c0.c2
    // z2 = c1.c0
    // z1 = c1.c1
    // z5 = c1.c2
    // -------------------------------------------------------------------------

    pub fn cyclotomic_squared(&self) -> Self {
        let xi = &self.params.fp6_params.non_residue;

        let z0 = self.c0.c0.clone();
        let z4 = self.c0.c1.clone();
        let z3 = self.c0.c2.clone();

        let z2 = self.c1.c0.clone();
        let z1 = self.c1.c1.clone();
        let z5 = self.c1.c2.clone();

        // a^2 = t0 + t1*y
        let tmp0 = z0.mul(&z1);

        let t0 = z0
            .add(&z1)
            .mul(&z0.add(&xi.mul(&z1)))
            .sub(&tmp0)
            .sub(&xi.mul(&tmp0));

        let t1 = tmp0.add(&tmp0);

        // b^2 = t2 + t3*y
        let tmp1 = z2.mul(&z3);

        let t2 = z2
            .add(&z3)
            .mul(&z2.add(&xi.mul(&z3)))
            .sub(&tmp1)
            .sub(&xi.mul(&tmp1));

        let t3 = tmp1.add(&tmp1);

        // c^2 = t4 + t5*y
        let tmp2 = z4.mul(&z5);

        let t4 = z4
            .add(&z5)
            .mul(&z4.add(&xi.mul(&z5)))
            .sub(&tmp2)
            .sub(&xi.mul(&tmp2));

        let t5 = tmp2.add(&tmp2);

        // A
        let new_z0 = three_mul_sub_two(&t0, &z0);
        let new_z1 = three_mul_add_two(&t1, &z1);

        // B
        let xi_t5 = xi.mul(&t5);
        let new_z2 = three_mul_add_two(&xi_t5, &z2);
        let new_z3 = three_mul_sub_two(&t4, &z3);

        // C
        let new_z4 = three_mul_sub_two(&t2, &z4);
        let new_z5 = three_mul_add_two(&t3, &z5);

        let c0 = Fp6::new(new_z0, new_z4, new_z3, self.params.fp6_params.clone());

        let c1 = Fp6::new(new_z2, new_z1, new_z5, self.params.fp6_params.clone());

        Self::new(c0, c1, self.params.clone())
    }

    // -------------------------------------------------------------------------
    // Sparse multiplication by 045.
    //
    // Input:
    //   ell_0
    //   ell_VW
    //   ell_VV
    //
    // This follows libff's optimized implementation.
    // -------------------------------------------------------------------------

    pub fn mul_by_045(&self, ell_0: &Fp2, ell_vw: &Fp2, ell_vv: &Fp2) -> Self {
        let z0 = self.c0.c0.clone();
        let z1 = self.c0.c1.clone();
        let z2 = self.c0.c2.clone();

        let z3 = self.c1.c0.clone();
        let z4 = self.c1.c1.clone();
        let z5 = self.c1.c2.clone();

        let x0 = ell_vw;
        let x4 = ell_0;
        let x5 = ell_vv;

        let xi = &self.params.fp6_params.non_residue;

        let tmp1 = xi.mul(x4);
        let tmp2 = xi.mul(x5);

        let t0 = x0.mul(&z0).add(&tmp1.mul(&z4)).add(&tmp2.mul(&z3));

        let t1 = x0.mul(&z1).add(&tmp1.mul(&z5)).add(&tmp2.mul(&z4));

        let t2 = x0.mul(&z2).add(&x4.mul(&z3)).add(&tmp2.mul(&z5));

        let t3 = x0.mul(&z3).add(&tmp1.mul(&z2)).add(&tmp2.mul(&z1));

        let t4 = x0.mul(&z4).add(&x4.mul(&z0)).add(&tmp2.mul(&z2));

        let t5 = x0.mul(&z5).add(&x4.mul(&z1)).add(&x5.mul(&z0));

        let c0 = Fp6::new(t0, t1, t2, self.params.fp6_params.clone());

        let c1 = Fp6::new(t3, t4, t5, self.params.fp6_params.clone());

        Self::new(c0, c1, self.params.clone())
    }

    // -------------------------------------------------------------------------
    // Sparse multiplication by 024.
    //
    // Input:
    //   ell_0
    //   ell_VW
    //   ell_VV
    //
    // This is the optimized sparse multiplication from libff.
    // -------------------------------------------------------------------------

    pub fn mul_by_024(&self, ell_0: &Fp2, ell_vw: &Fp2, ell_vv: &Fp2) -> Self {
        let mut z0 = self.c0.c0.clone();
        let mut z1 = self.c0.c1.clone();
        let mut z2 = self.c0.c2.clone();

        let mut z3 = self.c1.c0.clone();
        let mut z4 = self.c1.c1.clone();
        let mut z5 = self.c1.c2.clone();

        let x0 = ell_0;
        let x2 = ell_vv;
        let x4 = ell_vw;

        let xi = &self.params.fp6_params.non_residue;

        let d0 = z0.mul(x0);
        let d2 = z2.mul(x2);
        let d4 = z4.mul(x4);

        let t2 = z0.add(&z4);
        let mut t1 = z0.add(&z2);

        let s0 = z1.add(&z3).add(&z5);

        // For z0
        let mut s1 = z1.mul(x2);
        let t3 = s1.add(&d4);
        let t4 = xi.mul(&t3).add(&d0);

        z0 = t4;

        // For z1
        let t3 = z5.mul(x4);
        s1 = s1.add(&t3);

        let t3 = t3.add(&d2);
        let mut t4 = xi.mul(&t3);

        let t3 = z1.mul(x0);
        s1 = s1.add(&t3);
        t4 = t4.add(&t3);

        z1 = t4;

        // For z2
        let t0 = x0.add(x2);
        let t3 = t1.mul(&t0).sub(&d0).sub(&d2);

        let t4 = z3.mul(x4);
        s1 = s1.add(&t4);

        let t3 = t3.add(&t4);

        // For z3
        let t0 = z2.add(&z4);
        z2 = t3;

        t1 = x2.add(x4);

        let t3 = t0.mul(&t1).sub(&d2).sub(&d4);
        let mut t4 = xi.mul(&t3);

        let t3 = z3.mul(x0);
        s1 = s1.add(&t3);
        t4 = t4.add(&t3);

        z3 = t4;

        // For z4
        let t3 = z5.mul(x2);
        s1 = s1.add(&t3);

        let mut t4 = xi.mul(&t3);

        let t0 = x0.add(x4);
        let t3 = t2.mul(&t0).sub(&d0).sub(&d4);

        t4 = t4.add(&t3);
        z4 = t4;

        // For z5
        let t0 = x0.add(x2).add(x4);
        let t3 = s0.mul(&t0).sub(&s1);

        z5 = t3;

        let c0 = Fp6::new(z0, z1, z2, self.params.fp6_params.clone());

        let c1 = Fp6::new(z3, z4, z5, self.params.fp6_params.clone());

        Self::new(c0, c1, self.params.clone())
    }

    // -------------------------------------------------------------------------
    // Cyclotomic exponentiation.
    //
    // Only valid when self belongs to the cyclotomic subgroup.
    // -------------------------------------------------------------------------

    pub fn cyclotomic_exp(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());

        let mut started = false;

        for limb in exponent.iter().rev() {
            for bit in (0..64).rev() {
                if started {
                    result = result.cyclotomic_squared();
                }

                if ((limb >> bit) & 1) != 0 {
                    started = true;
                    result = result.mul(self);
                }
            }
        }

        result
    }

    // -------------------------------------------------------------------------
    // General exponentiation.
    // -------------------------------------------------------------------------

    pub fn pow(&self, exponent: &[u64]) -> Self {
        let mut result = Self::one(self.params.clone());
        let mut base = self.clone();

        for limb in exponent {
            let mut e = *limb;

            for _ in 0..64 {
                if (e & 1) != 0 {
                    result = result.mul(&base);
                }

                e >>= 1;
                base = base.squared();
            }
        }

        result
    }

    // -------------------------------------------------------------------------
    // Serialization.
    //
    // Fp6 = 3 * Fp2 = 24 u64 words.
    // Fp12 = 2 * Fp6 = 48 u64 words.
    // -------------------------------------------------------------------------

    pub fn to_words(&self) -> Vec<u64> {
        let mut words = self.c0.to_words();
        words.extend(self.c1.to_words());
        words
    }

    pub fn from_words(words: &[u64], params: Fp12Params) -> Option<Self> {
        if words.len() != 48 {
            return None;
        }

        let c0 = Fp6::from_words(&words[0..24], params.fp6_params.clone())?;

        let c1 = Fp6::from_words(&words[24..48], params.fp6_params.clone())?;

        Some(Self::new(c0, c1, params))
    }

    // -------------------------------------------------------------------------
    // Square root.
    //
    // Tonelli-Shanks in Fp12.
    //
    // p^12 is represented by 48 little-endian u64 limbs.
    // -------------------------------------------------------------------------

    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero(self.params.clone()));
        }

        let modulus = self.params.fp6_params.fp2_params.fp_params.modulus;

        // N = p^12
        let mut n = [0u64; 48];
        n[0] = 1;

        for _ in 0..12 {
            n = mul_48_by_4(&n, &modulus);
        }

        // Legendre exponent = (p^12 - 1) / 2
        let mut legendre_exp = n;
        sub_one_48(&mut legendre_exp);
        shr1_48(&mut legendre_exp);

        // N - 1 = q * 2^s
        let mut q = n;
        sub_one_48(&mut q);

        let mut s = 0usize;

        while !is_odd_48(&q) {
            shr1_48(&mut q);
            s += 1;
        }

        // x = a^((q+1)/2)
        let mut q_plus_one = q;
        add_one_48(&mut q_plus_one);
        shr1_48(&mut q_plus_one);

        let mut x = self.pow(&q_plus_one);
        let mut t = self.pow(&q);

        if t.is_one() {
            return Some(x);
        }

        // Find a quadratic non-residue.
        //
        // This bounded search is suitable as a reference implementation.
        // Production code should preferably receive a known non-residue
        // as part of Fp12Params.
        // Find a quadratic non-residue in Fp12.
        //
        // Elements embedded from Fp are always squares in Fp12,
        // because [Fp12 : Fp] = 12 is even.
        //
        // Therefore search elements with a nonzero Fp6 coefficient c1.
        // Legendre exponent = (p^12 - 1) / 2
        let mut legendre_exp = n;
        sub_one_48(&mut legendre_exp);
        shr1_48(&mut legendre_exp);

        let fp_params = self.params.fp6_params.fp2_params.fp_params.clone();

        let fp2_params = self.params.fp6_params.fp2_params.clone();

        let fp6_params = self.params.fp6_params.clone();

        let fp_zero = Fp::zero(fp_params.clone());
        let fp2_zero = Fp2::zero(fp2_params.clone());
        let fp6_zero = Fp6::zero(fp6_params.clone());

        let minus_one_fp = Fp::from_u64(modulus[0] - 1, fp_params.clone());

        let minus_one_fp2 = Fp2::new(minus_one_fp, fp_zero.clone(), fp2_params.clone());

        let minus_one_fp6 = Fp6::new(
            minus_one_fp2,
            fp2_zero.clone(),
            fp2_zero.clone(),
            fp6_params.clone(),
        );

        let minus_one_fp12 = Self::new(minus_one_fp6, fp6_zero.clone(), self.params.clone());

        let mut z: Option<Self> = None;

        for candidate in 0u64..4096 {
            let fp = Fp::from_u64(candidate, fp_params.clone());

            let fp2 = Fp2::new(fp, fp_zero.clone(), fp2_params.clone());

            let c0 = Fp6::new(fp2, fp2_zero.clone(), fp2_zero.clone(), fp6_params.clone());

            let candidate_z = Self::new(c0, Fp6::one(fp6_params.clone()), self.params.clone());

            if candidate_z.pow(&legendre_exp).equals(&minus_one_fp12) {
                z = Some(candidate_z);
                break;
            }
        }

        let z = z?;

        let mut c = z.pow(&q);
        let mut m = s;

        while !t.is_one() {
            let mut i = 1usize;
            let mut t2 = t.squared();

            while i < m && !t2.is_one() {
                t2 = t2.squared();
                i += 1;
            }

            if i == m {
                return None;
            }

            let mut exponent = [0u64; 48];
            exponent[0] = 1u64 << (m - i - 1);

            let b = c.pow(&exponent);
            let b2 = b.squared();

            x = x.mul(&b);
            t = t.mul(&b2);
            c = b2;
            m = i;
        }

        if x.squared().equals(self) {
            Some(x)
        } else {
            None
        }
    }
}

// =============================================================================
// Small Fp2 helpers
// =============================================================================

/// 3*a - 2*b
fn three_mul_sub_two(a: &Fp2, b: &Fp2) -> Fp2 {
    let t = a.sub(b);
    t.add(&t).add(a)
}

/// 3*a + 2*b
fn three_mul_add_two(a: &Fp2, b: &Fp2) -> Fp2 {
    let t = a.add(b);
    t.add(&t).add(a)
}

// =============================================================================
// 48-limb integer helpers for Tonelli-Shanks.
//
// All integers use little-endian u64 limbs.
// =============================================================================

fn is_odd_48(a: &[u64; 48]) -> bool {
    (a[0] & 1) != 0
}

fn shr1_48(a: &mut [u64; 48]) {
    let mut carry = 0u64;

    for i in (0..48).rev() {
        let next = a[i] & 1;
        a[i] = (a[i] >> 1) | (carry << 63);
        carry = next;
    }
}

fn add_one_48(a: &mut [u64; 48]) {
    for limb in a.iter_mut() {
        let (value, carry) = limb.overflowing_add(1);
        *limb = value;

        if !carry {
            break;
        }
    }
}

fn sub_one_48(a: &mut [u64; 48]) {
    for limb in a.iter_mut() {
        let (value, borrow) = limb.overflowing_sub(1);
        *limb = value;

        if !borrow {
            break;
        }
    }
}

/// Multiply a 48-limb number by a 4-limb number.
fn mul_48_by_4(a: &[u64; 48], b: &[u64; 4]) -> [u64; 48] {
    let mut out = [0u64; 48];

    for i in 0..48 {
        let mut carry = 0u128;

        for j in 0..4 {
            let k = i + j;

            if k >= 48 {
                break;
            }

            let value = (a[i] as u128) * (b[j] as u128) + (out[k] as u128) + carry;

            out[k] = value as u64;
            carry = value >> 64;
        }

        let mut k = i + 4;

        while carry != 0 && k < 48 {
            let value = (out[k] as u128) + carry;

            out[k] = value as u64;
            carry = value >> 64;
            k += 1;
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;

    use crate::algebra::fields::prime_base::fp::Fp;
    use crate::algebra::fields::prime_base::fp_montgomery::FieldParams;
    use crate::algebra::fields::prime_extension::fp2::{Fp2, Fp2Params};
    use crate::algebra::fields::prime_extension::fp6_3over2::{Fp6, Fp6Params};
    use crate::algebra::fields::prime_extension::fp12_2over3over2::Fp12Params;

    pub fn bn254_fp12_params() -> Fp12Params {
        let fp_params = FieldParams::new([
            0x3c20_8c16_d87c_fd47,
            0x9781_6a91_6871_ca8d,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ]);

        // Fp2: u^2 = -1
        let one = Fp::one(fp_params.clone());
        let minus_one = one.negate();

        let fp2_params = Fp2Params::new(
            fp_params.clone(),
            minus_one.clone(),
            [one.clone(), minus_one],
        );

        // xi = 9 + u
        let xi = Fp2::new(
            Fp::from_u64(9, fp_params.clone()),
            Fp::one(fp_params.clone()),
            fp2_params.clone(),
        );

        // Exponents are little-endian 64-bit limbs.
        let exp_q_minus_1_div_3 = [
            0x6960_2eb2_4829_a9c2,
            0xdd2b_2385_cd7b_4384,
            0xe81a_c1e7_8080_72c9,
            0x1021_6f7b_a065_e00d,
        ];

        let exp_2_q_minus_1_div_3 = [
            0xd2c0_5d64_9053_5384,
            0xba56_470b_9af6_8708,
            0xd035_83cf_0100_e593,
            0x2042_def7_40cb_c01b,
        ];

        let exp_q_minus_1_div_2 = [
            0x9e10_460b_6c3e_7ea3,
            0xcbc0_b548_b438_e546,
            0xdc28_22db_40c0_ac2e,
            0x1832_2739_7098_d014,
        ];

        // Fp6 Frobenius coefficients.
        let c1_first = xi.pow(&exp_q_minus_1_div_3);
        let c2_first = xi.pow(&exp_2_q_minus_1_div_3);

        let mut c1 = vec![Fp2::one(fp2_params.clone())];
        let mut c2 = vec![Fp2::one(fp2_params.clone())];

        c1.push(c1_first);
        c2.push(c2_first);

        for k in 1..5 {
            let next_c1 = c1[k].frobenius_map(1).mul(&c1[1]);

            let next_c2 = c2[k].frobenius_map(1).mul(&c2[1]);

            c1.push(next_c1);
            c2.push(next_c2);
        }

        let fp6_params = Fp6Params {
            fp2_params: fp2_params.clone(),
            non_residue: xi,
            frobenius_coeffs_c1: c1.try_into().expect("Fp6 c1 length"),
            frobenius_coeffs_c2: c2.try_into().expect("Fp6 c2 length"),
        };

        // Fp6 basis element w = (0, 1, 0).
        let zero2 = Fp2::zero(fp2_params.clone());
        let one2 = Fp2::one(fp2_params.clone());

        let w = Fp6::new(
            zero2.clone(),
            one2.clone(),
            zero2.clone(),
            fp6_params.clone(),
        );

        // Fp12 Frobenius coefficients.
        let d1 = w.pow(&exp_q_minus_1_div_2);

        let mut d = vec![Fp6::new(
            one2.clone(),
            zero2.clone(),
            zero2.clone(),
            fp6_params.clone(),
        )];

        d.push(d1);

        for k in 1..11 {
            let next = d[k].frobenius_map(1).mul(&d[1]);

            d.push(next);
        }

        Fp12Params {
            fp6_params,
            frobenius_coeffs_c1: d
                .try_into()
                .unwrap_or_else(|_| panic!("Fp12 coefficient length")),
        }
    }
    /// Test tower:
    ///
    /// Fp  = F7
    /// Fp2 = Fp[u] / (u^2 + 1)
    /// Fp6 = Fp2[w] / (w^3 - xi)
    /// Fp12 = Fp6[v] / (v^2 - w)
    ///
    /// xi = 1 + 2u
    ///
    /// Norm(xi) = 5, which is neither a square nor a cube in F7.
    fn params() -> Fp12Params {
        let fp_params = FieldParams::new([7, 0, 0, 0]);

        let one = Fp::one(fp_params.clone());
        let minus_one = Fp::from_u64(6, fp_params.clone());

        // Fp2: u^2 = -1
        let fp2_params = Fp2Params::new(
            fp_params.clone(),
            minus_one.clone(),
            [one.clone(), minus_one],
        );

        // xi = 1 + 2u
        let xi = Fp2::new(
            one.clone(),
            Fp::from_u64(2, fp_params.clone()),
            fp2_params.clone(),
        );

        // Fp6 Frobenius coefficients.
        //
        // c1[k] = xi^((p^k - 1) / 3)
        // c2[k] = xi^(2 * (p^k - 1) / 3)
        let frobenius_coeffs_c1 = std::array::from_fn(|k| {
            let p_k = 7u64.pow(k as u32);
            xi.pow(&[(p_k - 1) / 3])
        });

        let frobenius_coeffs_c2 = std::array::from_fn(|k| {
            let p_k = 7u64.pow(k as u32);
            xi.pow(&[2 * ((p_k - 1) / 3)])
        });

        let fp6_params = Fp6Params {
            fp2_params: fp2_params.clone(),
            non_residue: xi,
            frobenius_coeffs_c1,
            frobenius_coeffs_c2,
        };

        // Fp6 basis element w = (0, 1, 0).
        let w = Fp6::new(
            Fp2::zero(fp2_params.clone()),
            Fp2::one(fp2_params.clone()),
            Fp2::zero(fp2_params.clone()),
            fp6_params.clone(),
        );

        // Fp12 Frobenius coefficients:
        //
        // coeff[k] = w^((p^k - 1) / 2)
        //
        // Fp12 = Fp6[v] / (v^2 - w).
        let frobenius_coeffs_c1 = std::array::from_fn(|k| {
            let p_k = 7u64.pow(k as u32);
            let exponent = (p_k - 1) / 2;

            w.pow(&[exponent])
        });

        Fp12Params {
            fp6_params,
            frobenius_coeffs_c1,
        }
    }

    fn sample_element(params: Fp12Params) -> Fp12 {
        let fp2_params = params.fp6_params.fp2_params.clone();
        let fp6_params = params.fp6_params.clone();

        let fp = |x| Fp::from_u64(x, fp2_params.fp_params.clone());

        let fp2 = |a, b| Fp2::new(fp(a), fp(b), fp2_params.clone());

        let c0 = Fp6::new(fp2(1, 2), fp2(3, 4), fp2(5, 6), fp6_params.clone());

        let c1 = Fp6::new(fp2(2, 1), fp2(4, 3), fp2(6, 5), fp6_params);

        Fp12::new(c0, c1, params)
    }

    fn sample_nonzero_element(params: Fp12Params) -> Fp12 {
        let element = sample_element(params);

        assert!(!element.is_zero());

        element
    }

    #[test]
    fn test_fp12_add_sub() {
        let params = params();
        let a = sample_element(params.clone());
        let b = sample_element(params);

        assert!(a.add(&b).sub(&b).equals(&a));
    }

    #[test]
    fn test_fp12_mul_square() {
        let params = params();
        let a = sample_element(params);

        assert!(a.mul(&a).equals(&a.squared()));
    }

    #[test]
    fn test_fp12_two_square_algorithms() {
        let params = params();
        let a = sample_element(params);

        assert!(a.squared_complex().equals(&a.squared_karatsuba()));
    }

    #[test]
    fn test_fp12_inverse() {
        let params = params();
        let a = sample_nonzero_element(params.clone());

        let inv = a.inverse().expect("inverse exists");
        let one = Fp12::one(params);

        assert!(a.mul(&inv).equals(&one));
    }

    #[test]
    fn test_fp12_frobenius_period() {
        let params = params();
        let a = sample_element(params);

        assert!(a.frobenius_map(12).equals(&a));
    }

    #[test]
    fn test_fp12_serialization() {
        let params = params();
        let a = sample_element(params.clone());

        let words = a.to_words();

        assert_eq!(words.len(), 48);

        let decoded = Fp12::from_words(&words, params).unwrap();

        assert!(decoded.equals(&a));
    }

    #[test]
    fn test_fp12_sqrt() {
        let params = params();
        let a = sample_element(params);

        let square = a.squared();
        let root = square.sqrt().expect("square has a root");

        assert!(root.squared().equals(&square));
    }
    #[test]
    fn test_fp12_struct_sizes() {
        println!("Fp size: {}", std::mem::size_of::<Fp>());
        println!("Fp2Params size: {}", std::mem::size_of::<Fp2Params>());
        println!("Fp2 size: {}", std::mem::size_of::<Fp2>());
        println!("Fp6Params size: {}", std::mem::size_of::<Fp6Params>());
        println!("Fp6 size: {}", std::mem::size_of::<Fp6>());
        println!("Fp12Params size: {}", std::mem::size_of::<Fp12Params>());
        println!("Fp12 size: {}", std::mem::size_of::<Fp12>());
    }

    #[test]
    fn test_bn254_fp6_frobenius_cycle() {
        let params = bn254_fp12_params();
        let fp6_params = params.fp6_params;

        let x = Fp6::new(
            Fp2::from_u64(3, fp6_params.fp2_params.clone()),
            Fp2::from_u64(5, fp6_params.fp2_params.clone()),
            Fp2::from_u64(7, fp6_params.fp2_params.clone()),
            fp6_params,
        );

        assert!(x.frobenius_map(6).equals(&x));
    }
    #[test]
    fn test_bn254_fq12_frobenius_cycle() {
        use crate::algebra::curves::alt_bn128::alt_bn128_init::bn254_fp12_params;

        let params = bn254_fp12_params();
        let fp6_params = params.fp6_params.clone();

        let fp2_params = fp6_params.fp2_params.clone();

        let c0 = Fp6::new(
            Fp2::from_u64(3, fp2_params.clone()),
            Fp2::from_u64(5, fp2_params.clone()),
            Fp2::from_u64(7, fp2_params.clone()),
            fp6_params.clone(),
        );

        let c1 = Fp6::new(
            Fp2::from_u64(11, fp2_params.clone()),
            Fp2::from_u64(13, fp2_params.clone()),
            Fp2::from_u64(17, fp2_params),
            fp6_params,
        );

        let x = Fp12::new(c0, c1, params);

        assert!(
            x.frobenius_map(12).equals(&x),
            "BN254 Fp12 Frobenius^12 must be identity"
        );
    }
}
