use std::ops::Add;

/// 提供大整数的 bit 访问能力。
pub trait ScalarBits {
    /// 标量可访问的总 bit 数。
    fn max_bits(&self) -> usize;

    /// 判断指定 bit 是否为 1。
    fn test_bit(&self, index: usize) -> bool;
}

/// 提供群元素的基本运算。
pub trait GroupElement: Sized + Add<Output = Self> {
    /// 返回群的单位元。
    fn zero() -> Self;

    /// 返回当前元素的两倍。
    fn dbl(&self) -> Self;
}

/// 计算 base * scalar。
///
/// 使用从最高有效位到最低有效位的
/// left-to-right double-and-add 算法。
pub fn scalar_mul<GroupT, ScalarT>(
    base: &GroupT,
    scalar: &ScalarT,
) -> GroupT
where
    GroupT: GroupElement,
    ScalarT: ScalarBits,
{
    let mut result = GroupT::zero();
    let mut found_one = false;

    for i in (0..scalar.max_bits()).rev() {
        if found_one {
            result = result.dbl();
        }

        if scalar.test_bit(i) {
            found_one = true;
            result = result + base.clone();
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ops::Add;

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct TestGroup(u64);

    impl Add for TestGroup {
        type Output = Self;

        fn add(self, rhs: Self) -> Self {
            Self(self.0 + rhs.0)
        }
    }

    impl GroupElement for TestGroup {
        fn zero() -> Self {
            Self(0)
        }

        fn dbl(&self) -> Self {
            Self(self.0 * 2)
        }
    }

    struct TestScalar(u64);

    impl ScalarBits for TestScalar {
        fn max_bits(&self) -> usize {
            64
        }

        fn test_bit(&self, index: usize) -> bool {
            ((self.0 >> index) & 1) == 1
        }
    }

    #[test]
    fn test_scalar_mul_zero() {
        let base = TestGroup(5);
        let scalar = TestScalar(0);

        assert_eq!(scalar_mul(&base, &scalar), TestGroup(0));
    }

    #[test]
    fn test_scalar_mul_one() {
        let base = TestGroup(5);
        let scalar = TestScalar(1);

        assert_eq!(scalar_mul(&base, &scalar), TestGroup(5));
    }

    #[test]
    fn test_scalar_mul_13() {
        let base = TestGroup(5);
        let scalar = TestScalar(13);

        assert_eq!(scalar_mul(&base, &scalar), TestGroup(65));
    }

    #[test]
    fn test_scalar_mul_power_of_two() {
        let base = TestGroup(3);
        let scalar = TestScalar(8);

        assert_eq!(scalar_mul(&base, &scalar), TestGroup(24));
    }
}