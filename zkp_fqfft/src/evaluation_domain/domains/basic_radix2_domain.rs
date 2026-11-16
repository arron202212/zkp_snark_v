use crate::polynomial_arithmetic::basic_operations::PolynomialError;
use crate::polynomial_arithmetic::field::Radix2Field;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DomainError {
    InvalidSize,
    DomainTooLarge,
    InvalidInputLength { expected: usize, actual: usize },
    InvalidCoset,
    Polynomial(PolynomialError),
}

#[derive(Clone)]
pub struct BasicRadix2Domain<F: Radix2Field> {
    pub m: usize,
    pub omega: F,
    pub field: F,
}

impl<F: Radix2Field> BasicRadix2Domain<F> {
    pub fn new(m: usize, field: F) -> Result<Self, DomainError> {
        if m <= 1 || !m.is_power_of_two() {
            return Err(DomainError::InvalidSize);
        }

        let log_m = m.trailing_zeros();

        if log_m > field.two_adicity() {
            return Err(DomainError::DomainTooLarge);
        }

        let omega = field
            .root_of_unity(log_m)
            .ok_or(DomainError::DomainTooLarge)?;

        let one = omega.one_like();

        if omega.pow_usize(m) != one {
            return Err(DomainError::DomainTooLarge);
        }

        if omega.pow_usize(m / 2) == one {
            return Err(DomainError::DomainTooLarge);
        }

        Ok(Self { m, omega, field })
    }

    pub fn size(&self) -> usize {
        self.m
    }

    pub fn get_omega(&self) -> &F {
        &self.omega
    }

    pub fn get_domain_element(&self, index: usize) -> F {
        self.omega.pow_usize(index % self.m)
    }

    fn one(&self) -> F {
        self.field.one_like()
    }

    fn zero(&self) -> F {
        self.field.sub(&self.field)
    }

    fn check_length(&self, values: &[F]) -> Result<(), DomainError> {
        if values.len() != self.m {
            return Err(DomainError::InvalidInputLength {
                expected: self.m,
                actual: values.len(),
            });
        }

        Ok(())
    }

    /// Iterative radix-2 FFT.
    ///
    /// Input and output are in natural order.
    pub fn fft(&self, values: &mut [F]) -> Result<(), DomainError> {
        self.check_length(values)?;

        let n = self.m;
        let bits = n.trailing_zeros();

        // Bit-reversal permutation.
        for i in 0..n {
            let j = i.reverse_bits() >> (usize::BITS - bits);

            if i < j {
                values.swap(i, j);
            }
        }

        // Iterative butterfly stages.
        let mut len = 2usize;

        loop {
            let half = len / 2;
            let step = self.omega.pow_usize(n / len);

            for start in (0..n).step_by(len) {
                let mut w = self.one();

                for j in 0..half {
                    let u = values[start + j].clone();
                    let v = values[start + j + half].mul(&w);

                    values[start + j] = u.add(&v);
                    values[start + j + half] = u.sub(&v);

                    w = w.mul(&step);
                }
            }

            if len == n {
                break;
            }

            len *= 2;
        }

        Ok(())
    }

    /// Inverse radix-2 FFT.
    pub fn inverse_fft(&self, values: &mut [F]) -> Result<(), DomainError> {
        self.check_length(values)?;

        let omega_inv = self.omega.inverse().ok_or(DomainError::DomainTooLarge)?;

        let inverse_domain = Self {
            m: self.m,
            omega: omega_inv,
            field: self.field.clone(),
        };

        inverse_domain.fft(values)?;

        let n = self
            .field
            .from_u64_like(u64::try_from(self.m).map_err(|_| DomainError::DomainTooLarge)?);

        let n_inv = n.inverse().ok_or(DomainError::DomainTooLarge)?;

        for value in values.iter_mut() {
            *value = value.mul(&n_inv);
        }

        Ok(())
    }

    /// Coset FFT:
    /// FFT(a[i] * g^i).
    pub fn coset_fft(&self, values: &mut [F]) -> Result<(), DomainError> {
        self.check_length(values)?;

        let g = self.field.multiplicative_generator();
        let mut power = self.one();

        for value in values.iter_mut() {
            *value = value.mul(&power);
            power = power.mul(&g);
        }

        self.fft(values)
    }

    /// Inverse coset FFT:
    /// iFFT followed by multiplication by g^(-i).
    pub fn icoset_fft(&self, values: &mut [F]) -> Result<(), DomainError> {
        self.check_length(values)?;

        self.inverse_fft(values)?;

        let g_inv = self
            .field
            .multiplicative_generator()
            .inverse()
            .ok_or(DomainError::InvalidCoset)?;

        let mut power = self.one();

        for value in values.iter_mut() {
            *value = value.mul(&power);
            power = power.mul(&g_inv);
        }

        Ok(())
    }

    /// Evaluate all Lagrange basis polynomials at t.
    ///
    /// Returns:
    /// [L_0(t), L_1(t), ..., L_(m-1)(t)]
    pub fn evaluate_all_lagrange_polynomials(&self, t: &F) -> Vec<F> {
        let one = self.one();

        let z = t.pow_usize(self.m).sub(&one);

        // Special case: t is already a domain element.
        if z.is_zero() {
            let mut result = vec![self.zero(); self.m];

            for i in 0..self.m {
                if *t == self.get_domain_element(i) {
                    result[i] = one;
                    return result;
                }
            }

            // Mathematically unreachable for a valid field/domain.
            return result;
        }

        // l = (t^m - 1) / m.
        let m_element = self.field.from_u64_like(self.m as u64);
        let m_inv = m_element
            .inverse()
            .expect("domain size must be nonzero in the field");

        let mut l = z.mul(&m_inv);
        let mut omega_i = one.clone();

        let mut result = Vec::with_capacity(self.m);

        for _ in 0..self.m {
            let denominator = t.sub(&omega_i);

            let denominator_inv = denominator.inverse().expect("t is not a domain element");

            result.push(l.mul(&denominator_inv));

            l = l.mul(&self.omega);
            omega_i = omega_i.mul(&self.omega);
        }

        result
    }

    /// Vanishing polynomial Z(t) = t^m - 1.
    pub fn vanishing_polynomial(&self, t: &F) -> F {
        t.pow_usize(self.m).sub(&self.one())
    }

    /// Add coeff * (X^m - 1) to a polynomial.
    ///
    /// Coefficients are stored in ascending degree order.
    ///
    /// The vector is extended when necessary.
    pub fn add_poly_Z(&self, coefficients: &mut Vec<F>, coeff: &F) {
        if coefficients.len() <= self.m {
            coefficients.resize(self.m + 1, self.zero());
        }

        coefficients[self.m] = coefficients[self.m].add(coeff);
        coefficients[0] = coefficients[0].sub(coeff);
    }

    /// Divide a polynomial by Z on the multiplicative coset.
    ///
    /// On g * H, the vanishing polynomial is constant:
    /// Z(g * omega^i) = g^m - 1.
    ///
    /// Therefore:
    /// coset FFT -> multiply by (g^m - 1)^(-1) -> inverse coset FFT.
    pub fn divide_by_Z_on_coset(&self, coefficients: &mut [F]) -> Result<(), DomainError> {
        self.check_length(coefficients)?;

        let g = self.field.multiplicative_generator();

        let z_on_coset = g.pow_usize(self.m).sub(&self.one());

        let z_inv = z_on_coset.inverse().ok_or(DomainError::InvalidCoset)?;

        self.coset_fft(coefficients)?;

        for value in coefficients.iter_mut() {
            *value = value.mul(&z_inv);
        }

        self.icoset_fft(coefficients)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial_arithmetic::field::PolynomialField;
    use zkp_ff::algebra::fields::prime_base::fp::Fp;
    use zkp_ff::algebra::fields::prime_base::fp_montgomery::FieldParams;
    fn fr_params() -> FieldParams {
        FieldParams::new([
            0x43e1_f593_f000_0001,
            0x2833_e848_79b9_7091,
            0xb850_45b6_8181_585d,
            0x3064_4e72_e131_a029,
        ])
    }

    fn field() -> Fp {
        Fp::one(fr_params())
    }

    fn f(value: u64) -> Fp {
        Fp::from_u64(value, fr_params())
    }

    #[test]
    fn test_domain_creation() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        assert_eq!(domain.size(), 8);
        assert_eq!(domain.omega.pow_usize(8), f(1));
        assert_ne!(domain.omega.pow_usize(4), f(1));
    }

    #[test]
    fn test_invalid_domain_size() {
        assert!(matches!(
            BasicRadix2Domain::new(0, field()),
            Err(DomainError::InvalidSize)
        ));

        assert!(matches!(
            BasicRadix2Domain::new(1, field()),
            Err(DomainError::InvalidSize)
        ));

        assert!(matches!(
            BasicRadix2Domain::new(7, field()),
            Err(DomainError::InvalidSize)
        ));
    }

    #[test]
    fn test_fft_inverse_fft() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let original = vec![f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8)];

        let mut values = original.clone();

        domain.fft(&mut values).unwrap();
        domain.inverse_fft(&mut values).unwrap();

        assert_eq!(values, original);
    }

    #[test]
    fn test_fft_zero_polynomial() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let mut values = vec![f(0); 8];

        domain.fft(&mut values).unwrap();

        assert_eq!(values, vec![f(0); 8]);
    }

    #[test]
    fn test_fft_constant_polynomial() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let mut values = vec![f(9), f(0), f(0), f(0), f(0), f(0), f(0), f(0)];

        domain.fft(&mut values).unwrap();

        assert_eq!(values, vec![f(9); 8]);
    }

    #[test]
    fn test_coset_fft_inverse() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let original = vec![f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8)];

        let mut values = original.clone();

        domain.coset_fft(&mut values).unwrap();
        domain.icoset_fft(&mut values).unwrap();

        assert_eq!(values, original);
    }

    #[test]
    fn test_domain_elements_are_distinct() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        for i in 0..8 {
            for j in (i + 1)..8 {
                assert_ne!(domain.get_domain_element(i), domain.get_domain_element(j));
            }
        }
    }

    #[test]
    fn test_vanishing_polynomial_on_domain() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        for i in 0..8 {
            let x = domain.get_domain_element(i);

            assert!(domain.vanishing_polynomial(&x).is_zero());
        }
    }

    #[test]
    fn test_lagrange_at_domain_elements() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        for i in 0..8 {
            let x = domain.get_domain_element(i);

            let values = domain.evaluate_all_lagrange_polynomials(&x);

            assert_eq!(values.len(), 8);

            for j in 0..8 {
                if i == j {
                    assert_eq!(values[j], f(1));
                } else {
                    assert!(values[j].is_zero());
                }
            }
        }
    }

    #[test]
    fn test_lagrange_partition_of_unity() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let t = f(123);

        let values = domain.evaluate_all_lagrange_polynomials(&t);

        let sum = values.iter().fold(f(0), |acc, x| acc.add(x));

        assert_eq!(sum, f(1));
    }

    #[test]
    fn test_add_poly_Z() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let mut coefficients = vec![f(1), f(2), f(3)];
        let coeff = f(5);

        domain.add_poly_Z(&mut coefficients, &coeff);

        assert_eq!(coefficients.len(), 9);
        assert_eq!(coefficients[0], f(1).sub(&f(5)));
        assert_eq!(coefficients[8], f(5));
    }

    #[test]
    fn test_divide_by_Z_on_coset_roundtrip() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let original = vec![f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8)];

        let g = field().multiplicative_generator();
        let z = g.pow_usize(8).sub(&f(1));

        let mut multiplied = original.clone();

        for value in multiplied.iter_mut() {
            *value = value.mul(&z);
        }

        domain.divide_by_Z_on_coset(&mut multiplied).unwrap();

        assert_eq!(multiplied, original);
    }

    #[test]
    fn test_wrong_input_length() {
        let domain = BasicRadix2Domain::new(8, field()).unwrap();

        let mut values = vec![f(1); 4];

        assert!(matches!(
            domain.fft(&mut values),
            Err(DomainError::InvalidInputLength {
                expected: 8,
                actual: 4
            })
        ));
    }
}
