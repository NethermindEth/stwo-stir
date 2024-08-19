use num_traits::{One, Pow};
use crate::core::fields::cm31::CM31;
use crate::core::fields::{m31, ComplexConjugate};
use crate::core::fields::m31::M31;
use super::gaussian::Gaussian;
use super::{PowMod, RemEuclid, Xy};

/// A gaussian (complex) number, x is the real part, y is imaginary.
impl CM31 {
    pub fn new(x: u64, y: u64) -> Self {
        Self(M31(x.rem_euclid(m31::P as u64) as u32), M31(y.rem_euclid(m31::P as u64) as u32))
    }

    pub fn conj(self, parity: u64 /* default = 1 */) -> Self {
        if parity % 2 == 0 {
            self
        } else {
            self.complex_conjugate()
        }
    }
}

impl Pow<u32> for CM31 {
    type Output = Self;

    fn pow(self, exp: u32) -> Self::Output {
        let mut ans = Self::one();
        let mut v = self;
        let mut exp = exp;
        while exp != 0 {
            if exp % 2 == 1 {
                ans = ans * v;
            }
            v = v * v;
            exp = exp / 2;
        }
        ans
    }
}

impl PowMod for CM31 {
    fn pow_mod(self, exp: u32, modulus: u32) -> Self {
        assert_eq!(m31::P, modulus);
        self.pow(exp)
    }
}

impl RemEuclid for CM31 {
    fn rem_euclid(self, modulus: u32) -> Self {
        // No need to do anything
        assert_eq!(m31::P, modulus);
        self
    }
}

impl Xy for CM31 {
    fn x(&self) -> i128 { self.0.0 as i128 }
    fn y(&self) -> i128 { self.1.0 as i128 }
}

impl From<Gaussian> for CM31 {
    fn from(g: Gaussian) -> CM31 {
        CM31(M31(g.x.rem_euclid(m31::P as i128) as u32), M31(g.y.rem_euclid(m31::P as i128) as u32))
    }
}

impl Into<Gaussian> for CM31 {
    fn into(self) -> Gaussian {
        Gaussian::new(self.x(), self.y())
    }
}

