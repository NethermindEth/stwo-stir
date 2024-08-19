use std::ops::*;
use num_traits::{ConstOne, ConstZero, One, Pow, Zero};
use super::gaussian::Gaussian;
use super::{PowMod, RemEuclid, Xy};

/// A gaussian (complex) number, x is the real part, y is imaginary.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct GaussianF<const MOD: u64> {
    pub x: u64,
    pub y: u64,
}

impl<const MOD: u64> GaussianF<MOD> {
    pub fn new_signed(x: i64, y: i64) -> Self {
        Self { x: x.rem_euclid(MOD as i64) as u64, y: y.rem_euclid(MOD as i64) as u64 }
    }

    pub fn new(x: u64, y: u64) -> Self {
        Self { x: x.rem_euclid(MOD), y: y.rem_euclid(MOD) }
    }

    pub fn conj(self, parity: u64 /* default = 1 */) -> Self {
        if parity % 2 == 0 {
            self
        } else {
            GaussianF::new_signed(self.x as i64, -(self.y as i64))
        }
    }
}

impl<const MOD: u64> Zero for GaussianF<MOD> {
    fn zero() -> Self { Self::ZERO }
    fn is_zero(&self) -> bool { *self == Self::ZERO }
}

impl<const MOD: u64> One for GaussianF<MOD> {
    fn one() -> Self { Self::ONE }
}

impl<const MOD: u64> ConstZero for GaussianF<MOD> {
    const ZERO: GaussianF<MOD> = GaussianF { x: 0, y: 0 };
}

impl<const MOD: u64> ConstOne for GaussianF<MOD> {
    const ONE: GaussianF<MOD> = GaussianF { x: 1, y: 0 };
}

impl<const MOD: u64> Add for GaussianF<MOD> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(
            self.x + rhs.x,
            self.y + rhs.y,
        )
    }
}

impl<const MOD: u64> Sub for GaussianF<MOD> {
    type Output = Self;

    fn sub(self, rhs: GaussianF<MOD>) -> Self::Output {
        Self::new(
            self.x - rhs.x,
            self.y - rhs.y,
        )
    }
}

impl<const MOD: u64> Mul for GaussianF<MOD> {
    type Output = Self;

    fn mul(self, rhs: GaussianF<MOD>) -> Self::Output {
        let (x1, y1) = (self.x as i128, self.y as i128);
        let (x2, y2) = (rhs.x as i128, rhs.y as i128);
        Self::new(
            (x1 * x2 - y1 * y2).rem_euclid(MOD as i128) as u64,
            (x1 * y2 + y1 * x2).rem_euclid(MOD as i128) as u64,
        )
    }
}

impl<const MOD: u64> Div for GaussianF<MOD> {
    type Output = Self;

    fn div(self, _rhs: GaussianF<MOD>) -> Self::Output {
        unimplemented!("Division not implemented for GaussianF")
    }
}

impl<const MOD: u64> Rem for GaussianF<MOD> {
    type Output = Self;

    fn rem(self, _rhs: GaussianF<MOD>) -> Self::Output {
        unimplemented!("Remainder not implemented for GaussianF")
    }
}

impl<const MOD: u64> Pow<u32> for GaussianF<MOD> {
    type Output = Self;

    fn pow(self, exp: u32) -> Self::Output {
        let mut ans = GaussianF::ONE;
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

impl<const MOD: u64> PowMod for GaussianF<MOD> {
    fn pow_mod(self, exp: u32, modulus: u32) -> Self {
        assert_eq!(MOD, modulus as u64);
        self.pow(exp)
    }
}

impl<const MOD: u64> RemEuclid for GaussianF<MOD> {
    fn rem_euclid(self, modulus: u32) -> Self {
        // No need to do anything
        assert_eq!(MOD, modulus as u64);
        self
    }
}

impl<const MOD: u64> Xy for GaussianF<MOD> {
    fn x(&self) -> i128 { self.x as i128 }
    fn y(&self) -> i128 { self.y as i128 }
}

impl<const MOD: u64> From<Gaussian> for GaussianF<MOD> {
    fn from(g: Gaussian) -> GaussianF<MOD> {
        GaussianF::new_signed(g.x as i64, g.y as i64)
    }
}

impl<const MOD: u64> Into<Gaussian> for GaussianF<MOD> {
    fn into(self) -> Gaussian {
        Gaussian::new(self.x as i128, self.y as i128)
    }
}

