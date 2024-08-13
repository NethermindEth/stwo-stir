use std::ops::*;
use num_traits::{ConstOne, ConstZero, One, Zero};
use super::*;

/// A gaussian (complex) number, x is the real part, y is imaginary.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct Gaussian {
    pub x: i64,
    pub y: i64,
}

impl Gaussian {
    pub fn new(x: i64, y: i64) -> Self {
        Self { x, y }
    }

    pub fn conj(self, parity: u32 /* default = 1 */) -> Self {
        if parity % 2 == 0 {
            self
        } else {
            Gaussian { x: self.x, y: -self.y }
        }
    }
}

impl Zero for Gaussian {
    fn zero() -> Self { Self::ZERO }
    fn is_zero(&self) -> bool { *self == Self::ZERO }
}

impl One for Gaussian {
    fn one() -> Self { Self::ONE }
}

impl ConstZero for Gaussian {
    const ZERO: Gaussian = Gaussian { x: 0, y: 0 };
}

impl ConstOne for Gaussian {
    const ONE: Gaussian = Gaussian { x: 1, y: 0 };
}

impl Add for Gaussian {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Gaussian {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Mul<Gaussian> for Gaussian {
    type Output = Self;

    fn mul(self, rhs: Gaussian) -> Self::Output {
        Gaussian {
            x: self.x * rhs.x - self.y * rhs.y,
            y: self.x * rhs.y + self.y * rhs.x,
        }
    }
}

impl PowMod for Gaussian {
    fn pow_mod(self, exp: u32, modulus: u32) -> Self {
        let mut ans = Gaussian::new(1, 0);
        let mut v = self;
        let mut exp = exp;
        while exp != 0 {
            if exp % 2 == 1 {
                ans = (ans * v).rem_euclid(modulus);
            }
            v = (v * v).rem_euclid(modulus);
            exp = exp / 2;
        }
        ans
    }
}

impl RemEuclid for Gaussian {
    fn rem_euclid(self, modulus: u32) -> Self {
        Gaussian { x: self.x.rem_euclid(modulus as i64), y: self.y.rem_euclid(modulus as i64) }
    }
}
