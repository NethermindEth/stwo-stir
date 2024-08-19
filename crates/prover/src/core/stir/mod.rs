//! This module contains the core STIR implementation adapted from Python PoC made by Nethermind.

use blake2::{Blake2s256, Digest};
use num_bigint::{BigInt, Sign};
use num_traits::{Euclid, NumOps, One, Pow, ToPrimitive};
use crate::core::fields::cm31::CM31;
use crate::core::fields::{m31, ComplexConjugate};
use crate::core::fields::m31::M31;

#[allow(dead_code)]

mod stir;
mod fft;
mod merkle_trees;
mod poly_utils;

trait Xy {
    fn x(&self) -> i128;
    fn y(&self) -> i128;
}

trait KindaField: NumOps<Self> + PartialEq + Xy + Pow<u32, Output = Self> + One + Copy {}

impl<F> KindaField for F where F: NumOps<F> + PartialEq + Xy + Pow<u32, Output = Self> + One + Copy {}

fn to_32_be_bytes(x: i128) -> [u8; 32] {
    let mut res = [0; 32];
    res[16..].copy_from_slice(&x.to_be_bytes());
    res
}

fn blake(x: &[u8]) -> Vec<u8> {
    Blake2s256::digest(x).to_vec()
}

/// Extract pseudorandom indices from entropy
fn get_pseudorandom_indices(
    seed: &[u8],
    modulus: u32,
    count: usize,
    start: usize,
    exclude: &mut Vec<usize>,
) -> Vec<usize> {
    let mut ans = Vec::new();
    for c in start..(count + start) {
        exclude.sort();
        let exclude2: Vec<usize> = exclude.iter().enumerate().map(|(i, &x)| x - i).collect();
        let val = BigInt::from_bytes_be(
            Sign::Plus,
            &blake(&[&to_32_be_bytes(c as i128), seed].concat()),
        ).rem_euclid(
            &BigInt::from(modulus - exclude.len() as u32)
        ).to_usize().unwrap();
        let bisection_point = match exclude2.binary_search(&val) {
            Ok(mut i) => {
                // We need a last index, not the first one
                while i < exclude2.len() && exclude[i] == val {
                    i += 1;
                }
                i
            }
            Err(x) => x,
        };
        let val = val + bisection_point;
        ans.push(val);
        exclude.push(val);
    }
    ans
}

fn inv(a: i128, modulus: u32) -> i128 {
    let modulus = modulus as i128;
    let (mut lm, mut hm) = (1, 0);
    let (mut low, mut high) = (a.rem_euclid(modulus), modulus);
    if low == 0 {
        panic!("ZeroDivisionError");
    }
    while low > 1 {
        let r = high / low;
        let (nm, new) = (hm - lm * r, high - low * r);
        (lm, low, hm, high) = (nm, new, lm, low);
    }
    lm.rem_euclid(modulus)
}

fn get_power_cycle<F: KindaField>(r: F, offset: F) -> Vec<F> {
    let mut o = vec![offset];
    loop {
        let next = o[o.len() - 1] * r;
        if next == offset {
            break;
        }
        o.push(next);
    }

    o
}

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

impl Xy for CM31 {
    fn x(&self) -> i128 { self.0.0 as i128 }
    fn y(&self) -> i128 { self.1.0 as i128 }
}
