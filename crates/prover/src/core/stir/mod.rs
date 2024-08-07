//! This module contains the core STIR implementation adapted from Python PoC made by Nethermind.

use blake2::{Blake2s256, Digest};
use num_bigint::{BigInt, Sign};
use num_traits::{Euclid, ToPrimitive};
use gaussian::Gaussian;

mod circle_stir;
mod fft;
mod circle_fft;
mod merkle_trees;
mod poly_utils;
mod gaussian;

fn mul_mod(x: &[i64], modulus: u32) -> i64 {
    x.iter().fold(1, |acc, &y| (acc * y).rem_euclid(modulus as i64))
}

fn to_32_be_bytes(x: i64) -> [u8; 32] {
    let mut res = [0; 32];
    res[24..].copy_from_slice(&x.to_be_bytes());
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
    assert!(modulus < 2_u32.pow(24)); // inherited from Vitalik's code, not sure if this is needed.
    let mut ans = Vec::new();
    for c in start..(count + start) {
        exclude.sort();
        let exclude2: Vec<usize> = exclude.iter().enumerate().map(|(i, &x)| x - i).collect();
        let val = BigInt::from_bytes_be(
            Sign::Plus,
            &blake(&[&to_32_be_bytes(c as i64), seed].concat()),
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

pub fn inv(a: i64, modulus: u32) -> i64 {
    let modulus = modulus as i64;
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

pub fn get_power_cycle(r: Gaussian, modulus: u32, offset: Gaussian) -> Vec<Gaussian> {
    let mut o = vec![offset];
    loop {
        let next = (o[o.len() - 1] * r).rem_euclid(modulus);
        if next == offset {
            break;
        }
        o.push(next);
    }

    o
}

pub trait PowMod {
    fn pow_mod(self, exp: u32, modulus: u32) -> Self;
}

impl PowMod for i64 {
    fn pow_mod(self, exp: u32, modulus: u32) -> i64 {
        BigInt::from(self).modpow(&BigInt::from(exp), &BigInt::from(modulus)).to_i64().unwrap()
    }
}

pub trait RemEuclid {
    fn rem_euclid(self, modulus: u32) -> Self;
}

impl RemEuclid for i64 {
    fn rem_euclid(self, modulus: u32) -> i64 {
        Euclid::rem_euclid(&self, &(modulus as i64))
    }
}
