//! This module contains the core STIR implementation adapted from Python PoC made by Nethermind.

use blake2::{Blake2s256, Digest};
use num_bigint::{BigInt, Sign};
use num_traits::{Euclid, ToPrimitive};

mod stir;
mod fft;
mod merkle_trees;
mod poly_utils;

#[derive(Debug, Eq, PartialEq, Clone)]
struct Point {
    x: i64,
    y: i64,
}

fn pow_mod(base: i64, exp: u32, modulus: u32) -> i64 {
    BigInt::from(base).modpow(&BigInt::from(exp), &BigInt::from(modulus)).to_i64().unwrap()
}

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
    start: usize /* default=0 */,
) -> Vec<usize> {
    let mut exclude = vec![];
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
        let val = val + exclude2.binary_search(&val).unwrap_or_else(|x| x) as usize;
        ans.push(val);
        exclude.push(val);
    }
    ans
}
