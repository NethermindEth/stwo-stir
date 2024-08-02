//! This module contains the core STIR implementation adapted from Python PoC made by Nethermind.

use num_bigint::BigInt;
use num_traits::ToPrimitive;

mod fft;
mod merkle_trees;
mod poly_utils;

#[derive(Debug, Eq, PartialEq, Clone)]
struct Point {
    x: i64,
    y: i64,
}

fn pow_mod(base: u64, exp: u32, modulus: u64) -> u64 {
    BigInt::from(base).modpow(&BigInt::from(exp), &BigInt::from(modulus)).to_u64().unwrap()
}
