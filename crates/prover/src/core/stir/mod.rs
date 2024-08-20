//! This module contains the core STIR implementation adapted from Python PoC made by Nethermind.

use blake2::{Blake2s256, Digest};
use num_bigint::{BigInt, Sign};
use num_traits::{Euclid, ToPrimitive};
use crate::core::fields::cm31::CM31;
use crate::core::fields::{m31, Field};
use crate::core::fields::m31::M31;

#[allow(dead_code)]

mod stir;
mod fft;
mod merkle_trees;
mod poly_utils;

const MODULUS: u32 = m31::P;

trait Xy<F> {
    fn x(&self) -> F;
    fn y(&self) -> F;
}

trait StirField<B>: Field + Xy<B> {}

impl<B, F> StirField<B> for F where F: Field + Xy<B> {}

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

fn get_power_cycle<B, F: StirField<B>>(r: F, offset: F) -> Vec<F> {
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

impl Xy<M31> for CM31 {
    fn x(&self) -> M31 { self.0 }
    fn y(&self) -> M31 { self.1 }
}
