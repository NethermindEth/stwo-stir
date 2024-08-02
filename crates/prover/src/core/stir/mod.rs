//! This module contains the core STIR implementation adapted from Python PoC made by Nethermind.

mod merkle_trees;

fn pow_mod(base: u64, exp: u32, modulus: u64) -> u64 {
    // FIXME
    base.pow(exp) % modulus
}
