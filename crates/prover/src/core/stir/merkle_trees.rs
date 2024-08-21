//! Translated from Nethermind STIR's merkle_trees.py

use itertools::Itertools;
use super::*;

pub fn merkelize(l: &[M31]) -> Vec<Vec<Hash>> {
    assert!(l.len().is_power_of_two());
    let ilog2 = l.len().ilog2() as usize;
    let mut layers: Vec<Vec<Hash>> = (0..ilog2)
        .map(|layer_num| vec![[0; 32]; 2_usize.pow(layer_num as u32)])
        .collect_vec();
    layers.push(l.iter().map(|&x| to_32_be_bytes(x.0 as u128)).collect_vec());
    for layer_num in (0_usize..ilog2).rev() {
        for i in 0..layers[layer_num].len() {
            let prev = &layers[layer_num + 1][(i * 2)..=(i * 2 + 1)];
            let concat = prev.iter().flat_map(|v| *v).collect_vec();
            let hashed = blake(&concat);
            layers[layer_num][i] = hashed;
        }
    }
    layers
}

pub fn mk_branch(tree: &[Vec<Hash>], index: usize) -> Vec<Hash> {
    let mut index = index;
    let mut o = Vec::with_capacity(tree.len());
    let leaf_layer = tree.last().unwrap();
    o.push(leaf_layer[index]);
    o.push(leaf_layer[index ^ 1]);
    // Non-root intermediate layers
    for layer in tree[1..(tree.len() - 1)].iter().rev() {
        index /= 2;
        o.push(layer[index ^ 1]);
    }
    o
}

#[must_use]
pub fn verify_branch(root: &Hash, index: usize, proof: &[Hash]) -> bool {
    let mut index = index + 2_usize.pow(proof.len() as u32);
    let mut v = proof[0];
    for &p in &proof[1..] {
        if index % 2 == 0 {
            v = blake_pair(&v, &p);
        } else {
            v = blake_pair(&p, &v);
        }
        index /= 2;
    }
    &v == root
}

fn blake_pair(hash_1: &Hash, hash_2: &Hash) -> Hash {
    blake(&[*hash_1, *hash_2].iter().flat_map(|v| *v).collect_vec())
}
