//! Translated from Nethermind STIR's merkle_trees.py

use super::*;

pub fn merkelize(l: &[i128]) -> Vec<Vec<u8>> {
    let mut nodes: Vec<Vec<u8>> = vec![vec![0; 32]; l.len()];
    nodes.extend(l.iter().map(|&x| to_32_be_bytes(x).to_vec()));
    for i in (1..l.len()).rev() {
        let x1 = &nodes[i * 2];
        let x2 = &nodes[i * 2 + 1];
        let x3 = concat_clone(x1, x2);
        let x4 = blake(&x3);
        nodes[i] = x4; //blake(&concat_clone(&nodes[i * 2], &nodes[i * 2 + 1]));
    }
    nodes
}

pub fn mk_branch(tree: &[Vec<u8>], index: usize) -> Vec<Vec<u8>> {
    let mut index = index + tree.len() / 2;
    let mut o = vec![tree[index].clone()];
    while index > 1 {
        o.push(tree[index ^ 1].clone());
        index /= 2;
    }
    o
}

pub fn verify_branch(root: &[u8], index: usize, proof: &[&[u8]]) -> Vec<u8> {
    let mut index = index + 2_usize.pow(proof.len() as u32);
    let mut v = proof[0].to_vec();
    for &p in &proof[1..] {
        if index % 2 == 0 {
            v = blake(&concat_clone(&v, p));
        } else {
            v = blake(&concat_clone(p, &v));
        }
        index /= 2;
    }
    assert_eq!(v, root);
    proof[0].to_vec()
}
/// Byte length of a multi proof
fn _bin_length(proof: &[Vec<Vec<u8>>]) -> usize {
    proof.iter().map(|x| {
        let concat_len = x.iter().map(|y| y.len()).sum::<usize>();
        concat_len + x.len() / 8
    }).sum::<usize>() + proof.len() * 2
}

fn concat_clone<T: Clone>(slice_1: &[T], slice_2: &[T]) -> Vec<T> {
    let mut v = Vec::with_capacity(slice_1.len() + slice_2.len());
    v.extend_from_slice(slice_1);
    v.extend_from_slice(slice_2);
    v
}
