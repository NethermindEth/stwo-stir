use std::collections::HashMap;
use blake2::{Blake2s256, Digest};
use itertools::Itertools;

pub fn merkelize(l: &[u64]) -> Vec<Vec<u8>> {
    let mut nodes: Vec<Vec<u8>> = vec![vec![]; l.len()];
    nodes.extend(l.iter().map(|&x| {
        let mut res = vec![0; 32];
        res[24..].copy_from_slice(&x.to_be_bytes());
        res
    }));
    for i in (1..l.len()).rev() {
        nodes[i] = blake(&concat_clone(&nodes[i * 2], &nodes[i * 2 + 1]));
    }
    nodes
}

fn mk_branch(tree: &[Vec<u8>], index: usize) -> Vec<Vec<u8>> {
    let mut index = index + tree.len() / 2;
    let mut o = vec![tree[index].clone()];
    while index > 1 {
        o.push(tree[index ^ 1].clone());
        index /= 2;
    }
    o
}

fn verify_branch(root: &[u8], index: usize, proof: &[Vec<u8>]) -> Vec<u8> {
    let mut index = index + 2_usize.pow(proof.len() as u32);
    let mut v = proof[0].clone();
    for p in &proof[1..] {
        if index % 2 == 0 {
            v = blake(&concat_clone(&v, p));
        } else {
            v = blake(&concat_clone(p, &v));
        }
        index /= 2;
    }
    assert_eq!(v, root);
    proof[0].clone()
}

/// Make a compressed proof for multiple indices
fn mk_multi_branch(tree: &[Vec<u8>], indices: &[usize]) -> Vec<Vec<Vec<u8>>> {
    // Branches we are outputting
    let mut output = vec![];
    // Elements in the tree we can get from the branches themselves
    let mut calculable_indices = HashMap::new();
    for &i in indices {
        let new_branch = mk_branch(tree, i);
        let mut index = tree.len() / 2 + i;
        calculable_indices.insert(index, true);
        for _j in 1..new_branch.len() {
            calculable_indices.insert(index ^ 1, true);
            index /= 2;
        }
        output.push(new_branch);
    }

    // Fill in the calculable list: if we can get or calculate both children, we can calculate the parent
    let mut complete = false;
    while !complete {
        complete = true;
        let keys: Vec<_> = calculable_indices.keys().map(|&x| x / 2).sorted().rev().collect();
        for k in keys {
            if calculable_indices.contains_key(&(k * 2))
                && calculable_indices.contains_key(&(k * 2 + 1))
                && !calculable_indices.contains_key(&k)
            {
                calculable_indices.insert(k, true);
                complete = false;
            }
        }
    }

    // If for any branch node both children are calculable, or the node overlaps with a leaf, or the node
    // overlaps with a previously calculated one, elide it
    let mut scanned = HashMap::new();
    let tree_half_len = tree.len() / 2;
    for (&i, b) in indices.iter().zip(output.iter_mut()) {
        let mut index = tree_half_len + i;
        scanned.insert(index, true);
        for j in 1..b.len() {
            let index_xor_one = index ^ 1;
            if (calculable_indices.contains_key(&(index_xor_one * 2)) && calculable_indices.contains_key(&(index_xor_one * 2 + 1)))
                || (index_xor_one >= tree_half_len) && indices.contains(&(index_xor_one - tree_half_len))
                || scanned.contains_key(&(index_xor_one))
            {
                b[j] = vec![];
            }
            scanned.insert(index_xor_one, true);
            index /= 2;
        }
    }
    output
}

/// Verify a compressed proof
fn verify_multi_branch(root: &[u8], indices: &[usize], proof: &[Vec<Vec<u8>>]) -> Vec<Vec<u8>> {
    // The values in the Merkle tree we can fill in
    let mut partial_tree = HashMap::new();
    // Fill in elements from the branches
    for (&i, b) in indices.iter().zip(proof.iter()) {
        let half_tree_size = 2_usize.pow((b.len() - 1) as u32);
        let mut index = half_tree_size + i;
        partial_tree.insert(index, b[0].clone());
        for j in 1..b.len() {
            if !b[j].is_empty() {
                partial_tree.insert(index ^ 1, b[j].clone());
            }
            index /= 2;
        }
    }

    // If we can calculate or get both children, we can calculate the parent
    let mut complete = false;
    while !complete {
        complete = true;
        let keys: Vec<_> = partial_tree.keys().map(|&x| x / 2).sorted().rev().collect();
        for k in keys {
            if partial_tree.contains_key(&(k * 2))
                && partial_tree.contains_key(&(k * 2 + 1))
                && !partial_tree.contains_key(&k)
            {
                partial_tree.insert(k, blake(&concat_clone(&partial_tree[&(k * 2)], &partial_tree[&(k * 2 + 1)])));
                complete = false;
            }
        }
    }

    // If any branch node is missing, we can calculate it
    let mut proof = proof.to_vec();
    for (&i, b) in indices.iter().zip(proof.iter_mut()) {
        let half_tree_size = 2_usize.pow((b.len() - 1) as u32);
        let mut index = half_tree_size + i;
        for j in 1..b.len() {
            if b[j].is_empty() {
                b[j] = partial_tree[&(index ^ 1)].clone();
            }
            partial_tree.insert(index ^ 1, b[j].clone());
            index /= 2;
        }
    }
    indices.iter().zip(proof.iter()).map(|(&i, b)| {
        verify_branch(root, i, b)
    }).collect()
}

/// Byte length of a multi proof
fn _bin_length(proof: &[Vec<Vec<u8>>]) -> usize {
    proof.iter().map(|x| {
        let concat_len = x.iter().map(|y| y.len()).sum::<usize>();
        concat_len + x.len() / 8
    }).sum::<usize>() + proof.len() * 2
}

fn blake(x: &[u8]) -> Vec<u8> {
    Blake2s256::digest(x).to_vec()
}

fn concat_clone<T: Clone>(slice_1: &[T], slice_2: &[T]) -> Vec<T> {
    let mut v = Vec::with_capacity(slice_1.len() + slice_2.len());
    v.extend_from_slice(slice_1);
    v.extend_from_slice(slice_2);
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_merkle_tree() {
        let tree = merkelize(&(0..16).collect::<Vec<u64>>());
        for i in 0..65536 {
            let indices: Vec<usize> = (0..16).filter(|&j| (i >> j) % 2 == 1).collect();
            let branch = mk_multi_branch(&tree, &indices);
            assert_eq!(
                verify_multi_branch(&tree[1], &indices, &branch),
                indices.iter().map(|&j| tree[16 + j].clone()).collect::<Vec<_>>()
            );
            if i % 1024 == 1023 {
                println!("{} of 65536 16-element proofs checked", i + 1);
            }
        }
        println!("Multi Merkle tree test passed");
    }
}
