//! Translated from Nethermind STIR's fft.py

use super::*;

pub fn fft(vals: &[u64], modulus: u64, root_of_unity: u64, inv: bool) -> Vec<u64> {
    let rootz = expand_root_of_unity(root_of_unity, modulus);
    let mut vals = vals.to_vec();
    // Fill in vals with zeroes if needed
    if rootz.len() > vals.len() + 1 {
        vals.extend(vec![0; rootz.len() - vals.len() - 1]);
    }
    if inv {
        // Inverse FFT
        let invlen = pow_mod(vals.len() as u64, (modulus - 2) as u32, modulus);
        let rootz_inv: Vec<_> = rootz[1..].iter().rev().cloned().collect();
        let fft = _fft(&vals, modulus, &rootz_inv);
        fft
            .iter()
            .map(|&x| x * invlen % modulus)
            .collect()
    } else {
        _fft(&vals, modulus, &rootz[..rootz.len() - 1])
    }
}

fn _fft(vals: &[u64], modulus: u64, roots_of_unity: &[u64]) -> Vec<u64> {
    if vals.len() <= 4 {
        return _simple_ft(vals, modulus, roots_of_unity);
    }
    let roots_2 = roots_of_unity.iter().step_by(2).cloned().collect::<Vec<_>>();
    let l = _fft(&vals.iter().step_by(2).cloned().collect::<Vec<_>>(), modulus, &roots_2);
    let r = _fft(&vals.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), modulus, &roots_2);
    let mut o = vec![0; vals.len()];
    for (i, (&x, &y)) in l.iter().zip(r.iter()).enumerate() {
        let y_times_root = y * roots_of_unity[i] % modulus;
        o[i] = (x + y_times_root) % modulus;
        o[i + l.len()] = (x + modulus - y_times_root) % modulus;
    }
    o
}

fn _simple_ft(vals: &[u64], modulus: u64, roots_of_unity: &[u64]) -> Vec<u64> {
    let l = roots_of_unity.len();
    let mut o = Vec::with_capacity(l);
    for i in 0..l {
        let mut last = 0;
        for j in 0..l {
            last = (last + vals[j] * roots_of_unity[(i * j) % l]) % modulus;
        }
        o.push(last);
    }
    o
}

fn expand_root_of_unity(root_of_unity: u64, modulus: u64) -> Vec<u64> {
    let mut rootz = vec![1, root_of_unity];
    while *rootz.last().unwrap() != 1 {
        rootz.push((rootz.last().unwrap() * root_of_unity) % modulus);
    }
    rootz
}
