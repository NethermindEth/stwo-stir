//! Translated from Nethermind STIR's fft.py
#![allow(dead_code, unused)]
use super::*;

pub fn fft(vals: &[i64], modulus: u32, root_of_unity: i64, inv: bool /* default=false */) -> Vec<i64> {
    let rootz = expand_root_of_unity(root_of_unity, modulus);
    let mut vals = vals.to_vec();
    // Fill in vals with zeroes if needed
    if rootz.len() > vals.len() + 1 {
        vals.extend(vec![0; rootz.len() - vals.len() - 1]);
    }
    if inv {
        // Inverse FFT
        let invlen = (vals.len() as i64).pow_mod(modulus - 2, modulus);
        let rootz_inv: Vec<_> = rootz[1..].iter().rev().cloned().collect();
        let fft = _fft(&vals, modulus, &rootz_inv);
        fft
            .iter()
            .map(|&x| (x * invlen).rem_euclid(modulus as i64))
            .collect()
    } else {
        _fft(&vals, modulus, &rootz[..rootz.len() - 1])
    }
}

fn _fft(vals: &[i64], modulus: u32, roots_of_unity: &[i64]) -> Vec<i64> {
    if vals.len() <= 4 {
        return _simple_ft(vals, modulus, roots_of_unity);
    }
    let roots_2 = roots_of_unity.iter().step_by(2).cloned().collect::<Vec<_>>();
    let l = _fft(&vals.iter().step_by(2).cloned().collect::<Vec<_>>(), modulus, &roots_2);
    let r = _fft(&vals.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), modulus, &roots_2);
    let mut o = vec![0; vals.len()];
    let modulus = modulus as i64;
    for (i, (&x, &y)) in l.iter().zip(r.iter()).enumerate() {
        let y_times_root = (y * roots_of_unity[i]).rem_euclid(modulus);
        o[i] = (x + y_times_root).rem_euclid(modulus);
        o[i + l.len()] = (x + modulus - y_times_root).rem_euclid(modulus);
    }
    o
}

fn _simple_ft(vals: &[i64], modulus: u32, roots_of_unity: &[i64]) -> Vec<i64> {
    let l = roots_of_unity.len();
    let mut o = Vec::with_capacity(l);
    for i in 0..l {
        let mut last = 0;
        for j in 0..l {
            last = (last + vals[j] * roots_of_unity[(i * j) % l]).rem_euclid(modulus as i64);
        }
        o.push(last);
    }
    o
}

fn expand_root_of_unity(root_of_unity: i64, modulus: u32) -> Vec<i64> {
    let mut rootz = vec![1, root_of_unity];
    while *rootz.last().unwrap() != 1 {
        rootz.push((rootz.last().unwrap() * root_of_unity).rem_euclid(modulus as i64));
    }
    rootz
}

pub fn shift_domain(vals: &[i64], modulus: u32, root_of_unity: i64, factor: i64, expand: u32) -> Vec<i64> {
    if vals.len() == 1 {
        return vec![vals[0]; expand as usize];
    }
    // 1/2 in the field
    let half = (modulus + 1) / 2;
    let rt = root_of_unity.pow_mod(expand, modulus);
    let rootz = {
        let mut v = expand_root_of_unity(rt, modulus);
        let _ = v.pop();
        v
    };
    assert_eq!(rootz.len(), vals.len());
    let rootz2 = {
        let mut v = expand_root_of_unity(root_of_unity, modulus);
        let _ = v.pop();
        v
    };
    let half_length = vals.len() / 2;
    // f(-x) in evaluation form
    let f_of_minus_x_vals: Vec<i64> = vals[half_length..].iter().chain(&vals[..half_length]).cloned().collect();
    // e(x) = (f(x) + f(-x)) / 2 in evaluation form
    let evens: Vec<i64> = vals.iter().zip(&f_of_minus_x_vals).map(|(&f, &g)| {
        mul_mod(&[f + g, half as i64], modulus)
    }).collect();
    // o(x) = (f(x) - f(-x)) / 2x in evaluation form
    let odds: Vec<i64> = vals.iter().zip(&f_of_minus_x_vals).enumerate().map(|(i, (&f, &g))| {
        let root = if i == 0 { rootz[0] } else { rootz[rootz.len() - i] };
        mul_mod(&[f - g, half as i64, root], modulus)
    }).collect();
    let shifted_evens = shift_domain(&evens[..half_length], modulus, root_of_unity.pow_mod(2, modulus), factor.pow_mod(2, modulus), expand);
    let shifted_odds = shift_domain(&odds[..half_length], modulus, root_of_unity.pow_mod(2, modulus), factor.pow_mod(2, modulus), expand);
    (0..2 * shifted_evens.len()).map(|i| {
        (
            shifted_evens[i % shifted_evens.len()]
                + mul_mod(&[factor, rootz2[i % rootz2.len()], shifted_odds[i % shifted_odds.len()]], modulus)
        ).rem_euclid(modulus as i64)
    }).collect()
}

/// Evaluates f(x) for f in evaluation form
pub fn inv_fft_at_point(vals: &[i64], modulus: u32, root_of_unity: i64, x: i64) -> i64 {
    if vals.len() == 1 {
        return vals[0];
    }
    // 1/2 in the field
    let half = (modulus as i64 + 1) / 2;
    // 1/w
    let inv_root = root_of_unity.pow_mod(vals.len() as u32 - 1, modulus);
    // f(-x) in evaluation form
    let f_of_minus_x_vals: Vec<i64> = vals[vals.len() / 2..].iter().chain(&vals[..vals.len() / 2]).cloned().collect();
    // e(x) = (f(x) + f(-x)) / 2 in evaluation form
    let evens: Vec<i64> = vals.iter().zip(&f_of_minus_x_vals).map(|(&f, &g)| ((f + g) * half).rem_euclid(modulus as i64)).collect();
    // o(x) = (f(x) - f(-x)) / 2 in evaluation form
    let odds: Vec<i64> = vals.iter().zip(&f_of_minus_x_vals).map(|(&f, &g)| ((f - g) * half).rem_euclid(modulus as i64)).collect();
    // e(x^2) + coordinate * x * o(x^2) in evaluation form
    let comb: Vec<i64> = odds.iter().zip(&evens).enumerate().map(|(i, (&o, &e))| {
        (mul_mod(&[o, x, inv_root.pow_mod(i as u32, modulus)], modulus) + e).rem_euclid(modulus as i64)
    }).collect();
    inv_fft_at_point(&comb[..comb.len() / 2], modulus, root_of_unity.pow_mod(2, modulus), x.pow_mod(2, modulus))
}

pub fn inv_shift_poly(poly: &[i64], modulus: u32, factor: i64) -> Vec<i64> {
    let mut factor_power = 1;
    let mut o = Vec::with_capacity(poly.len());
    for &p in poly {
        o.push((p * factor_power).rem_euclid(modulus as i64));
        factor_power = (factor_power * factor).rem_euclid(modulus as i64);
    }
    o
}
