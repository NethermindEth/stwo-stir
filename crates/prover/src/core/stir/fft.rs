//! Translated from Nethermind STIR's circle_fft.py

use num_traits::One;
use crate::core::fields::FieldExpOps;
use crate::core::fields::m31::{M31};
use super::*;

/// Uses twin coset
pub fn fft_inv<F: StirField<M31>>(coefs: &[M31], root_of_unity: F, offset: F) -> Vec<M31> {
    let rootz = get_power_cycle(root_of_unity, offset);
    let xs: Vec<M31> = rootz.iter().map(|r| r.x()).collect();
    let even = x_fft_inv(&coefs.iter().step_by(2).cloned().collect::<Vec<_>>(), &xs);
    let odd = x_fft_inv(&coefs.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), &xs);
    let first_cycle: Vec<M31> = even.iter().zip(&odd).zip(&rootz).map(|((&a, &b), r)| a + b * r.y()).collect();
    let second_cycle: Vec<M31> = even.iter().zip(&odd).zip(&rootz).map(|((&a, &b), r)| a - b * r.y()).collect();
    [first_cycle, second_cycle].concat()
}

fn interleave(l1: &[M31], l2: &[M31]) -> Vec<M31> {
    l1.iter().zip(l2.iter()).flat_map(|(&a, &b)| vec![a, b]).collect()
}

fn x_fft(vals: &[M31], xs: &[M31]) -> Vec<M31> {
    assert_eq!(vals.len(), xs.len());
    if vals.len() == 1 {
        return vals.to_vec();
    }
    let new_len = vals.len() / 2;
    let half = M31::from_u32_unchecked((MODULUS + 1) / 2);
    let first_half = &vals[..new_len];
    let second_half = &vals[new_len..];
    let even: Vec<_> = first_half.iter().zip(second_half).map(|(&a, &b)| {
        (a + b) * half
    }).collect();
    let odd: Vec<_> = first_half.iter().zip(second_half).zip(xs).map(|((&a, &b), &x)| {
        (a - b) * half * x.inverse()
    }).collect();
    let new_xs: Vec<_> = xs[..new_len].iter().map(|x| x.square().double() - M31::one()).collect();
    interleave(&x_fft(&even, &new_xs),
               &x_fft(&odd, &new_xs))
}

fn x_fft_inv(coefs: &[M31], xs: &[M31]) -> Vec<M31> {
    if coefs.len() == 1 {
        return vec![coefs[0]; xs.len()];
    }
    let new_xs_len = if xs.len() == 1 { 1 } else { xs.len() / 2 };
    let new_xs: Vec<_> = xs[..new_xs_len].iter().map(|x| x.square().double() - M31::one()).collect();
    let even = x_fft_inv(&coefs.iter().step_by(2).cloned().collect::<Vec<_>>(), &new_xs);
    let odd = x_fft_inv(&coefs.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), &new_xs);
    let even_twice = [even.clone(), even].concat();
    let odd_twice = [odd.clone(), odd].concat();
    even_twice.iter().zip(&odd_twice).zip(xs).map(|((&a, &b), &x)| {
        a + b * x
    }).collect()
}

/// Uses twin coset
pub fn fft<F: StirField<M31>>(vals: &[M31], root_of_unity: F, offset: F) -> Vec<M31> {
    let rootz = get_power_cycle(root_of_unity, offset);
    let xs: Vec<M31> = rootz.iter().map(|r| r.x()).collect();
    let half = M31::from_u32_unchecked((MODULUS + 1) / 2);
    let new_len = vals.len() / 2;
    assert_eq!(new_len, rootz.len());
    let first_half = &vals[..new_len];
    let second_half = &vals[new_len..];
    let even: Vec<M31> = first_half.iter().zip(second_half).map(|(&a, &b)| {
        (a + b) * half
    }).collect();
    let odd: Vec<M31> = first_half.iter().zip(second_half).zip(&rootz).map(|((&a, &b), r)| {
        (a - b) * half * r.y().inverse()
    }).collect();
    interleave(&x_fft(&even, &xs), &x_fft(&odd, &xs))
}

fn _fft(vals: &[i128], modulus: u32, roots_of_unity: &[i128]) -> Vec<i128> {
    if vals.len() <= 4 {
        return _simple_ft(vals, modulus, roots_of_unity);
    }
    let roots_2 = roots_of_unity.iter().step_by(2).cloned().collect::<Vec<_>>();
    let l = _fft(&vals.iter().step_by(2).cloned().collect::<Vec<_>>(), modulus, &roots_2);
    let r = _fft(&vals.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), modulus, &roots_2);
    let mut o = vec![0; vals.len()];
    let modulus = modulus as i128;
    for (i, (&x, &y)) in l.iter().zip(r.iter()).enumerate() {
        let y_times_root = (y * roots_of_unity[i]).rem_euclid(modulus);
        o[i] = (x + y_times_root).rem_euclid(modulus);
        o[i + l.len()] = (x + modulus - y_times_root).rem_euclid(modulus);
    }
    o
}

fn _simple_ft(vals: &[i128], modulus: u32, roots_of_unity: &[i128]) -> Vec<i128> {
    let l = roots_of_unity.len();
    let mut o = Vec::with_capacity(l);
    for i in 0..l {
        let mut last = 0;
        for j in 0..l {
            last = (last + vals[j] * roots_of_unity[(i * j) % l]).rem_euclid(modulus as i128);
        }
        o.push(last);
    }
    o
}

pub fn shift_domain<F: StirField<M31>>(
    vals: &[M31],
    root_of_unity: F,
    offset: F,
    factor: F,
    expand: u32, /* default = 1*/
) -> Vec<M31> {
    fft_inv(&fft(vals, root_of_unity.pow(expand as u128), offset), root_of_unity, factor)
}

/// Evaluates f(x) for f in evaluation form
pub fn inv_fft_at_point<F: StirField<M31>>(vals: &[M31], root_of_unity: F, offset: F, x: F) -> M31 {
    fft_inv(&fft(vals, root_of_unity, offset), F::one(), x)[0]
}
