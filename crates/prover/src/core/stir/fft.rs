//! Translated from Nethermind STIR's circle_fft.py
use super::*;

/// Uses twin coset
pub fn fft_inv<F: KindaField>(coefs: &[i128], modulus: u32, root_of_unity: F, offset: F) -> Vec<i128> {
    let rootz = get_power_cycle(root_of_unity, offset);
    let xs: Vec<i128> = rootz.iter().map(|r| r.x()).collect();
    let even = x_fft_inv(&coefs.iter().step_by(2).cloned().collect::<Vec<_>>(), modulus, &xs);
    let odd = x_fft_inv(&coefs.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), modulus, &xs);
    let first_cycle: Vec<i128> = even.iter().zip(&odd).zip(&rootz).map(|((&a, &b), r)| (a + b * r.y()).rem_euclid(modulus as i128)).collect();
    let second_cycle: Vec<i128> = even.iter().zip(&odd).zip(&rootz).map(|((&a, &b), r)| (a - b * r.y()).rem_euclid(modulus as i128)).collect();
    [first_cycle, second_cycle].concat()
}

fn interleave(l1: &[i128], l2: &[i128]) -> Vec<i128> {
    l1.iter().zip(l2.iter()).flat_map(|(&a, &b)| vec![a, b]).collect()
}

fn x_fft(vals: &[i128], modulus: u32, xs: &[i128]) -> Vec<i128> {
    assert_eq!(vals.len(), xs.len());
    if vals.len() == 1 {
        return vals.to_vec();
    }
    let modulus = modulus as i128;
    let new_len = vals.len() / 2;
    let half = (modulus + 1) / 2;
    let first_half = &vals[..new_len];
    let second_half = &vals[new_len..];
    let even: Vec<_> = first_half.iter().zip(second_half).map(|(&a, &b)| {
        ((a + b) * half).rem_euclid(modulus)
    }).collect();
    let odd: Vec<_> = first_half.iter().zip(second_half).zip(xs).map(|((&a, &b), &x)| {
        ((a - b) * half * inv(x, modulus as u32)).rem_euclid(modulus)
    }).collect();
    let new_xs: Vec<_> = xs[..new_len].iter().map(|x| (2 * x * x - 1).rem_euclid(modulus)).collect();
    interleave(&x_fft(&even, modulus as u32, &new_xs),
               &x_fft(&odd, modulus as u32, &new_xs))
}

fn x_fft_inv(coefs: &[i128], modulus: u32, xs: &[i128]) -> Vec<i128> {
    if coefs.len() == 1 {
        return vec![coefs[0]; xs.len()];
    }
    let modulus = modulus as i128;
    let new_xs_len = if xs.len() == 1 { 1 } else { xs.len() / 2 };
    let new_xs: Vec<_> = xs[..new_xs_len].iter().map(|x| (2 * x * x - 1).rem_euclid(modulus)).collect();
    let even = x_fft_inv(&coefs.iter().step_by(2).cloned().collect::<Vec<_>>(), modulus as u32, &new_xs);
    let odd = x_fft_inv(&coefs.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(), modulus as u32, &new_xs);
    let even_twice = [even.clone(), even].concat();
    let odd_twice = [odd.clone(), odd].concat();
    even_twice.iter().zip(&odd_twice).zip(xs).map(|((&a, &b), &x)| {
        (a + b * x).rem_euclid(modulus)
    }).collect()
}

/// Uses twin coset
pub fn fft<F: KindaField>(vals: &[i128], modulus: u32, root_of_unity: F, offset: F) -> Vec<i128> {
    let rootz = get_power_cycle(root_of_unity, offset);
    let xs: Vec<i128> = rootz.iter().map(|r| r.x()).collect();
    let half = (modulus as i128 + 1) / 2;
    let new_len = vals.len() / 2;
    assert_eq!(new_len, rootz.len());
    let first_half = &vals[..new_len];
    let second_half = &vals[new_len..];
    let even: Vec<i128> = first_half.iter().zip(second_half).map(|(&a, &b)| {
        ((a + b) * half).rem_euclid(modulus as i128)
    }).collect();
    let odd: Vec<i128> = first_half.iter().zip(second_half).zip(&rootz).map(|((&a, &b), r)| {
        ((a - b) * half * inv(r.y(), modulus)).rem_euclid(modulus as i128)
    }).collect();
    interleave(&x_fft(&even, modulus, &xs), &x_fft(&odd, modulus, &xs))
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

pub fn shift_domain<F: KindaField>(
    vals: &[i128],
    modulus: u32,
    root_of_unity: F,
    offset: F,
    factor: F,
    expand: u32, /* default = 1*/
) -> Vec<i128> {
    fft_inv(&fft(vals, modulus, root_of_unity.pow(expand as u128), offset), modulus, root_of_unity, factor)
}

/// Evaluates f(x) for f in evaluation form
pub fn inv_fft_at_point<F: KindaField>(vals: &[i128], modulus: u32, root_of_unity: F, offset: F, x: F) -> i128 {
    fft_inv(&fft(vals, modulus, root_of_unity, offset), modulus, F::one(), x)[0]
}
