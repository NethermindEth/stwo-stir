//! Translated from Nethermind STIR's poly_utils.py

use num_traits::{One, Zero};
use super::*;

// can be simplified using the formula
pub fn geom_sum(x: M31, p: u64) -> M31 {
    let mut ans = M31::one();
    let mut prod = M31::one();
    for _ in 0..p {
        prod = prod * x;
        ans += prod;
    }
    ans
}

/// Evaluate a polynomial at a point
pub fn eval_poly_at(p: &[M31], x: M31) -> M31 {
    let mut y = M31::zero();
    let mut power_of_x = M31::one();
    for &p_coeff in p {
        y = y + power_of_x * p_coeff;
        power_of_x *= x;
    }
    y
}

/// Arithmetic for polynomials
fn add_polys(a: &[M31], b: &[M31]) -> Vec<M31> {
    (0..std::cmp::max(a.len(), b.len()))
        .map(|i| {
            (if i < a.len() { a[i] } else { M31::zero() }) + (if i < b.len() { b[i] } else { M31::zero() })
        })
        .collect()
}

fn sub_polys(a: &[M31], b: &[M31]) -> Vec<M31> {
    (0..std::cmp::max(a.len(), b.len()))
        .map(|i| {
            (if i < a.len() { a[i] } else { M31::zero() }) - (if i < b.len() { b[i] } else { M31::zero() })
        })
        .collect()
}

fn mul_polys(a: &[M31], b: &[M31]) -> Vec<M31> {
    let mut o = vec![M31::zero(); a.len() + b.len() - 1];
    for (i, &aval) in a.iter().enumerate() {
        for (j, &bval) in b.iter().enumerate() {
            o[i + j] += aval * bval;
        }
    }
    o
}

// Arithmetic for circular polynomials

fn add_circ_polys(a: &(Vec<M31>, Vec<M31>), b: &(Vec<M31>, Vec<M31>)) -> (Vec<M31>, Vec<M31>) {
    (add_polys(&a.0, &b.0), add_polys(&a.1, &b.1))
}

fn sub_circ_polys(a: &(Vec<M31>, Vec<M31>), b: &(Vec<M31>, Vec<M31>)) -> (Vec<M31>, Vec<M31>) {
    (sub_polys(&a.0, &b.0), sub_polys(&a.1, &b.1))
}

/// Multiply two circular polynomials
fn mul_circ_polys(a: &(Vec<M31>, Vec<M31>), b: &(Vec<M31>, Vec<M31>)) -> (Vec<M31>, Vec<M31>) {
    let a1b1 = mul_polys(&a.1, &b.1);
    (
        sub_polys(
            &add_polys(&mul_polys(&a.0, &b.0), &a1b1),
            &[M31::zero(), M31::zero()].iter().chain(a1b1.iter()).cloned().collect::<Vec<_>>(),
        ),
        add_polys(&mul_polys(&a.0, &b.1), &mul_polys(&a.1, &b.0)),
    )
}

/// Evaluate a circular polynomial at a point
pub fn eval_circ_poly_at<F: Xy<M31>>(p: &(Vec<M31>, Vec<M31>), pt: &F) -> M31 {
    eval_poly_at(&p.0, pt.x()) + eval_poly_at(&p.1, pt.x()) * pt.y()
}

/// Create a line polynomial between two points
fn line<F: Xy<M31>>(pt1: &F, pt2: &F) -> (Vec<M31>, Vec<M31>) {
    let dx = pt1.x() - pt2.x();
    let max_minus_one = -M31::one();
    if M31::is_zero(&dx) {
        return (vec![pt1.x(), max_minus_one], vec![]);
    }
    let slope = (pt1.y() - pt2.y()) / dx;
    (
        vec![pt1.y() - slope * pt1.x(), slope],
        vec![max_minus_one],
    )
}

/// Build a circular polynomial that returns 0 at all specified points
pub fn circ_zpoly<F: Xy<M31>>(pts: &[F], nzero: Option<&F> /* default = None */) -> (Vec<M31>, Vec<M31>) {
    let mut ans = (vec![M31::one()], vec![]);
    for i in 0..pts.len() / 2 {
        ans = mul_circ_polys(&ans, &line(&pts[2 * i], &pts[2 * i + 1]));
    }
    let max_minus_one = -M31::one();
    if pts.len() % 2 == 1 {
        match nzero {
            Some(nzero) if nzero.x() == pts[pts.len() - 1].x() => {
                ans = mul_circ_polys(&ans, &(vec![pts[pts.len() - 1].y()], vec![max_minus_one]));
            }
            _ => {
                ans = mul_circ_polys(&ans, &(vec![pts[pts.len() - 1].x(), max_minus_one], vec![]));
            }
        }
    }
    ans
}

/// Circular Lagrange interpolation
pub fn circ_lagrange_interp<F: Xy<M31> + Copy>(pts: &[F], vals: &[M31], normalize: bool /* default = false */) -> (Vec<M31>, Vec<M31>) {
    let mul_by_const = |a: &[M31], c: M31| -> Vec<M31> {
        a.iter().map(|&x| x * c).collect()
    };
    let mul_circ_by_const = |a: &(Vec<M31>, Vec<M31>), c: M31| -> (Vec<M31>, Vec<M31>) {
        (mul_by_const(&a.0, c),
         mul_by_const(&a.1, c))
    };

    assert_eq!(pts.len(), vals.len());
    let mut ans = (vec![], vec![]);
    for i in 0..pts.len() {
        let pol = circ_zpoly(&[&pts[..i], &pts[i + 1..]].concat(), Some(&pts[i]));
        let scale = vals[i] / eval_circ_poly_at(&pol, &pts[i]);
        ans = add_circ_polys(&ans, &mul_circ_by_const(&pol, scale));
    }
    if normalize && pts.len() % 2 == 0 {
        let d = pts.len() / 2;
        let zpol = circ_zpoly(pts, None);
        let coef_a = if ans.1.len() >= d { ans.1[d - 1] } else { M31::zero() };
        let scale = coef_a / zpol.1[d - 1];
        ans = sub_circ_polys(&ans, &mul_circ_by_const(&zpol, scale));
    }
    ans
}
