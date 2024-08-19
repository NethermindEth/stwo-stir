//! Translated from Nethermind STIR's poly_utils.py

use super::*;

/// An object that includes convenience operations for numbers
/// and polynomials in some prime field
pub struct PrimeField {
    modulus: u32,
}

impl PrimeField {
    pub fn new(modulus: u32) -> Self {
        Self { modulus }
    }

    pub fn add(&self, x: i128, y: i128) -> i128 {
        (x + y).rem_euclid(self.modulus as i128)
    }

    pub fn sub(&self, x: i128, y: i128) -> i128 {
        (x - y).rem_euclid(self.modulus as i128)
    }

    // can be simplified using the formula
    pub fn geom_sum(&self, x: i128, p: u64) -> i128 {
        let mut ans = 1;
        let mut prod = 1;
        for _ in 0..p {
            prod = prod * x;
            ans = self.add(ans, prod);
        }
        ans
    }

    /// Modular inverse using the extended Euclidean algorithm
    pub fn inv(&self, a: i128) -> i128 {
        inv(a, self.modulus)
    }

    pub fn div(&self, x: i128, y: i128) -> i128 {
        (x * self.inv(y)).rem_euclid(self.modulus as i128)
    }

    /// Evaluate a polynomial at a point
    pub fn eval_poly_at(&self, p: &[i128], x: i128) -> i128 {
        let mut y = 0;
        let mut power_of_x = 1;
        for &p_coeff in p {
            y = y + power_of_x * p_coeff;
            power_of_x = (power_of_x * x).rem_euclid(self.modulus as i128);
        }
        y.rem_euclid(self.modulus as i128)
    }

    /// Arithmetic for polynomials
    fn add_polys(&self, a: &[i128], b: &[i128]) -> Vec<i128> {
        (0..std::cmp::max(a.len(), b.len()))
            .map(|i| {
                ((if i < a.len() { a[i] } else { 0 }) + (if i < b.len() { b[i] } else { 0 }))
                    .rem_euclid(self.modulus as i128)
            })
            .collect()
    }

    fn sub_polys(&self, a: &[i128], b: &[i128]) -> Vec<i128> {
        (0..std::cmp::max(a.len(), b.len()))
            .map(|i| {
                ((if i < a.len() { a[i] } else { 0 }) - (if i < b.len() { b[i] } else { 0 }))
                    .rem_euclid(self.modulus as i128)
            })
            .collect()
    }

    fn mul_polys(&self, a: &[i128], b: &[i128]) -> Vec<i128> {
        let mut o = vec![0; a.len() + b.len() - 1];
        for (i, &aval) in a.iter().enumerate() {
            for (j, &bval) in b.iter().enumerate() {
                o[i + j] += aval * bval;
            }
        }
        o.iter().map(|&x| x.rem_euclid(self.modulus as i128)).collect()
    }

    // Arithmetic for circular polynomials

    fn add_circ_polys(&self, a: &[Vec<i128>], b: &[Vec<i128>]) -> Vec<Vec<i128>> {
        vec![self.add_polys(&a[0], &b[0]), self.add_polys(&a[1], &b[1])]
    }

    fn sub_circ_polys(&self, a: &[Vec<i128>], b: &[Vec<i128>]) -> Vec<Vec<i128>> {
        vec![self.sub_polys(&a[0], &b[0]), self.sub_polys(&a[1], &b[1])]
    }

    /// Multiply two circular polynomials
    fn mul_circ_polys(&self, a: &[Vec<i128>], b: &[Vec<i128>]) -> Vec<Vec<i128>> {
        let a1b1 = self.mul_polys(&a[1], &b[1]);
        vec![
            self.sub_polys(
                &self.add_polys(&self.mul_polys(&a[0], &b[0]), &a1b1),
                &[0, 0].iter().chain(a1b1.iter()).cloned().collect::<Vec<_>>(),
            ),
            self.add_polys(&self.mul_polys(&a[0], &b[1]), &self.mul_polys(&a[1], &b[0])),
        ]
    }

    /// Evaluate a circular polynomial at a point
    pub fn eval_circ_poly_at<F: KindaField>(&self, p: &[Vec<i128>], pt: F) -> i128 {
        self.add(self.eval_poly_at(&p[0], pt.x()), self.eval_poly_at(&p[1], pt.x()) * pt.y())
    }

    /// Create a line polynomial between two points
    fn line<F: KindaField>(&self, pt1: F, pt2: F) -> Vec<Vec<i128>> {
        let dx = self.sub(pt1.x(), pt2.x());
        if dx == 0 {
            return vec![vec![pt1.x(), self.modulus as i128 - 1], vec![]];
        }
        let slope = self.div(pt1.y() - pt2.y(), dx);
        vec![
            vec![(pt1.y() - slope * pt1.x()).rem_euclid(self.modulus as i128), slope],
            vec![self.modulus as i128 - 1],
        ]
    }

    /// Build a circular polynomial that returns 0 at all specified points
    pub fn circ_zpoly<F: KindaField>(&self, pts: &[F], nzero: Option<&F> /* default = None */) -> Vec<Vec<i128>> {
        let mut ans = vec![vec![1], vec![]];
        for i in 0..pts.len() / 2 {
            ans = self.mul_circ_polys(&ans, &self.line(pts[2 * i], pts[2 * i + 1]));
        }
        if pts.len() % 2 == 1 {
            match nzero {
                Some(nzero) if nzero.x() == pts[pts.len() - 1].x() => {
                    ans = self.mul_circ_polys(&ans, &vec![vec![pts[pts.len() - 1].y()], vec![self.modulus as i128 - 1]]);
                }
                _ => {
                    ans = self.mul_circ_polys(&ans, &vec![vec![pts[pts.len() - 1].x(), self.modulus as i128 - 1], vec![]]);
                }
            }
        }
        ans
    }

    /// Circular Lagrange interpolation
    pub fn circ_lagrange_interp<F: KindaField>(&self, pts: &[F], vals: &[i128], normalize: bool /* default = false */) -> Vec<Vec<i128>> {
        let mul_by_const = |a: &[i128], c: i128| -> Vec<i128> {
            a.iter().map(|&x| (x * c).rem_euclid(self.modulus as i128)).collect()
        };
        let mul_circ_by_const = |a: &[Vec<i128>], c: i128| -> Vec<Vec<i128>> {
            vec![mul_by_const(&a[0], c),
                 mul_by_const(&a[1], c)]
        };

        assert_eq!(pts.len(), vals.len());
        let mut ans = vec![vec![], vec![]];
        for i in 0..pts.len() {
            let pol = self.circ_zpoly(&[&pts[..i], &pts[i + 1..]].concat(), Some(&pts[i]));
            let scale = self.div(vals[i], self.eval_circ_poly_at(&pol, pts[i]));
            ans = self.add_circ_polys(&ans, &mul_circ_by_const(&pol, scale));
        }
        if normalize && pts.len() % 2 == 0 {
            let d = pts.len() / 2;
            let zpol = self.circ_zpoly(pts, None);
            let coef_a = if ans[1].len() >= d { ans[1][d - 1] } else { 0 };
            let scale = self.div(coef_a, zpol[1][d - 1]);
            ans = self.sub_circ_polys(&ans, &mul_circ_by_const(&zpol, scale));
        }
        ans
    }
}
