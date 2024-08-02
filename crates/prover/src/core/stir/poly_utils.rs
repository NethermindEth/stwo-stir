//! Translated from Nethermind STIR's poly_utils.py

use std::collections::HashMap;
use super::*;

/// An object that includes convenience operations for numbers
/// and polynomials in some prime field
pub struct PrimeField {
    modulus: u32,
}

#[allow(dead_code)]
impl PrimeField {
    pub fn new(modulus: u32) -> Self {
        assert_eq!(pow_mod(2_u64, modulus, modulus as u64), 2);
        Self { modulus }
    }

    pub fn modulus(&self) -> u32 {
        self.modulus
    }

    fn add(&self, x: i64, y: i64) -> i64 {
        (x + y) % (self.modulus as i64)
    }

    fn sub(&self, x: i64, y: i64) -> i64 {
        (x + y) % (self.modulus as i64)
    }

    fn mul(&self, x: i64, y: i64) -> i64 {
        (x * y) % (self.modulus as i64)
    }

    fn prod(&self, factors: &[i64]) -> i64 {
        factors.iter().fold(1, |acc, &f| self.mul(acc, f))
    }

    // can be simplified using the formula
    fn geom_sum(&self, x: i64, p: u64) -> i64 {
        let mut ans = 1;
        let mut prod = 1;
        for _ in 0..p {
            prod = self.mul(prod, x);
            ans = self.add(ans, prod);
        }
        ans
    }

    fn exp(&self, x: u64, p: u32) -> u64 {
        pow_mod(x, p, self.modulus as u64)
    }

    /// Modular inverse using the extended Euclidean algorithm
    fn inv(&self, a: i64) -> i64 {
        let (mut lm, mut hm) = (1, 0);
        let (mut low, mut high) = (a % (self.modulus as i64), self.modulus as i64);
        if low == 0 {
            panic!("ZeroDivisionError");
        }
        while low > 1 {
            let r = high / low;
            let (nm, new) = (hm - lm * r, high - low * r);
            lm = nm;
            low = new;
            hm = lm;
            high = low;
        }
        lm % (self.modulus as i64)
    }

    fn multi_inv(&self, values: &[i64]) -> Vec<i64> {
        let mut partials = vec![1];
        for &value in values {
            partials.push(self.mul(*partials.last().unwrap(), if value != 0 { value } else { 1 }));
        }
        let mut inv = self.inv(*partials.last().unwrap());
        let mut outputs = vec![0; values.len()];
        // for i in range(len(values), 0, -1):
        for i in (0..values.len()).rev() {
            outputs[i] = if values[i] != 0 { self.mul(partials[i], inv) } else { 0 };
            inv = self.mul(inv, if values[i] != 0 { values[i] } else { 1 });
        }
        outputs
    }

    fn div(&self, x: i64, y: i64) -> i64 {
        self.mul(x, self.inv(y))
    }

    /// Evaluate a polynomial at a point
    fn eval_poly_at(&self, p: &[i64], x: i64) -> i64 {
        let mut y = 0;
        let mut power_of_x = 1;
        for &p_coeff in p {
            y += power_of_x * p_coeff;
            power_of_x = (power_of_x * x) % (self.modulus as i64);
        }
        y % (self.modulus as i64)
    }

    /// Arithmetic for polynomials
    fn add_polys(&self, a: &[i64], b: &[i64]) -> Vec<i64> {
        (0..std::cmp::max(a.len(), b.len()))
            .map(|i| {
                ((if i < a.len() { a[i] } else { 0 }) + (if i < b.len() { b[i] } else { 0 }))
                    % (self.modulus as i64)
            })
            .collect()
    }

    fn sub_polys(&self, a: &[i64], b: &[i64]) -> Vec<i64> {
        (0..std::cmp::max(a.len(), b.len()))
            .map(|i| {
                ((if i < a.len() { a[i] } else { 0 }) - (if i < b.len() { b[i] } else { 0 }))
                    % (self.modulus as i64)
            })
            .collect()
    }

    fn mul_by_const(&self, a: &[i64], c: i64) -> Vec<i64> {
        a.iter().map(|&x| (x * c) % (self.modulus as i64)).collect()
    }

    fn mul_polys(&self, a: &[i64], b: &[i64]) -> Vec<i64> {
        let mut o = vec![0; a.len() + b.len() - 1];
        for (i, &aval) in a.iter().enumerate() {
            for (j, &bval) in b.iter().enumerate() {
                o[i + j] += aval * bval;
            }
        }
        o.iter().map(|&x| x % (self.modulus as i64)).collect()
    }

    fn div_polys(&self, a: &[i64], b: &[i64]) -> Vec<i64> {
        assert!(a.len() >= b.len());
        let mut a = a.to_vec();
        let mut o = vec![];
        let mut apos = a.len() - 1;
        let bpos = b.len() - 1;
        let mut diff = apos as isize - bpos as isize;
        while diff >= 0 {
            let quot = self.div(a[apos], b[bpos]);
            o.insert(0, quot);
            for i in (0..=bpos).rev() {
                a[diff as usize + i] -= b[i] * quot;
            }
            apos -= 1;
            diff -= 1;
        }
        o.iter().map(|&x| x % (self.modulus as i64)).collect()
    }

    fn mod_polys(&self, a: &[i64], b: &[i64]) -> Vec<i64> {
        self.sub_polys(a, &self.mul_polys(b, &self.div_polys(a, b)))[..b.len() - 1].to_vec()
    }

    /// Build a polynomial from a few coefficients
    fn sparse(&self, coeff_dict: &HashMap<usize, i64>) -> Vec<i64> {
        let mut o = vec![0; *coeff_dict.keys().max().unwrap() + 1];
        for (&k, &v) in coeff_dict {
            o[k] = v % (self.modulus as i64);
        }
        o
    }

    fn zpoly(&self, xs: &[i64]) -> Vec<i64> {
        let mut root = vec![1];
        for &x in xs {
            root.insert(0, 0);
            for j in 0..root.len() - 1 {
                root[j] -= root[j + 1] * x;
            }
        }
        root.iter().map(|&x| x % (self.modulus as i64)).collect()
    }

    /// Given p+1 y values and x values with no errors, recovers the original
    /// p+1 degree polynomial.
    /// Lagrange interpolation works roughly in the following way.
    /// 1. Suppose you have a set of points, eg. x = [1, 2, 3], y = [2, 5, 10]
    /// 2. For each x, generate a polynomial which equals its corresponding
    ///    y coordinate at that point and 0 at all other points provided.
    /// 3. Add these polynomials together.
    fn lagrange_interp(&self, xs: &[i64], ys: &[i64]) -> Vec<i64> {
        // Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
        let root = self.zpoly(xs);
        assert_eq!(root.len(), ys.len() + 1);
        // Generate per-value numerator polynomials, eg. for x=x2,
        // (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
        // polynomial back by each x coordinate
        let nums: Vec<Vec<_>> = xs.iter().map(|&x| self.div_polys(&root, &[-x, 1])).collect();
        // Generate denominators by evaluating numerator polys at each x
        let denoms: Vec<_> = xs.iter().enumerate().map(|(i, &x_i)| self.eval_poly_at(&nums[i], x_i)).collect();
        let invdenoms = self.multi_inv(&denoms);
        // Generate output polynomial, which is the sum of the per-value numerator
        // polynomials rescaled to have the right y values
        let mut b = vec![0; ys.len()];
        for i in 0..xs.len() {
            let yslice = self.mul(ys[i], invdenoms[i]);
            for j in 0..ys.len() {
                if nums[i][j] != 0 && ys[i] != 0 {
                    b[j] += nums[i][j] * yslice;
                }
            }
        }
        b.iter().map(|&x| x % (self.modulus as i64)).collect()
    }

    /// Optimized poly evaluation for degree 4
    fn eval_quartic(&self, p: &[i64], x: i64) -> i64 {
        let xsq = x * x % (self.modulus as i64);
        let xcb = xsq * x % (self.modulus as i64);
        (p[0] + p[1] * x + p[2] * xsq + p[3] * xcb) % (self.modulus as i64)
    }

    /// Optimized version of the above restricted to deg-4 polynomials
    fn lagrange_interp_4(&self, _xs: &[u64], _ys: &[u64]) -> Vec<u64> {
        /*
        x01, x02, x03, x12, x13, x23 = \
            xs[0] * xs[1], xs[0] * xs[2], xs[0] * xs[3], xs[1] * xs[2], xs[1] * xs[3], xs[2] * xs[3]
        m = self.modulus
        eq0 = [-x12 * xs[3] % m, (x12 + x13 + x23), -xs[1]-xs[2]-xs[3], 1]
        eq1 = [-x02 * xs[3] % m, (x02 + x03 + x23), -xs[0]-xs[2]-xs[3], 1]
        eq2 = [-x01 * xs[3] % m, (x01 + x03 + x13), -xs[0]-xs[1]-xs[3], 1]
        eq3 = [-x01 * xs[2] % m, (x01 + x02 + x12), -xs[0]-xs[1]-xs[2], 1]
        e0 = self.eval_poly_at(eq0, xs[0])
        e1 = self.eval_poly_at(eq1, xs[1])
        e2 = self.eval_poly_at(eq2, xs[2])
        e3 = self.eval_poly_at(eq3, xs[3])
        e01 = e0 * e1
        e23 = e2 * e3
        invall = self.inv(e01 * e23)
        inv_y0 = ys[0] * invall * e1 * e23 % m
        inv_y1 = ys[1] * invall * e0 * e23 % m
        inv_y2 = ys[2] * invall * e01 * e3 % m
        inv_y3 = ys[3] * invall * e01 * e2 % m
        return [(eq0[i] * inv_y0 + eq1[i] * inv_y1 + eq2[i] * inv_y2 + eq3[i] * inv_y3) % m for i in range(4)]
         */
        /*
        let x01 = xs[0] * xs[1] % self.modulus;
        let x02 = xs[0] * xs[2] % self.modulus;
        let x03 = xs[0] * xs[3] % self.modulus;
        let x12 = xs[1] * xs[2] % self.modulus;
        let x13 = xs[1] * xs[3] % self.modulus;
        let x23 = xs[2] * xs[3] % self.modulus;
        let m = self.modulus;
        let eq0 = [(-x12 * xs[3] % m + m) % m, (x12 + x13 + x23) % m, (m - xs[1] - xs[2] - xs[3]) % m, 1];
        let eq1 = [(-x02 * xs[3] % m + m) % m, (x02 + x03 + x23) % m, (m - xs[0] - xs[2] - xs[3]) % m, 1];
        let eq2 = [(-x01 * xs[3] % m + m) % m, (x01 + x03 + x13) % m, (m - xs[0] - xs[1] - xs[3]) % m, 1];
        let eq3 = [(-x01 * xs[2] % m + m) % m, (x01 + x02 + x12) % m, (m - xs[0] - xs[1] - xs[2]) % m, 1];
        let e0 = self.eval_poly_at(&eq0, xs[0]);
        let e1 = self.eval_poly_at(&eq1, xs[1]);
        let e2 = self.eval_poly_at(&eq2, xs[2]);
        let e3 = self.eval_poly_at(&eq3, xs[3]);
        let e01 = e0 * e1 % m;
        let e23 = e2 * e3 % m;
        let invall = self.inv(e01 * e23 % m);
        let inv_y0 = ys[0] * invall % m * e1 % m * e23 % m;
        let inv_y1 = ys[1] * invall % m * e0 % m * e23 % m;
        let inv_y2 = ys[2] * invall % m * e01 % m * e3 % m;
        let inv_y3 = ys[3] * invall % m * e01 % m * e2 % m;
        (0..4).map(|i| (eq0[i] * inv_y0 + eq1[i] * inv_y1 + eq2[i] * inv_y2 + eq3[i] * inv_y3) % m).collect()
        */
        unimplemented!()
    }

    // Optimized version of the above restricted to deg-2 polynomials
    fn lagrange_interp_2(&self, _xs: &[u64], _ys: &[u64]) -> Vec<u64> {
        /*
        m = self.modulus
        eq0 = [-xs[1] % m, 1]
        eq1 = [-xs[0] % m, 1]
        e0 = self.eval_poly_at(eq0, xs[0])
        e1 = self.eval_poly_at(eq1, xs[1])
        invall = self.inv(e0 * e1)
        inv_y0 = ys[0] * invall * e1
        inv_y1 = ys[1] * invall * e0
        return [(eq0[i] * inv_y0 + eq1[i] * inv_y1) % m for i in range(2)]
        */
        /*let m = self.modulus;
        let eq0 = [(-xs[1] % m + m) % m, 1];
        let eq1 = [(-xs[0] % m + m) % m, 1];
        let e0 = self.eval_poly_at(&eq0, xs[0]);
        let e1 = self.eval_poly_at(&eq1, xs[1]);
        let invall = self.inv(e0 * e1 % m);
        let inv_y0 = ys[0] * invall % m * e1 % m;
        let inv_y1 = ys[1] * invall % m * e0 % m;
        (0..2).map(|i| (eq0[i] * inv_y0 + eq1[i] * inv_y1) % m).collect()*/
        unimplemented!()
    }

    // Optimized version of the above restricted to deg-4 polynomials
    fn multi_interp_4(&self, _xsets: &[Vec<u64>], _ysets: &[Vec<u64>]) -> Vec<Vec<u64>> {
        /*
        data = []
        invtargets = []
        for xs, ys in zip(xsets, ysets):
            x01, x02, x03, x12, x13, x23 = \
                xs[0] * xs[1], xs[0] * xs[2], xs[0] * xs[3], xs[1] * xs[2], xs[1] * xs[3], xs[2] * xs[3]
            m = self.modulus
            eq0 = [-x12 * xs[3] % m, (x12 + x13 + x23), -xs[1]-xs[2]-xs[3], 1]
            eq1 = [-x02 * xs[3] % m, (x02 + x03 + x23), -xs[0]-xs[2]-xs[3], 1]
            eq2 = [-x01 * xs[3] % m, (x01 + x03 + x13), -xs[0]-xs[1]-xs[3], 1]
            eq3 = [-x01 * xs[2] % m, (x01 + x02 + x12), -xs[0]-xs[1]-xs[2], 1]
            e0 = self.eval_quartic(eq0, xs[0])
            e1 = self.eval_quartic(eq1, xs[1])
            e2 = self.eval_quartic(eq2, xs[2])
            e3 = self.eval_quartic(eq3, xs[3])
            data.append([ys, eq0, eq1, eq2, eq3])
            invtargets.extend([e0, e1, e2, e3])
        invalls = self.multi_inv(invtargets)
        o = []
        for (i, (ys, eq0, eq1, eq2, eq3)) in enumerate(data):
            invallz = invalls[i*4:i*4+4]
            inv_y0 = ys[0] * invallz[0] % m
            inv_y1 = ys[1] * invallz[1] % m
            inv_y2 = ys[2] * invallz[2] % m
            inv_y3 = ys[3] * invallz[3] % m
            o.append([(eq0[i] * inv_y0 + eq1[i] * inv_y1 + eq2[i] * inv_y2 + eq3[i] * inv_y3) % m for i in range(4)])
        # assert o == [self.lagrange_interp_4(xs, ys) for xs, ys in zip(xsets, ysets)]
        return o
         */
        /*let mut data = vec![];
        let mut invtargets = vec![];
        for (xs, ys) in xsets.iter().zip(ysets.iter()) {
            let x01 = xs[0] * xs[1] % self.modulus;
            let x02 = xs[0] * xs[2] % self.modulus;
            let x03 = xs[0] * xs[3] % self.modulus;
            let x12 = xs[1] * xs[2] % self.modulus;
            let x13 = xs[1] * xs[3] % self.modulus;
            let x23 = xs[2] * xs[3] % self.modulus;
            let m = self.modulus;
            let eq0 = [(-x12 * xs[3] % m + m) % m, (x12 + x13 + x23) % m, (m - xs[1] - xs[2] - xs[3]) % m, 1];
            let eq1 = [(-x02 * xs[3] % m + m) % m, (x02 + x03 + x23) % m, (m - xs[0] - xs[2] - xs[3]) % m, 1];
            let eq2 = [(-x01 * xs[3] % m + m) % m, (x01 + x03 + x13) % m, (m - xs[0] - xs[1] - xs[3]) % m, 1];
            let eq3 = [(-x01 * xs[2] % m + m) % m, (x01 + x02 + x12) % m, (m - xs[0] - xs[1] - xs[2]) % m, 1];
            let e0 = self.eval_quartic(&eq0, xs[0]);
            let e1 = self.eval_quartic(&eq1, xs[1]);
            let e2 = self.eval_quartic(&eq2, xs[2]);
            let e3 = self.eval_quartic(&eq3, xs[3]);
            data.push((ys.clone(), eq0, eq1, eq2, eq3));
            invtargets.extend_from_slice(&[e0, e1, e2, e3]);
        }
        let invalls = self.multi_inv(&invtargets);
        let mut o = vec![];
        for (i, (ys, eq0, eq1, eq2, eq3)) in data.into_iter().enumerate() {
            let invallz = &invalls[i * 4..i * 4 + 4];
            let inv_y0 = ys[0] * invallz[0] % self.modulus;
            let inv_y1 = ys[1] * invallz[1] % self.modulus;
            let inv_y2 = ys[2] * invallz[2] % self.modulus;
            let inv_y3 = ys[3] * invallz[3] % self.modulus;
            o.push((0..4).map(|i| (eq0[i] * inv_y0 + eq1[i] * inv_y1 + eq2[i] * inv_y2 + eq3[i] * inv_y3) % self.modulus).collect());
        }
        o*/
        unimplemented!()
    }

    // Arithmetic for circular polynomials

    fn add_circ_polys(&self, a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
        vec![self.add_polys(&a[0], &b[0]), self.add_polys(&a[1], &b[1])]
    }

    fn sub_circ_polys(&self, a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
        vec![self.sub_polys(&a[0], &b[0]), self.sub_polys(&a[1], &b[1])]
    }

    fn mul_circ_by_const(&self, a: &[Vec<i64>], c: i64) -> Vec<Vec<i64>> {
        vec![self.mul_by_const(&a[0], c), self.mul_by_const(&a[1], c)]
    }

    /// Multiply two circular polynomials
    fn mul_circ_polys(&self, a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
        /*
        a1b1=self.mul_polys(a[1],b[1])
        return [self.sub_polys(self.add_polys(self.mul_polys(a[0],b[0]), a1b1), [0,0] + a1b1),
                self.add_polys(self.mul_polys(a[0],b[1]), self.mul_polys(a[1],b[0]))]
         */
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
    fn eval_circ_poly_at(&self, p: &[Vec<i64>], pt: &Point) -> i64 {
        self.add(self.eval_poly_at(&p[0], pt.x), self.eval_poly_at(&p[1], pt.x) * pt.y)
    }

    /// Create a line polynomial between two points
    fn line(&self, pt1: &Point, pt2: &Point) -> Vec<Vec<i64>> {
        let dx = self.sub(pt1.x, pt2.x);
        if dx == 0 {
            return vec![vec![pt1.x, self.modulus as i64 - 1], vec![]];
        }
        let slope = self.div(pt1.y - pt2.y, dx);
        vec![
            vec![(pt1.y - slope * pt1.x) % self.modulus as i64, slope],
            vec![self.modulus as i64 - 1],
        ]
    }

    /// Build a circular polynomial that returns 0 at all specified points
    fn circ_zpoly(&self, pts: &[Point], nzero: Option<&Point> /* default = None */) -> Vec<Vec<i64>> {
        let mut ans = vec![vec![1], vec![]];
        for i in 0..pts.len() / 2 {
            ans = self.mul_circ_polys(&ans, &self.line(&pts[2 * i], &pts[2 * i + 1]));
        }
        if pts.len() % 2 == 1 {
            match nzero {
                Some(nzero) if nzero.x == pts[pts.len() - 1].x => {
                    ans = self.mul_circ_polys(&ans, &vec![vec![pts[pts.len() - 1].y], vec![self.modulus as i64 - 1]]);
                }
                _ => {
                    ans = self.mul_circ_polys(&ans, &vec![vec![pts[pts.len() - 1].x, self.modulus as i64 - 1], vec![]]);
                }
            }
        }
        ans
    }

    /// Circular Lagrange interpolation
    fn circ_lagrange_interp(&self, pts: &[Point], vals: &[i64], normalize: bool /* default = false */) -> Vec<Vec<i64>> {
        assert_eq!(pts.len(), vals.len());
        let mut ans = vec![vec![], vec![]];
        for i in 0..pts.len() {
            let pol = self.circ_zpoly(&[&pts[..i], &pts[i + 1..]].concat(), Some(&pts[i]));
            let scale = self.div(vals[i], self.eval_circ_poly_at(&pol, &pts[i]));
            ans = self.add_circ_polys(&ans, &self.mul_circ_by_const(&pol, scale));
        }
        if normalize && pts.len() % 2 == 0 {
            let d = pts.len() / 2;
            let zpol = self.circ_zpoly(pts, None);
            let coef_a = if ans[1].len() >= d { ans[1][d - 1] } else { 0 };
            let scale = self.div(coef_a, zpol[1][d - 1]);
            ans = self.sub_circ_polys(&ans, &self.mul_circ_by_const(&zpol, scale));
        }
        ans
    }
}
