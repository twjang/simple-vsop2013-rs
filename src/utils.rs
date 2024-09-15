use std::{ops::{Add, Div, Mul, Sub}, simd::{f64x4, num::SimdFloat, StdFloat}};


pub fn u64_as_f64_slice(data: &[u64]) -> &[f64] {
    unsafe {
        std::slice::from_raw_parts(
            data.as_ptr() as *const f64,
            data.len(),
        )
    }
}

#[inline]
pub fn mod_float(x: f64, modval:f64) -> f64 {
    x - (x / modval).floor() * modval
}

pub const DPI: f64 = std::f64::consts::PI * 2.0;
const DPI_VEC: f64x4 = f64x4::from_array([DPI, DPI, DPI, DPI]);

#[inline]
fn _accumulate_block(data: &[f64], idx: usize, tj: f64) -> f64 {
    let ss_vec: f64x4 = f64x4::from_slice(&data[(idx+1)..(idx+5)]);
    let cc_vec: f64x4 = f64x4::from_slice(&data[(idx+5)..(idx+9)]);
    let aa_vec: f64x4 = f64x4::from_slice(&data[(idx+9)..(idx+13)]);
    let bb_vec: f64x4 = f64x4::from_slice(&data[(idx+13)..(idx+17)]);
    let tj_vec = f64x4::from_array([tj, tj, tj, tj]);

    let arg = aa_vec + bb_vec * tj_vec;
    let arg = arg - (arg / DPI_VEC).floor() * DPI_VEC;
    (arg.sin() * ss_vec + arg.cos() * cc_vec).reduce_sum()
}

pub fn accumulate(data: &[f64], tj: f64, tjp: f64, tol: f64) -> f64 {
    let n_blocks = data.len() / 17;
    let mut res: f64 = 0.0;
    let mut idx: usize = 0;

    let mut blk_start_idx= 0;
    let mut blk_end_idx = n_blocks;

    loop {
        let blk_mid_idx = (blk_start_idx + blk_end_idx) / 2;
        if data[blk_mid_idx * 17] > tol {
            blk_end_idx = blk_mid_idx;
        } else {
            blk_start_idx = blk_mid_idx;
        }
        if blk_end_idx - blk_start_idx < 2 {
            break
        }
    }

    idx = blk_start_idx * 17;

    loop {
        let update = _accumulate_block(data, idx, tj);
        res += update;
        idx += 17;
        if idx >= data.len() {
            break
        }
    }

    res * tjp
}

#[derive(Clone)]
struct Complex(f64, f64);

impl Add for Complex {
    type Output=Complex;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Complex(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Sub for Complex {
    type Output=Complex;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Complex(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl Mul for Complex {
    type Output=Complex;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Complex(
            self.0 * rhs.0 - self.1 * rhs.1, 
            self.0 * rhs.1 + self.1 * rhs.0, 
        )
    }
}

impl Div for Complex {
    type Output=Complex;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let d = rhs.0 * rhs.0 + rhs.1 * rhs.1;
        Complex(
            (self.0 * rhs.0 + self.1 * rhs.1) / d, 
            (self.1 * rhs.0 - self.0 * rhs.1) / d, 
        )
    }
}

impl Complex {
    #[inline]
    fn abs(&self) -> f64 {
        (self.0 * self.0 + self.1 * self.1).sqrt()
    }

    #[inline]
    fn imag(&self) -> f64 {
        self.1
    }

    #[inline]
    fn real(&self) -> f64 {
        self.0
    }

    #[inline] 
    fn conj(&self) -> Self {
        Complex(self.0, -self.1)
    }

    #[inline]
    fn exp(&self) -> Complex {
        let r = self.0.exp();
        let c = self.1.cos();
        let s = self.1.sin();
        Complex(r*c, r*s)
    }
}




pub fn elliptic_vars_to_cartesian(ep: &(f64, f64, f64, f64, f64, f64), rgm: f64) -> (f64, f64, f64, f64, f64, f64) {
    let xa = ep.0;
    let xl = ep.1;
    let xk = ep.2;
    let xh = ep.3;
    let xq = ep.4;
    let xp = ep.5;

    let xfi = (1.0 - xk*xk - xh*xh).sqrt();
    let xki = (1.0 - xq*xq - xp*xp).sqrt();
    let u   =  1.0/(1.0+xfi);

    let z = Complex(xk, xh);
    let ex = z.abs();
    let ex2 = ex*ex;
    let ex3 = ex2*ex;

    let mut z1 = z.conj();
    let mut z2: Complex;
    let mut z3: Complex;
    let mut zteta: Complex;

    let gl = mod_float(xl, DPI);
    let gm = gl - f64::atan2(xh, xk);
    let mut e = gl 
        + (ex-0.125 * ex3) * gm.sin() 
        + 0.5 * ex2 * (2.0*gm).sin() 
        + 0.375 * ex3 * (3.0*gm).sin();
    let mut rsa:f64; 
    
    loop {
        z2 = Complex(0.0, e);
        zteta = z2.exp();
        z3 = z1.clone() * zteta.clone();
        let dl = gl - e + z3.imag();
        rsa = 1.0 - z3.real();
        e = e + dl/rsa;
        if dl.abs() < 1e-15 {
            break
        }
    }
    
    z1 = z.clone() * Complex(u * z3.imag(), 0.0);
    z2 = Complex(z1.imag(), -z1.real());
    let zto = (zteta + z2 - z) / Complex(rsa, 0.0);
    let xcw = zto.real();
    let xsw = zto.imag();
    let xm = xp * xcw - xq * xsw;
    let xr = xa * rsa;

    let w0 = xr * (xcw - 2.0 * xp * xm);
    let w1 = xr * (xsw + 2.0 * xq * xm);
    let w2 = -2.0 * xr * xki * xm;

    let xms = xa * (xh+xsw) / xfi;
    let xmc = xa * (xk+xcw) / xfi;
    let xn = rgm / (xa.powf(1.5));

    let w3 = xn * ((2.0*xp*xp-1.0)*xms + 2.0*xp*xq*xmc);
    let w4 = xn * ((1.0-2.0*xq*xq)*xmc - 2.0*xp*xq*xms);
    let w5 = 2.0*xn*xki*(xp*xms + xq*xmc);

    (w0, w1, w2, w3, w4, w5)
}
