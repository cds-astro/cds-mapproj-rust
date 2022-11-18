//! Implementation of the SIP standard.

use std::ops::RangeInclusive;
use crate::{CustomFloat, ImgXY};

/// SIP Polynomial coefficient.
/// In the polynomial, coefficient must be ordered like this:
/// * `0_0, 0_1, 0_2, 0_3, 1_0, 1_1, 1_2, 2_0, 2_1, 3_0`
/// * in which `p_q` correspond to the polynomial part `coeff_p_q * u^p * v^q` 
/// * Given an order `n`, the size of the array must be `n(n+1)/2`.
#[derive(Clone)]
pub struct SipCoeff {
  /// Computed order of the polynomial
  order: u16,
  /// Polynomials coefficient matrix
  c: Box<[f64]>,
}

impl SipCoeff {
  /// # Param
  /// * `c`: array polynomial coefficients of size `n(n+1)/2`
  pub fn new(c: Box<[f64]>) -> Self {
    let order = Self::order_from_n_coeff(c.len());
    debug_assert_eq!(order * (order + 1), (c.len() as u16) << 1);
    Self { order, c }
  }
  
  /// Returns the order of a bivariate polynomial from its number of coefficient.
  /// Sum of k for k = 1 to n equals n(n+1)/2 = l
  /// Thus n = (sqrt(8*l + 1) - 1) / 2
  /// # Param 
  /// * `n_coeff`: number of coefficient of the bivariate polynomial
  fn order_from_n_coeff(n_coeff: usize) -> u16 {
    ((((n_coeff << 3) + 1) as f64).sqrt()  as u16 - 1) >> 1
  }
  
  /// Returns the value of the polynomial, evaluated in `(u, v)`.
  pub fn p(&self, u: f64, v: f64) -> f64 {
    // Probably not the most efficient way to implement this. TODO: think twice...
    let mut k = 0;
    let mut p = 0_f64;
    let mut x = u;
    for i in 0..self.order {
      let l = self.order - i;
      let mut y = v;
      for _ in 0..l {
        p += x * y * self.c[k];
        k += 1;
        y *= y; // Not optimal for numerical stability :o/
      }
      x *= x; // Not optimal for numerical stability :o/
    }
    debug_assert_eq!(k, self.c.len());
    p
  }

  /// Returns the value of the `dp/du`, evaluated in `(u, v)`.
  pub fn dpdu(&self, u: f64, v: f64) -> f64 {
    // Probably not the most efficient way to implement this. TODO: think twice...
    let mut k = 0;
    let mut p = 0_f64;
    let mut x = u;
    for i in 1..self.order {
      let l = self.order - i;
      let mut y = v;
      for _ in 0..l {
        p += i as f64 * x * y * self.c[k];
        k += 1;
        y *= y; // Not optimal for numerical stability :o/
      }
      x *= x; // Not optimal for numerical stability :o/
    }
    debug_assert_eq!(k, self.c.len());
    p
  }

  /// Returns the value of the `dp/dv`, evaluated in `(u, v)`.
  pub fn dpdv(&self, u: f64, v: f64) -> f64 {
    // Probably not the most efficient way to implement this. TODO: think twice...
    let mut k = 0;
    let mut p = 0_f64;
    let mut x = u;
    for i in 0..self.order {
      let l = self.order - i;
      let mut y = v;
      for j in 1..l {
        p += j as f64 * x * y * self.c[k];
        k += 1;
        y *= y; // Not optimal for numerical stability :o/
      }
      x *= x; // Not optimal for numerical stability :o/
    }
    debug_assert_eq!(k, self.c.len());
    p
  }
  
}


/// SIP (un)projection coefficients for 1st and 2nd axis
#[derive(Clone)]
pub struct SipAB {
  /// Polynomials coefficient matrix on the 1st axis.
  a: SipCoeff,
  /// Polynomials coefficient matrix on the2ndt axis.
  b: SipCoeff
}

impl SipAB {
  /// # Params
  /// * `a`: 1st axis SIP coefficients
  /// * `b`: 2nd axis SIP coefficients
  pub fn new(a: SipCoeff, b: SipCoeff) -> Self {
    Self { a, b }
  }
}

/// For the SIP convention, see
/// "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
/// in the proceedings of ADASS XIV (2005).
#[derive(Clone)]
pub struct Sip {
  /// Projection coefficient.
  ab_proj: SipAB,
  /// Unprojection coefficients (if any).
  ab_deproj: Option<SipAB>, // Unprojection coefficients are optional :o/
  /// Approximatve bounds of the 1st axis domain of validity.
  /// * `start`: `-(CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  /// * `end`: `(NAXIS1 - CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  u: RangeInclusive<f64>,
  /// Approximatve bounds of the 2nd axis domain of validity.
  /// * `start`: `-(CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  /// * `end`: `(NAXIS2 - CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  v: RangeInclusive<f64>,
  fuv: RangeInclusive<f64>,
  guv: RangeInclusive<f64>,
  /// Number of iteration of the mutli-variate Newton-Raphson method (if no unproj polynomial).
  n_iter: u8, // = 20;
  /// Precision used in the mutli-variate Newton-Raphson method (if no unproj polynomial).
  eps: f64,   // = 1e-9;
}

impl Sip {
  
  /// Implements the SIP convention with the given polynomial coefficients.
  /// # Params
  /// * `ab_proj`: SIP coefficients for the projection on the 1st and 2nd axis 
  /// * `ab_deproj`: SIP coefficients for the deprojection on the 1st and 2nd axis (if any)
  /// * `u`: 1st axis domain of validity, e.g. `[-CRPIX1..NAXIS1 - CRPIX1]` 
  /// * `v`: 2nd axis domain of validity, e.g. `[-CRPIX2..NAXIS2 - CRPIX2]`
  pub fn new(
    ab_proj: SipAB, 
    ab_deproj: Option<SipAB>, 
    u: RangeInclusive<f64>, 
    v: RangeInclusive<f64>,
  ) -> Self {
    let t = ab_proj.a.p(*u.start(), *v.start()).min(ab_proj.a.p(*u.start(), *v.end()));
    let fuv_min = ab_proj.a.p(*u.start(), 0.0).min(t);
    let t = ab_proj.a.p(*u.end(), *v.start()).max(ab_proj.a.p(*u.end(), *v.end()));
    let fuv_max = ab_proj.a.p(*u.end(), 0.0).max(t);
    
    let t =  ab_proj.b.p(*u.start(), *v.start()).min( ab_proj.b.p(*u.end(), *v.start()));
    let guv_min =  ab_proj.b.p(0.0, *v.start()).min(t);
    let t =  ab_proj.b.p(*u.start(), *v.end()).max( ab_proj.b.p(*u.end(), *v.end()));
    let guv_max =  ab_proj.b.p(0.0, *v.end()).max(t);

    Self { 
      ab_proj, 
      ab_deproj, 
      u, v,
      fuv: fuv_min..=fuv_max,
      guv: guv_min..=guv_max,
      n_iter: 20,
      eps: 1.0e-9
    }
  }
  
  pub fn has_polynomial_deproj(&self) -> bool {
    self.ab_deproj.is_some()
  }
  
  pub fn f(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.p(u, v)
  }

  pub fn g(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.p(u, v)
  }

  pub fn dfdu(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.dpdu(u, v)
  }

  pub fn dfdv(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.dpdv(u, v)
  }

  pub fn dgdu(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.dpdu(u, v)
  }

  pub fn dgdv(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.dpdv(u, v)
  }


  pub fn u(&self, fuv: f64, guv: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|ab| ab.a.p(fuv, guv))
  }

  pub fn v(&self, fuv: f64, guv: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|ab| ab.b.p(fuv, guv))
  }

  pub fn inverse(&self, fuv: f64, guv: f64) -> Option<ImgXY> { // uv
    if self.has_polynomial_deproj() {
      let u = self.u(fuv, guv).unwrap();
      let v = self.v(fuv, guv).unwrap();
      Some(ImgXY::new(u, v))
    } else {
      // Make a grid and a 2-d tree to find the starting point (then multi-variate Newton) 
      None
    }
  }
  
  /// Mutli-variate Newton-Raphson:
  /// f1(x1, ..., xn) = 0es006500
  /// ...
  /// fn(x1, ..., xn) = 0
  ///  
  /// x = (x1, ..., xn)
  /// f = (f1(x), ... fn(x))
  ///  
  /// x = x - J^-1 f
  ///
  ///  With J = df1/dx1 ... df1/dxn
  /// ...     ... ...
  /// dfn/dx1 ... dfn/dxn
  ///
  /// 2d case: M = ab => M^-1 = 1/(ad-bc)  d -b 
  /// cd                     -c  a
  ///
  pub fn bivariate_newton(&self, fuv :f64, guv: f64) -> Option<ImgXY> {
    // Check input values are in the domain of validity
    if self.fuv.contains(&fuv) && self.guv.contains(&guv) {
      // Initial guess
      let mut u = fuv;
      let mut v = guv;
      // Initial values
      let mut f = self.f(u, v) - fuv;
      let mut g = self.g(u, v) - guv;
      // Bivariate Newton's method
      let eps2 = self.eps.pow2();
      let mut norm2 = f.pow2() + g.pow2();
      let mut i = 0;
      while i < self.n_iter && norm2 < eps2 {
        let a = self.dfdu(u, v);
        let b = self.dfdv(u, v);
        let c = self.dgdu(u, v);
        let d = self.dgdv(u, v);
        let det = 1.0 / (a * d - b * c);
        u -= det * (f * d - g * b);
        v -= det * (g * a - f * c);
        f = self.f(u, v) - fuv;
        g = self.g(u, v) - guv;
        norm2 = f.pow2() + g.pow2();
        i += 1;
      }
      // Check that the result is in the domain of validity
      if self.u.contains(&u) && self.v.contains(&v) {
        // All good :)
        Some(ImgXY::new(u, v))
      } else {
        // TODO: look for a different initial guess by making a grid of 300x300 an put results
        // in a kd-tree. Then look at the nearest neighbour: it is the initial guess.
        // And redo Newton.
        None
      }
    } else {
      None
    }
  }
}