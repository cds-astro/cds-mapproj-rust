//! Implementation of the SIP standard.
//!

use crate::{CustomFloat, ImgXY};
use std::ops::RangeInclusive;

/// SIP Polynomial coefficient.
/// In the polynomial, coefficient must be ordered like this:
///
/// `0_0, 1_0, 2_0, 3_0, 0_1, 1_1, 2_1, 0_2, 1_2, 0_3`
/// in which `p_q` correspond to the polynomial part `coeff_p_q * u^p * v^q`.
///
/// Note on "order"/degree: this implementation stores an internal `order` value
/// equal to (degree + 1). For an internal `order = n` the flattened coefficient
/// array length must be `n * (n + 1) / 2`. Equivalently, for polynomial degree
/// `D` the number of coefficients is `(D+1)*(D+2)/2`.
#[derive(Debug, Clone)]
pub struct SipCoeff {
  /// Computed order of the polynomial
  order: u16,

  /// Polynomials coefficient matrix
  coefficients: Box<[f64]>,
}

impl SipCoeff {
  /// Create a new SIP coefficient polynomial from the provided coefficients.
  ///
  /// # Params
  /// * `coefficients`: flattened polynomial coefficients. Coefficients are
  ///   ordered with the power of `u` (p) increasing fastest for each power
  ///   of `v` (q). For example, for degree D = 3 the first coefficients are
  ///   listed as: `0_0, 1_0, 2_0, 3_0, 0_1, 1_1, 2_1, 0_2, 1_2, 0_3` where
  ///   `p_q` corresponds to the polynomial term `coeff_p_q * u^p * v^q`.
  ///   In general for polynomial degree `D` the number of coefficients is
  ///   `(D+1)*(D+2)/2` and this function expects a boxed slice of that size.
  #[must_use]
  pub fn new(coefficients: Box<[f64]>) -> Self {
    let n_coeff = coefficients.len();

    let order = ((((n_coeff * 8) + 1) as f64).sqrt() as u16 - 1) / 2;

    debug_assert_eq!(
      order * (order + 1) / 2,
      (coefficients.len() as u16),
      "Invalid number of coefficients for a SIP polynomial"
    );
    Self {
      order,
      coefficients,
    }
  }

  /// Returns the value of the polynomial, evaluated in `(u, v)`.
  ///
  /// We define u and v as relative pixel coordinates with origin at CRPIX1, CRPIX2.
  /// x' and y' are “intermediate world coordinates” in degrees, with origin at CRVAL1, CRVAL2
  #[must_use]
  pub fn p(&self, u: f64, v: f64) -> f64 {
    let mut k = 0;
    let mut p = 0_f64;
    let mut v_pow = 1.0;
    let mut u_pow: f64;
    // loop over v^i
    for i in 0..self.order {
      u_pow = 1.0;
      // loop over u^j
      for _ in 0..self.order - i {
        p += u_pow * v_pow * self.coefficients[k];
        k += 1;
        u_pow *= u;
      }
      v_pow *= v;
    }
    debug_assert_eq!(
      k,
      self.coefficients.len(),
      "Should have iterated over all of c"
    );
    p
  }

  /// Returns the value of the `dp/du`, evaluated in `(u, v)`.
  #[must_use]
  pub fn dpdu(&self, u: f64, v: f64) -> f64 {
    let mut k = 0;
    let mut p = 0_f64;
    let mut v_pow = 1.0;
    let mut u_pow: f64;
    for i in 0..self.order {
      u_pow = 1.0;
      // The flattened coefficient array stores coefficients with the
      // power of `u` (p) increasing fastest for each power of `v` (q).
      // For dp/du we skip the p=0 term for each q (its derivative
      // wrt u is zero), so advance `k` once to move past c_{0,q}.
      k += 1;
      for j in 1..self.order - i {
        p += u_pow * v_pow * self.coefficients[k] * f64::from(j);
        k += 1;
        u_pow *= u;
      }
      v_pow *= v;
    }
    p
  }

  /// Returns the value of the `dp/dv`, evaluated in `(u, v)`.
  #[must_use]
  pub fn dpdv(&self, u: f64, v: f64) -> f64 {
    // For dp/dv we skip the entire q=0 block (which contains `order`
    // coefficients: p=0..order-1 for q=0). Start `k` at `order` to
    // position at the first coefficient with q=1.
    let mut k = (self.order) as usize;
    let mut p = 0_f64;
    let mut v_pow: f64 = 1.0;
    let mut u_pow: f64;
    for i in 1..self.order {
      u_pow = 1.0;
      for _ in 0..self.order - i {
        p += u_pow * v_pow * self.coefficients[k] * f64::from(i);
        k += 1;
        u_pow *= u;
      }
      v_pow *= v;
    }
    p
  }
}

/// SIP (un)projection coefficients for 1st and 2nd axis
#[derive(Debug, Clone)]
pub struct SipAB {
  /// Polynomials coefficient matrix on the 1st axis.
  a: SipCoeff,

  /// Polynomials coefficient matrix on the2ndt axis.
  b: SipCoeff,
}

impl SipAB {
  /// # Params
  /// * `a`: 1st axis SIP coefficients
  /// * `b`: 2nd axis SIP coefficients
  #[must_use]
  pub fn new(a: SipCoeff, b: SipCoeff) -> Self {
    Self { a, b }
  }
}

/// For the SIP convention, see
/// "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
/// in the proceedings of ADASS XIV (2005).
#[derive(Debug, Clone)]
pub struct Sip {
  /// Projection coefficient.
  ab_proj: SipAB,

  /// De-projection coefficients (if any).
  ab_deproj: Option<SipAB>,

  /// Approximate bounds of the 1st axis domain of validity.
  /// * `start`: `-(CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  /// * `end`: `(NAXIS1 - CRPIX1 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  u: RangeInclusive<f64>,

  /// Approximate bounds of the 2nd axis domain of validity.
  /// * `start`: `-(CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  /// * `end`: `(NAXIS2 - CRPIX2 + EPS)`, with EPS a number of pixels allowing to enlarge the image bounds
  v: RangeInclusive<f64>,

  /// Approximate bounds of image function.
  fuv: RangeInclusive<f64>,
  guv: RangeInclusive<f64>,

  n_iter: u8,

  /// Precision used in the mutli-variate Newton-Raphson method (if no unproj polynomial).
  eps: f64,
}

impl Sip {
  /// Implements the SIP convention with the given polynomial coefficients.
  ///
  /// # Params
  /// * `ab_proj`: SIP coefficients for the projection on the 1st and 2nd axis
  /// * `ab_deproj`: SIP coefficients for the deprojection on the 1st and 2nd axis (if any)
  /// * `u`: 1st axis domain of validity, e.g. `[-CRPIX1..NAXIS1 - CRPIX1]`
  /// * `v`: 2nd axis domain of validity, e.g. `[-CRPIX2..NAXIS2 - CRPIX2]`
  #[must_use]
  pub fn new(
    ab_proj: SipAB,
    ab_deproj: Option<SipAB>,
    u: RangeInclusive<f64>,
    v: RangeInclusive<f64>,
  ) -> Self {
    // Compute approximate range of f(u,v) and g(u,v) over the provided
    // (u,v) domain by evaluating the projection polynomials at the
    // four corners and a couple of midpoints (u at 0 and v at 0). This
    // is explicit and easier to reason about than chained min/max calls.
    let u_min = *u.start();
    let u_max = *u.end();
    let v_min = *v.start();
    let v_max = *v.end();

    // Build a consistent set of sample (u,v) points used to compute
    // exact (sampled) ranges of the distortion functions f(u,v) and g(u,v).
    let sample_points = [
      (u_min, v_min),
      (u_min, v_max),
      (u_max, v_min),
      (u_max, v_max),
      (u_min, 0.0),
      (u_max, 0.0),
      (0.0, v_min),
      (0.0, v_max),
    ];

    let mut fuv_min = f64::INFINITY;
    let mut fuv_max = f64::NEG_INFINITY;
    let mut guv_min = f64::INFINITY;
    let mut guv_max = f64::NEG_INFINITY;

    for (su, sv) in sample_points.iter().copied() {
      let a_val = ab_proj.a.p(su, sv);
      let b_val = ab_proj.b.p(su, sv);
      if a_val < fuv_min {
        fuv_min = a_val;
      }
      if a_val > fuv_max {
        fuv_max = a_val;
      }
      if b_val < guv_min {
        guv_min = b_val;
      }
      if b_val > guv_max {
        guv_max = b_val;
      }
    }

    Self {
      ab_proj,
      ab_deproj,
      u,
      v,
      fuv: fuv_min..=fuv_max,
      guv: guv_min..=guv_max,
      n_iter: 20,
      eps: 1.0e-9,
    }
  }

  /// Does SIP contain polynomial deprojection
  #[must_use]
  pub fn has_polynomial_deproj(&self) -> bool {
    self.ab_deproj.is_some()
  }

  #[must_use]
  pub(crate) fn f(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.p(u, v)
  }

  #[must_use]
  pub(crate) fn g(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.p(u, v)
  }

  #[must_use]
  pub(crate) fn dfdu(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.dpdu(u, v)
  }

  #[must_use]
  pub(crate) fn dfdv(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.a.dpdv(u, v)
  }

  #[must_use]
  pub(crate) fn dgdu(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.dpdu(u, v)
  }

  #[must_use]
  pub(crate) fn dgdv(&self, u: f64, v: f64) -> f64 {
    self.ab_proj.b.dpdv(u, v)
  }

  #[must_use]
  pub(crate) fn u(&self, fuv: f64, guv: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|ab| ab.a.p(fuv, guv))
  }

  #[must_use]
  pub(crate) fn v(&self, fuv: f64, guv: f64) -> Option<f64> {
    self.ab_deproj.as_ref().map(|ab| ab.b.p(fuv, guv))
  }

  #[must_use]
  pub(crate) fn inverse(&self, fuv: f64, guv: f64) -> Option<ImgXY> {
    // Quick containment check: compute conservative observed-coordinate
    // bounds from the stored sampled distortion ranges `fuv`/`guv` and
    // the u/v domain. Observed coordinates are approximately
    // u + f(u,v) and v + g(u,v), so we form bounds as the sum of the
    // domain and sampled distortion ranges.
    let obs_u_min = *self.u.start() + *self.fuv.start();
    let obs_u_max = *self.u.end() + *self.fuv.end();
    let obs_v_min = *self.v.start() + *self.guv.start();
    let obs_v_max = *self.v.end() + *self.guv.end();
    if fuv < obs_u_min || fuv > obs_u_max || guv < obs_v_min || guv > obs_v_max {
      return None;
    }

    // If polynomial deprojection is available, use it; otherwise fall
    // back to the iterative Newton solver.
    if self.has_polynomial_deproj() {
      let u = self.u(fuv, guv).unwrap();
      let v = self.v(fuv, guv).unwrap();
      Some(ImgXY::new(u, v))
    } else {
      self.bivariate_newton(fuv, guv)
    }
  }

  /// Multi-variate Newton-Raphson:
  /// f1(x1, ..., xn) = 0
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
  #[must_use]
  pub fn bivariate_newton(&self, fuv: f64, guv: f64) -> Option<ImgXY> {
    // Try several small perturbations to the initial guess when solving
    // the additive distortion equations u + f(u,v) = fuv, v + g(u,v) =
    // guv. This helps avoid singular Jacobians or poor starting points.
    let attempts = [
      (0.0_f64, 0.0_f64),
      (1e-3, 0.0),
      (-1e-3, 0.0),
      (0.0, 1e-3),
      (0.0, -1e-3),
      (1e-2, 1e-2),
      (-1e-2, -1e-2),
    ];

    let eps2 = self.eps.pow2();

    for (du, dv) in attempts.iter().copied() {
      let mut u = fuv + du;
      let mut v = guv + dv;

      // Initial residuals for the additive system
      let mut f = (u + self.f(u, v)) - fuv;
      let mut g = (v + self.g(u, v)) - guv;

      let mut norm2 = f.pow2() + g.pow2();
      let mut i = 0;
      while i < self.n_iter && norm2 > eps2 {
        // Jacobian of the system [u + f(u,v), v + g(u,v)]
        let a = 1.0 + self.dfdu(u, v);
        let b = self.dfdv(u, v);
        let c = self.dgdu(u, v);
        let d = 1.0 + self.dgdv(u, v);
        let jac = a * d - b * c;
        // If Jacobian is unusable, break and try the next initial guess.
        // Use a practical tolerance rather than testing against
        // `MIN_POSITIVE` (which is too small for numerical contexts).
        // Scale epsilon by the diagonal terms to get a relative
        // threshold that adapts to the magnitude of the problem.
        let jac_tol = f64::EPSILON * (a.abs() + d.abs()).max(1.0);
        if jac.is_nan() || jac.is_infinite() || jac.abs() <= jac_tol {
          break;
        }
        let inv_jac = 1.0 / jac;
        // Take the Newton step. If the step is excessively large, damp it.
        let step_u = inv_jac * (f * d - g * b);
        let step_v = inv_jac * (g * a - f * c);
        // Damping: prevent explosion by limiting maximum step magnitude
        let max_step = 1e6_f64; // very large, defensive
        let step_norm = step_u.abs().max(step_v.abs());
        let (du_step, dv_step) = if step_norm > max_step {
          (
            step_u * (max_step / step_norm),
            step_v * (max_step / step_norm),
          )
        } else {
          (step_u, step_v)
        };
        u -= du_step;
        v -= dv_step;
        f = (u + self.f(u, v)) - fuv;
        g = (v + self.g(u, v)) - guv;
        norm2 = f.pow2() + g.pow2();
        i += 1;
      }

      // Check for convergence and domain validity
      if norm2 <= eps2 && self.u.contains(&u) && self.v.contains(&v) {
        return Some(ImgXY::new(u, v));
      }
    }
    // If we exhausted attempts without returning above, no convergence
    // was reached that satisfies the domain and tolerance.
    None
  }
}

#[cfg(test)]
mod tests {
  use crate::sip::{Sip, SipAB, SipCoeff};

  #[test]
  fn test_sip() {
    // taken from table 1 of
    // "The SIP convention for Representing Distortion in FITS Image Headers"
    // by David L. Shupe et al. in the proceedings of ADASS XIV (2005).
    //
    // CTYPE1 'RA---TAN-SIP'
    // CTYPE2 'DEC--TAN-SIP'
    // CRPIX1 2048.0
    // CRPIX2 1024.0
    // CRVAL1 5.6260667398471
    // CRVAL2 -72.076963036772
    // CD1_1 -7.8481866550866E-06
    // CD2_1 1.1406694624771E-05
    // CD1_2 1.0939720432379E-05
    // CD2_2 8.6942510845452E-06

    // A_0_0 0.0
    // A_0_1 0.0
    // A_0_2 2.1634068532689E-06
    // A_0_3 1.0622437604068E-11
    // A_0_4 1.4075878614807E-14
    // A_1_0 0.0
    // A_1_1 -5.194753640575E-06
    // A_1_2 -5.2797808038221E-10
    // A_1_3 -1.9317154005522E-14
    // A_2_0 8.543473309812E-06
    // A_2_1 -4.4012735467525E-11
    // A_2_2 3.767898933666E-14
    // A_3_0 -4.7518233007536E-10
    // A_3_1 5.0860953083043E-15
    // A_4_0 2.5776347115304E-14

    // B_0_0 0.0
    // B_0_1 0.0
    // B_0_2 -7.2299995118730E-06
    // B_0_3 -4.2102920235938E-10
    // B_0_4 6.5531313110898E-16
    // B_1_0 0.0
    // B_1_1 6.1778338717084E-06
    // B_1_2 -6.7603466821178E-11
    // B_1_3 1.3892905568706E-14
    // B_2_0 -1.7442694174934E-06
    // B_2_1 -5.1333879897858E-10
    // B_2_2 -2.9648166208490E-14
    // B_3_0 8.5722142612681E-11
    // B_3_1 -2.0749495718513E-15
    // B_4_0 -1.812610418272E-14

    // A_ORDER 4
    // B_ORDER 4

    // In the polynomial, coefficients must be flattened with the power of
    // `u` (p) increasing fastest for each power of `v` (q). For example
    // for degree 3 the ordering used here is:
    // `0_0, 1_0, 2_0, 3_0, 0_1, 1_1, 2_1, 0_2, 1_2, 0_3`.
    let a_coef = [
      0.0,
      0.0,
      2.1634068532689E-06,
      1.0622437604068E-11,
      1.4075878614807E-14,
      0.0,
      -5.194753640575E-06,
      -5.2797808038221E-10,
      -1.9317154005522E-14,
      8.543473309812E-06,
      -4.4012735467525E-11,
      3.767898933666E-14,
      -4.7518233007536E-10,
      5.0860953083043E-15,
      2.5776347115304E-14,
    ];

    let a_coef = SipCoeff::new(a_coef.into());

    let b_coef = [
      0.0,
      0.0,
      -7.2299995118730E-06,
      -4.2102920235938E-10,
      6.5531313110898E-16,
      0.0,
      6.1778338717084E-06,
      -6.7603466821178E-11,
      1.3892905568706E-14,
      -1.7442694174934E-06,
      -5.1333879897858E-10,
      -2.9648166208490E-14,
      8.5722142612681E-11,
      -2.0749495718513E-15,
      -1.812610418272E-14,
    ];

    let b_coef = SipCoeff::new(b_coef.into());
    let ab_coef = SipAB::new(a_coef, b_coef);

    let sip = Sip::new(ab_coef, None, -2048.0..=2048.0, -1024.0..=1024.0);

    let dfdu = sip.dfdu(3.0, 3.0);
    let dfdv = sip.dfdv(3.0, 3.0);

    // Approximate dfdu and dfdv to ensure we calculated it correctly.
    // As eps goes to 0, our approx derivative gets better, but numerical
    // noise creeps in, the 1e-8 is working to get us near a sweet spot of
    // mathematical precision and numerical precision for this comparison.
    // Chosen by manual tweaking.
    let eps = 1e-8;
    let da = sip.f(3.0, 3.0);
    let db = sip.f(3.0 + eps, 3.0);
    let dc = sip.f(3.0, 3.0 + eps);
    let dfdu_approx = (db - da) / eps;
    let dfdv_approx = (dc - da) / eps;

    assert!((dfdu - dfdu_approx).abs() < 1e-11);
    assert!((dfdv - dfdv_approx).abs() < 1e-11);

    // dfdu and dfdv pass our sanity checks, but we do not know if f and g are
    // correct

    let tmp = sip.f(3.0, 3.0);
    let tmp2 = sip.g(3.0, 3.0);

    // These values were computed with astropy

    // from astropy.wcs import WCS
    // import numpy as np
    // cr_pix = np.array([2048.0, 1024.0])
    // wcs = WCS(header)
    // wcs.sip_pix2foc(np.array([cr_pix + 3]), 1) - 3

    assert!((tmp - 4.95811569e-05).abs() < 1e-11);
    assert!((tmp2 - -2.51926572e-05).abs() < 1e-11);

    // test inverse
    let uv = sip.inverse(tmp + 3.0, tmp2 + 3.0).unwrap();
    assert!((uv.x() - 3.0).abs() < 1e-7);
    assert!((uv.y() - 3.0).abs() < 1e-7);
  }

  #[test]
  fn test_sipcoeff_ordering() {
    // Degree D = 2 -> coefficients count = (D+1)*(D+2)/2 = 6
    // Coefficient order (flattened): 0_0, 1_0, 2_0, 0_1, 1_1, 0_2
    // Define polynomial p(u,v) = c00 + c10*u + c20*u^2 + c01*v + c11*u*v + c02*v^2
    let coeffs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let sc = SipCoeff::new(coeffs.into());
    let u = 2.0_f64;
    let v = 3.0_f64;
    let expected = 1.0 + 2.0 * u + 3.0 * u * u + 4.0 * v + 5.0 * u * v + 6.0 * v * v;
    assert!((sc.p(u, v) - expected).abs() < 1e-12);

    // dp/du = 2 + 6*u + 5*v
    let expected_du = 2.0 + 6.0 * u + 5.0 * v;
    // dp/dv = 4 + 5*u + 12*v
    let expected_dv = 4.0 + 5.0 * u + 12.0 * v;
    assert!((sc.dpdu(u, v) - expected_du).abs() < 1e-12);
    assert!((sc.dpdv(u, v) - expected_dv).abs() < 1e-12);
  }

  #[test]
  fn test_inverse_with_deproj() {
    // Projection polynomials (degree 1): f(u,v) = u, g(u,v) = v
    let a_proj = [0.0, 1.0, 0.0]; // c00, c10, c01
    let b_proj = [0.0, 0.0, 1.0];
    let a_proj = SipCoeff::new(a_proj.into());
    let b_proj = SipCoeff::new(b_proj.into());
    let ab_proj = SipAB::new(a_proj, b_proj);

    // Deprojection polynomials: u = 2*f + 3*g + 5, v = -f + 4*g - 1
    // degree 1 flattened: c00, c10, c01
    let a_deproj = [5.0, 2.0, 3.0];
    let b_deproj = [-1.0, -1.0, 4.0];
    let a_deproj = SipCoeff::new(a_deproj.into());
    let b_deproj = SipCoeff::new(b_deproj.into());
    let ab_deproj = SipAB::new(a_deproj, b_deproj);

    let sip = Sip::new(ab_proj, Some(ab_deproj), -10.0..=10.0, -10.0..=10.0);

    let fuv = 1.5;
    let guv = -2.0;
    let uv = sip.inverse(fuv, guv).unwrap();
    // Expect u = 2*fuv + 3*guv + 5
    assert!((uv.x() - (2.0 * fuv + 3.0 * guv + 5.0)).abs() < 1e-12);
    // Expect v = -fuv + 4*guv - 1
    assert!((uv.y() - (-fuv + 4.0 * guv - 1.0)).abs() < 1e-12);
  }

  #[test]
  fn test_inverse_linear_newton() {
    // Small linear distortion: f = a*u, g = a*v with a = 1e-3
    let a = 1e-3;
    let a_proj = [0.0, a, 0.0];
    let b_proj = [0.0, 0.0, a];
    let sip = Sip::new(
      SipAB::new(SipCoeff::new(a_proj.into()), SipCoeff::new(b_proj.into())),
      None,
      -1000.0..=1000.0,
      -1000.0..=1000.0,
    );

    let u_true = 10.0;
    let v_true = -7.5;
    // Observed coordinates: u' = u + a*u = u*(1+a)
    let fuv = u_true * (1.0 + a);
    let guv = v_true * (1.0 + a);

    let uv = sip.inverse(fuv, guv).unwrap();
    assert!((uv.x() - u_true).abs() < 1e-9);
    assert!((uv.y() - v_true).abs() < 1e-9);
  }

  #[test]
  fn test_singular_jacobian_returns_none() {
    // Create projection f = -u, g = -v so that Jacobian at any (u,v) makes
    // 1 + df/du = 0 and 1 + dg/dv = 0 => determinant = 0
    let a_proj = [0.0, -1.0, 0.0];
    let b_proj = [0.0, 0.0, -1.0];
    let sip = Sip::new(
      SipAB::new(SipCoeff::new(a_proj.into()), SipCoeff::new(b_proj.into())),
      None,
      -10.0..=10.0,
      -10.0..=10.0,
    );

    // For arbitrary observed coords, inverse should be None (no solution)
    assert!(sip.inverse(1.0, 1.0).is_none());
  }

  #[test]
  fn test_inverse_out_of_bounds_observed_returns_none() {
    // Small linear distortion: f = a*u, g = a*v with a = 1e-3
    let a = 1e-3;
    let sip = Sip::new(
      SipAB::new(
        SipCoeff::new([0.0, a, 0.0].into()),
        SipCoeff::new([0.0, 0.0, a].into()),
      ),
      None,
      -1000.0..=1000.0,
      -1000.0..=1000.0,
    );

    // Use an observed coordinate far outside the sampled observed bounds
    assert!(sip.inverse(1.0e6, 1.0e6).is_none());
  }

  #[test]
  fn test_sampled_fuv_guv_ranges_for_linear() {
    // For a_proj(u,v) = u and b_proj(u,v) = v, sampled fuv/guv should
    // match the u/v sample bounds used in construction.
    let a_proj = [0.0, 1.0, 0.0]; // p(u,v) = u
    let b_proj = [0.0, 0.0, 1.0]; // q(u,v) = v
    let sip = Sip::new(
      SipAB::new(SipCoeff::new(a_proj.into()), SipCoeff::new(b_proj.into())),
      None,
      -10.0..=10.0,
      -5.0..=5.0,
    );

    let fuv_min = *sip.fuv.start();
    let fuv_max = *sip.fuv.end();
    let guv_min = *sip.guv.start();
    let guv_max = *sip.guv.end();

    assert!((fuv_min - -10.0).abs() < 1e-12);
    assert!((fuv_max - 10.0).abs() < 1e-12);
    assert!((guv_min - -5.0).abs() < 1e-12);
    assert!((guv_max - 5.0).abs() < 1e-12);
  }
}
