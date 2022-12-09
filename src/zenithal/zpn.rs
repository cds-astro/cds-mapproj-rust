//! Zenithal polynomlial projection.

use std::{
  ops::RangeInclusive,
  f64::consts::PI,
};

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Zenithal polynomlial projection.
#[derive(Debug, Clone)]
pub struct Zpn {
  /// Coefficients of the polynomial: a0 + a1 x^2 + a2 x^3 + ...
  p: Box<[f64]>,
  // Step to look at the domain of validity of the projection.
  // domain_step: f64, // = toRadians(1.0 / 60);
  // Epslion to look at the domain of validity of the projection
  // domain_eps: f64, //Epsilon = epsilon;
  /// Computed domain validity, min angular dist
  ang_dist: RangeInclusive<f64>, // for projection
  /// Computed domain validity, max angular dist
  euc_dist: RangeInclusive<f64>, // for deprojection
  /// Max number of iteration for the Newton-Raphson iterative method.
  n_iter: u8,
  /// Precision for the Newton-Raphson iterative method
  eps: f64,
  /// Step for the dichotomy (when newton fails)
  step: f64,
  // Proj bounds
  proj_bounds: ProjBounds,
}

impl Zpn {
  
  /// # Params
  /// * `coeffs`: polynomial coefficients provided by keywords PVi_1a, PVi_2a, ..., PVi_na
  /// # Return
  /// * `None` if negative polynomial en `[0, pi]`
  pub fn from_params(coeffs: Vec<f64>) -> Option<Self> {
    let eps = (1.0_f64 / (60.0_f64 * 60.0 * 1000.0)).to_radians(); // <=> 1 mas ~= 5e-9 radians
    let step = (1.0_f64 / 60.0).to_radians(); // <=> 1 arcmin
    Self::from_params_custom(coeffs, step, eps)
  }

  /// # Params
  /// * `coeffs`: polynomial coefficients provided by keywords PVi_1a, PVi_2a, ..., PVi_na
  /// * `domain_step`: step used to determine the domain of validity of the function
  ///                  default value: 1 arcmin converted in radians
  /// * `domain_eps`: epsilon used in dichotomy to end the convergence process
  ///                 default value: 1 mas converted into radians
  ///                 (this value is also assigned to the `epsilon` used in deprojection.
  /// # Return
  /// * `None` if negative polynomial en `[0, pi]`
  pub fn from_params_custom(coeffs: Vec<f64>, domain_step: f64, domain_eps: f64) -> Option<Self> {
    let p = coeffs.into_boxed_slice();
    // Set constants
    let n_iter = 20;
    // Find the smallest angular distance leading to a positive Euclidean distance
    let mut ang_dist_min = 0.0_f64;
    let mut euc_dist_min = polynomial(ang_dist_min, &p);
    if euc_dist_min < 0.0 { // In case a[0] < 0
      // By dichotomy
      // - starts by looking at the first positive value
      let mut min = 0.0_f64;
      let mut max = domain_step;
      let mut max_val = polynomial(max, &p);
      while max_val < 0.0 && max < PI {
        max += domain_step;
        max_val = polynomial(max, &p);
      }
      if max_val < 0.0 {
        return None;
      }
      // - then compute the smallest value leading to a positive Euclidean distance
      while max - min > domain_eps {
        let med = (max + min).half();
        let med_val = polynomial(med, &p);
        if med_val < 0.0 {
          min = med;
        } else {
          max = med;
          max_val = med_val;
        }
      }
      ang_dist_min = max;
      euc_dist_min = max_val;
    }
    // Find the smallest angular distance (ang_dist_max) reaching a local maximum such that the
    // polynomial is bijective on [ang_dist_min, ang_dist_max].
    // To do so, I use the smallest positive value of 'x' at which the sign of the derivative
    // changes.
    // 
    // The more robust solution would have been to look at the smallest positive real root
    // of the polynomial derivative but it seems too much code to implements.
    // Doing so, there will be no need for the two parameters 'domain_step' and 'domain_eps'.
    //
    // Two solutions would have been retained:
    // - eigen-decomposition of a matrix to find all roots (see e.g. Numerical recipes,
    //   Eigenvalue method in Roots of Polynomials).
    // - Usage of the Sturm's theorem to check there is no smaller positive real root:
    //   https://math.stackexchange.com/questions/1214897/how-to-upper-bound-the-smallest-positive-root-of-a-polynomial
    //   https://en.wikipedia.org/wiki/Sturm%27s_theorem
    //   http://www.sciencedirect.com/science/article/pii/S0377042703007271
    // 
    // Here, I am using the fact that the polynomial represents a distance from a center, so:
    // - I assume it is increasing on the bijective part from the projection center to a maximum
    // - I assume the smallest value of x is 0 and the largest is \pi
    // - I assume the highest 'frequency' of the derivative is < self.domain_step
    // - I assume the derivative is positive (increasing distance) on the bijective part
    // We do not use Newton's method since it does not guarantee to provide the smallest positive root
    let mut ang_dist_max = ang_dist_min;
    let euc_dist_deriv = dpolynomial(ang_dist_max, &p);
    if euc_dist_deriv < 0.0 {
      // "Negative derivative of the polynomial at distance"
      return None;
    } else {
      // By dichotomy
      // - find the first value for which the derivative become negative
      let mut min = ang_dist_min;
      let mut max = ang_dist_min + domain_step;
      let mut max_val = dpolynomial(max, &p);
      while max_val > 0.0 && max < PI {
        max += domain_step;
        max_val = dpolynomial(max, &p);
      }
      if max_val >= 0.0 { // No change of sign found
        min = PI;
      } else { // Dichotomy to find the point in which the derivative becomes null
        while (max - min) > domain_eps {
          let med = (max + min).half();
          let med_val = dpolynomial(med, &p);
          if med_val <= 0.0 {
            max = med;
          } else {
            min = med;
          }
        }
      }
      ang_dist_max = min;
    }
    let euc_dist_max = polynomial(ang_dist_max, &p);
    // Build struct
    Some(Self {
      p,
      ang_dist: ang_dist_min..=ang_dist_max,
      euc_dist: euc_dist_min..=euc_dist_max,
      n_iter,
      eps: domain_eps,
      step: domain_step,
      proj_bounds: ProjBounds::new(
        Some(-euc_dist_max..=euc_dist_max),
        Some(-euc_dist_max..=euc_dist_max)
      )
    })
  }
  
  
  /// Set the max number of iteration for the Newton-Raphson iterative method.
  pub fn set_n_iter(&mut self, n_iter: u8) {
    self.n_iter = n_iter;
  }

  /// Set precision for the Newton-Raphson iterative method
  pub fn set_eps(&mut self, eps: f64) {
    self.eps = eps;
  }

  /// Solve equation:
  ///  `polynomial(x) - p = 0`
  /// using Newton's method.
  /// We guess the polynomial function is smooth with no local extremum, ...
  fn solve(&self, p: f64) -> f64 {
    // First guess: Euclidean distance probably similar to angular distance
    let mut r = p;
    // Second guess linearizing the polynomial: p[0] + p[1] * x = p => x = (p - p[0]) / p[1]
    if self.p.len() > 1 && self.p[1] != 0.0 {
      r = (p - self.p[0]) / self.p[1];
    }
    // Solve by newton's method
    // If pathological, solve by the hybrid method
    self.newton_solve(p, r)
      .unwrap_or_else(|| self.hybrid_solve(p))
  }

  fn newton_solve(&self, p: f64, first_guess: f64) -> Option<f64> {
    // Newton's method
    let mut r = first_guess;
    let mut f = self.polynomial(r) - p;
    let mut i = 0;
    while f.abs() > self.eps && i < self.n_iter {
      r -= f / self.dpolynomial(r);
      f = self.polynomial(r) - p;
      i += 1;
    }
    // Deal with the result
    if f.abs() <= self.eps { // Algorithm converges
      if self.ang_dist.contains(&r) { // Result in domain of validity
        Some(r)
      } else if (*self.ang_dist.start() - self.eps..*self.ang_dist.start()).contains(&r) { // Valid at -epsilon
        Some(*self.ang_dist.start())
      } else if (*self.ang_dist.end()..=*self.ang_dist.end() + self.eps).contains(&r) {  // Valid at +epsilon
        Some(*self.ang_dist.end())
      } else { // Out of the domain of validity, find a better solution starting by dichotomy (slower)
        None
      }
    } else { // Algorithm did not converged
      if self.ang_dist.contains(&r) { // Result in domain of validity
        // Could be pathological, but I assume it is not, and problem due to a too small iteration
        // number or a too small epsilon.
        Some(r)
      } else { // Pathological situation
        None
      }
    }
  }

  /// Star by a dichotomy till max - min > step
  /// Then re-perform a Newton.
  fn hybrid_solve(&self, p: f64) -> f64 {
    // Make the first guess by dichotomy.
    let mut min = *self.ang_dist.start();
    let mut max = *self.ang_dist.end();
    while max - min > self.step {
      let med = 0.5 * (max + min);
      let med_val = self.polynomial(med);
      if med_val < p {
        min = med;
      } else {
        max = med;
      }
    }
    // Solve by Newton method starting at the center of the reduced domain.
    // If still pathological, solve by dichotomy.
   self.newton_solve(p, 0.5 * (max + min))
     .unwrap_or_else(
       || self.dichotomy_solve(p, min, max)
     )
  }

  fn dichotomy_solve(&self, y: f64, mut xmin: f64, mut xmax: f64) -> f64 {
    // Assuming a range of [0, pi], an epsilon of 5e-9, then the dichotomy will need
    // (pi - 0) / 2^n = 5e-9 => n = ln(pi/5e-9)/ln(2) => n = 30 iterations
    debug_assert!(xmin <= xmax);
    let mut min_val = self.polynomial(xmin);
    let mut max_val = self.polynomial(xmax);
    debug_assert!(min_val <= max_val);
    while max_val - min_val > self.eps {
      let xmed = 0.5 * (xmax + xmin);
      let med_val = self.polynomial(xmed);
      if med_val < y {
        xmin = xmed;
        min_val = med_val;
      } else {
        xmax = xmed;
        max_val = med_val;
      }
    }
    0.5 * (xmax + xmin)
  }

  fn polynomial(&self, x: f64) -> f64 {
    polynomial(x, &self.p)
  }
  
  fn dpolynomial(&self, x: f64) -> f64 {
    dpolynomial(x, &self.p)
  }
  
}


impl CanonicalProjection for Zpn {

  const NAME: &'static str = "Zenithal polynomlial";
  const WCS_NAME: &'static str = "ZPN";

  fn bounds(&self) -> &ProjBounds {
    &self.proj_bounds
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let r = (xyz.y.pow2() + xyz.z.pow2()).sqrt();
    // Compute the angular distance from the projection center.
    if r == 0.0 { // => y == 0 and z == 0
      let a = if xyz.x > 0.0 { // => x == 1 => angular distance = 0
        debug_assert!((1.0 - xyz.x).abs() < self.eps);
        self.p[0]
      } else { // x == -1 => angular distance = PI
        debug_assert!((1.0 + xyz.x).abs() < self.eps);
        self.polynomial(PI)
      };
      // Annulus of radius a
      Some(ProjXY::new(a, 0.0)) 
    } else {
      let a = if xyz.x > 0.0 {
        r.asin()
      } else {
        xyz.x.acos()
      };
      // Check the domain of validity
      if self.ang_dist.contains(&a) {
        // Compute projection
        let a = self.polynomial(a);
        if a >= 0.0 {
          // a /= r; we do not do this because of precision when r, y and z are smalls
          Some(ProjXY::new(a * (xyz.y / r), a * (xyz.z / r)))
        } else {
          // Problem in the study of the domain of validity!
          None
        }
      } else {
        // Out of the validity domain
        None
      }
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    // All positions valid, just opposite pole
    let r = (pos.x.pow2() + pos.y.pow2()).sqrt();
    // Check domain of validity
    if r == 0.0 {
      Some(XYZ::new(1.0, 0.0, 0.0))
    } else if self.euc_dist.contains(&r) {
      // Solve equation: polynomial(angular distance) = r
      let a = self.solve(r); // cos(a) = x => y^2+z^2 = 1 - cos^2(a)
      let x = a.cos();
      let a = (1.0 - x.pow2()).sqrt();
      // r /= a; => we do not do this because of precision when r, y and z are smalls
      Some(XYZ::new(x, a * (pos.x / r), a * (pos.y / r)))
    } else {
      None
    }
  }
}


/// Compute:
///   `p0 + x(p1 + x(p2 + x(p3 + x p4)))`
/// which is equivalent to:
///   `p0 + p1 x + p2 x^2 + p3 x^3 + p4 x^4`
fn polynomial(x: f64, p: &[f64]) -> f64 {
  let mut it = p.iter().rev();
  let init = it.next().cloned().unwrap_or(0.0);
  it.fold(init, move |acc, coeff| acc * x + *coeff)
}

/// Compute:
///   `p1 + x(2 p2 + x(3 p3 + x 4 p4)`
/// which is equivalent to
///   `p1 + 2 p2 x + 3 p3 x^2 + 4 p4 x^3`
fn dpolynomial(x: f64, p: &[f64]) -> f64 {
  let mut it = p.iter().enumerate().skip(1).rev();
  let init = it.next().map(|(i, coeff)| (i as f64) * coeff).unwrap_or(0.0);
  it.fold(init, move |acc, (i, coeff)| acc * x + (i as f64) * *coeff)
}


#[cfg(test)]
mod tests {
  use crate::{LonLat, Projection};
  use super::*;

  #[test]
  fn test_poly_1() {
    let coeffs = vec![0.0, 1.0, 0.0, -50.0];
    let zpn = Zpn::from_params(coeffs).unwrap();
    println!("{:?}", &zpn);
    let p0 = zpn.polynomial(0.0,);
    let p1 = zpn.polynomial(0.1);
    let p2 = zpn.polynomial(0.2);
    assert_eq!(p0, 0.0);
    assert_eq!(p1, 0.05);
    assert_eq!(p2, -0.2);
    assert_eq!(zpn.proj_lonlat(&LonLat::new(0.0, 0.0)), Some(ProjXY::new(0.0, 0.0)));
    assert_eq!(zpn.proj_lonlat(&LonLat::new(1.5e-3, 1.0e-3)), Some(ProjXY::new(0.0014997557501373168, 9.99837874976561E-4)));
    assert_eq!(zpn.unproj_lonlat(&ProjXY::new(0.0, 0.0)), Some(LonLat::new(0.0, 0.0)));
    assert_eq!(zpn.unproj_lonlat(&ProjXY::new(0.0014997557501373168, 9.99837874976561E-4)), Some(LonLat::new(0.0014999999999823648, 0.0009999999999882431)));
  }
  #[test]
  fn test_dpolyv1() {
    let coeffs = vec![0.0, 1.0, 0.0, -50.0];
    let zpn = Zpn::from_params(coeffs).unwrap();
    let p0 = zpn.dpolynomial(0.0);
    let p1 = zpn.dpolynomial(0.1);
    let p2 = zpn.dpolynomial(0.2);
    assert_eq!(p0, 1.0);
    assert_eq!(p1, -0.5);
    assert_eq!(p2, -5.0);
  }
  
  #[test]
  fn test_poly_2() {
    let coeffs = vec![0.050, 0.975, -0.807, 0.337, -0.065, 0.010, 0.003, -0.001];
    let zpn = Zpn::from_params(coeffs).unwrap();
    println!("{:?}", &zpn);
    let p0 = zpn.polynomial(0.0);
    let p1 = zpn.polynomial(0.1);
    let p2 = zpn.polynomial(0.2);
    assert_eq!(p0, 0.05);
    assert_eq!(p1, 0.1397606029);
    assert_eq!(p2, 0.2153153792);
    assert_eq!(zpn.proj_lonlat(&LonLat::new(0.0, 0.0)), Some(ProjXY::new(0.05, 0.0)));
    assert_eq!(zpn.proj_lonlat(&LonLat::new(1.5e-3, 1.0e-3)), Some(ProjXY::new(0.04306282454554517, 0.028708570032263056)));
    assert_eq!(zpn.unproj_lonlat(&ProjXY::new(0.05, 0.0)), Some(LonLat::new(0.0, 0.0)));
    assert_eq!(zpn.unproj_lonlat(&ProjXY::new(0.04306282454554517, 0.028708570032263056)), Some(LonLat::new(0.0014999999950119666,9.99999996674649E-4)));

  }
  
  #[test]
  fn test_dpoly_2() {
    let coeffs = vec![0.050, 0.975, -0.807, 0.337, -0.065, 0.010, 0.003, -0.001];
    let zpn = Zpn::from_params(coeffs).unwrap();
    let p0 = zpn.dpolynomial(0.0);
    let p1 = zpn.dpolynomial(0.1);
    let p2 = zpn.dpolynomial(0.2);
    assert_eq!(p0, 0.975);
    assert_eq!(p1, 0.8234551729999999);
    assert_eq!(p2, 0.690645312);
  }
}
