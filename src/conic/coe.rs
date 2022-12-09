//! Conic Equal Area projection.

use std::f64::consts::PI;

use crate::{CustomFloat, CanonicalProjection, ProjXY, XYZ, conic::Conic, ProjBounds};

/// Conic Equal Area projection.
pub struct Coe {
  conic: Conic,
  one_plus_sint1_sint2: f64,
  gamma: f64,
  c: f64,
  c2: f64,
  y0: f64,
  r2_min: f64,
  r2_max: f64,
  proj_bounds: ProjBounds,
}

impl Default for Coe {
  fn default() -> Self {
    Self::new()
  }
}

impl Coe {

  pub fn new() -> Self {
    Self::from_params(0.0, 0.0)
  }

  pub fn from_params(theta_a: f64, nu: f64) -> Self {
    let conic = Conic::from_params(theta_a, nu);
    let sin_t1 = conic.theta1.sin();
    let sin_t2 = conic.theta2.sin();
    let sin_ta = conic.ta.sin();
    let one_plus_sint1_sint2 = 1.0 + sin_t1 * sin_t2;
    let gamma = sin_t1 + sin_t2;
    let c = 0.5 * gamma;
    let c2 = c.pow2();
    let y0 = (one_plus_sint1_sint2 - gamma * sin_ta).sqrt() /  c;
    let (r2_min, r2_max) = if gamma >= 0.0 {
      (
        (one_plus_sint1_sint2 - gamma) /  c2,
        (one_plus_sint1_sint2 + gamma) / c2
      )
    } else {
      (
        (one_plus_sint1_sint2 + gamma) / c2,
        (one_plus_sint1_sint2 - gamma) / c2
      )
    };
    debug_assert!(r2_min <= r2_max);
    let r_max = r2_max.sqrt();
    Self {
      conic,
      one_plus_sint1_sint2, gamma, c, c2, y0,
      r2_min, r2_max,
      proj_bounds: ProjBounds::new(
        Some(-r_max..=r_max),
        Some(y0 - r_max..=y0 + r_max)
      )
    }
  }
}


impl CanonicalProjection for Coe {

  const NAME: &'static str = "Conic Equal Area";
  const WCS_NAME: &'static str = "Coe";

  fn bounds(&self) -> &ProjBounds {
    &self.proj_bounds
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let r = self.one_plus_sint1_sint2 - self.gamma * xyz.z;
    if r >= 0.0 {
      let r = r.sqrt() / self.c;
      let lon = xyz.y.atan2(xyz.x);
      let (sin, cos) = (self.c * lon).sin_cos();
      Some(ProjXY::new(r * sin, self.y0 - r * cos))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    const EPS: f64 = 1.0e-14;
    let x2d = pos.x;
    let y2d = self.y0 - pos.y;
    let r2 = x2d.pow2() + y2d.pow2();
    if (self.r2_min - EPS..self.r2_max + EPS).contains(&r2) {
      let r = if self.conic.negative_ta { -r2.sqrt() } else { r2.sqrt() };
      let lon =  y2d.atan2(x2d) / self.c; // no need to divide both y2d and x2d by r
      if (-PI - EPS..PI + EPS).contains(&lon) {
        let z = (self.one_plus_sint1_sint2 - self.c2 * r) / self.gamma;
        let r = (1.0 - z.pow2()).sqrt();
        let (sinl, cosl) = lon.sin_cos();
        Some(XYZ::new(r * cosl, r * sinl, z))
      } else {
        None
      }
    } else {
      None
    }
  }
}
