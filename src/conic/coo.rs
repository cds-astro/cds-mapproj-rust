//! Conic Orthomorphic projection.

use std::f64::consts::PI;

use crate::{CustomFloat, CanonicalProjection, ProjXY, XYZ, math::HALF_PI, conic::Conic, ProjBounds};

/// Conic Orthomorphic projection.
pub struct Coo {
  conic: Conic,
  c: f64,
  one_over_c: f64,
  y0: f64,
  psi: f64,
}

impl Default for Coo {
  fn default() -> Self {
    Self::new()
  }
}

impl Coo {

  // default theta1 = theta2 = 45 deg
  pub fn new() -> Self {
    Self::from_params(HALF_PI.half(), 0.0)
  }

  pub fn from_params(theta_a: f64, nu: f64) -> Self {
    let conic = Conic::from_params(theta_a, nu);
    let cos_t1 = conic.theta1.cos();
    let tan_ft1 = (HALF_PI - conic.theta1).half().tan();
    let c = if conic.nu == 0.0 {
      debug_assert_eq!(conic.theta1, conic.theta2);
      conic.theta1.sin()
    } else {
      (conic.theta2.cos() / cos_t1).ln() / ((HALF_PI - conic.theta2).half().tan() / tan_ft1).ln()
    };
    let psi = cos_t1 / (c * tan_ft1.powf(c));
    let y0 = psi * (HALF_PI - conic.ta).half().tan().powf(c);
    let one_over_c = 1.0 / c;
    Self {
      conic,
      c,
      one_over_c,
      y0,
      psi
    }
  }
}


impl CanonicalProjection for Coo {

  const NAME: &'static str = "Conic Orthomorphic";
  const WCS_NAME: &'static str = "Coo";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      None,
      None
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let lon = xyz.y.atan2(xyz.x);
    // Use something else than asin?
    let r = self.psi * (HALF_PI - xyz.z.asin()).half().tan().powf(self.c);
    let (sinc, cosc) = (self.c * lon).sin_cos();
    Some(ProjXY::new(r * sinc, self.y0 - r * cosc))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    const EPS: f64 = 1.0e-14;
    let x2d = pos.x;
    let y2d = self.y0 - pos.y;
    let r2 = x2d.pow2() + y2d.pow2();
    let r = if self.conic.negative_ta { -(r2.sqrt()) } else { r2.sqrt() };
    let lon =  (x2d / r).atan2(y2d / r) / self.c; // / r important because of its sign
    if (-PI - EPS..PI + EPS).contains(&lon) {
      let lat = HALF_PI - (r / self.psi).powf(self.one_over_c).atan().twice();
      let (sinb, cosb) = lat.sin_cos();
      let (sinl, cosl) = lon.sin_cos();
      Some(XYZ::new(cosb * cosl, cosb * sinl, sinb))
    } else {
      None
    }
  }
}
