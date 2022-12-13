//! Conic Equidistant projection.
 
use std::f64::consts::PI;

use crate::{CustomFloat, CanonicalProjection, ProjXY, XYZ, math::HALF_PI, conic::Conic, ProjBounds};

/// Conic Equidistant projection.
#[derive(Debug, Clone)]
pub struct Cod {
  conic: Conic,
  c: f64,
  y0: f64,
  r2_min: f64,
  r2_max: f64,
  ta_plus_y0: f64,
  proj_bounds: ProjBounds,
}

impl Default for Cod {
  fn default() -> Self {
    Self::new()
  }
}

impl Cod {

  // default theta1 = theta2 = 45 deg
  pub fn new() -> Self {
    Self::from_params(HALF_PI.half(), 0.0)
  }

  pub fn from_params(theta_a: f64, nu: f64) -> Self {
    let conic = Conic::from_params(theta_a, nu);
    let sin_ta = conic.ta.sin();
    let cot_ta = 1.0 / conic.ta.tan();
    let (c, y0) = if conic.nu == 0.0 {
      debug_assert_eq!(conic.theta1, conic.theta2);
      (sin_ta, cot_ta)
    } else {
      let sin_nu = conic.nu.sin();
      ((sin_ta * sin_nu) / conic.nu, conic.nu * cot_ta / conic.nu.tan())
    };
    let ta_plus_y0 = conic.ta + y0;
    let (r_min, r_max) = if ta_plus_y0 >= 0.0 {
      (ta_plus_y0 - HALF_PI, ta_plus_y0 + HALF_PI)
    } else {
      (ta_plus_y0 + HALF_PI, ta_plus_y0 - HALF_PI)
    };
    let yrange = if conic.negative_ta {
      Some(y0 - r_max * (PI * c).cos().abs()..=y0 + r_max)
    } else {
      Some(y0 - r_max..=y0 + r_max * (PI * c).cos().abs())
    };
    Self {
      conic,
      c,
      y0,
      r2_min: r_min.pow2(),
      r2_max: r_max.pow2(),
      ta_plus_y0,
      proj_bounds: ProjBounds::new(
        Some(-r_max..=r_max),
        yrange
      )
    }
  }
}


impl CanonicalProjection for Cod {

  const NAME: &'static str = "Conic Equidistant";
  const WCS_NAME: &'static str = "Cod";

  fn bounds(&self) -> &ProjBounds {
    &self.proj_bounds
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let lon = xyz.y.atan2(xyz.x);
    // more computations but ore accurate than r = self.ta_plus_y0 - z.asin()
    let r = (xyz.x.pow2() + xyz.y.pow2()).sqrt();
    let lat = xyz.z.atan2(r);
    let r = self.ta_plus_y0 - lat;
    let (sinc, cosc) = (self.c * lon).sin_cos();
    Some(ProjXY::new(r * sinc, self.y0 - r * cosc))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    const EPS: f64 = 1.0e-14;
    let x2d = pos.x;
    let y2d = self.y0 - pos.y;
    let r2 = x2d.pow2() + y2d.pow2();
    if (self.r2_min..=self.r2_max).contains(&r2) {
      let r = if self.conic.negative_ta { -(r2.sqrt()) } else { r2.sqrt() };
      let lon =  (x2d / r).atan2(y2d / r) / self.c; // / r important because of its sign
      if (-PI - EPS..PI + EPS).contains(&lon) {
        let (sinb, cosb) = (self.ta_plus_y0 - r).sin_cos();
        let (sinl, cosl) = lon.sin_cos();
        Some(XYZ::new(cosb * cosl, cosb * sinl, sinb))
      } else {
        None
      }
    } else {
      None
    }
  }
}
