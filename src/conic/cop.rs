//! Conic perspective projection.

use std::f64::consts::PI;

use crate::{CustomFloat, CanonicalProjection, ProjXY, XYZ, math::HALF_PI, conic::Conic, ProjBounds};

/// Conic perspective projection.
pub struct Cop {
  conic: Conic,
  c: f64,
  y0: f64,
  cos_nu: f64,
  cos_ta: f64,
  sin_ta: f64,
  tan_ta: f64,
  cotan_ta: f64,
  z_min: f64,
  z_max: f64,
}

impl Default for Cop {
  fn default() -> Self {
    Self::new()
  }
}

impl Cop {

  pub fn new() -> Self {
    Self::from_params(0.0, 0.0)
  }

  pub fn from_params(theta_a: f64, nu: f64) -> Self {
    let conic = Conic::from_params(theta_a, nu);
    let (sin_ta, cos_ta) = conic.ta.sin_cos();
    let tan_ta = conic.ta.tan();
    let cotan_ta = 1.0 / tan_ta;
    let c = sin_ta;
    let cos_nu = conic.nu.cos();
    let y0 = cos_nu * cotan_ta;
    let (z_min, z_max) = if conic.ta >= 0.0 {
      ((-HALF_PI + conic.ta).sin(), 1.0)
    } else {
      (-1.0, (HALF_PI + conic.ta).sin())
    };
    Self {
      conic, c, y0, cos_nu,
      cos_ta, sin_ta, tan_ta, cotan_ta,
      z_min, z_max
    }
  }
  
  fn is_lat_in_domain_of_validity(&self, z: f64) -> bool {
    self.z_min < z && z < self.z_max
  }
}


impl CanonicalProjection for Cop {

  const NAME: &'static str = "Conic perspective";
  const WCS_NAME: &'static str = "Cop";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      None,
      None
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if self.is_lat_in_domain_of_validity(xyz.z) {
      let lon = xyz.y.atan2(xyz.z);
      // r(dec) depends on tan(dec - theta_a) 
      // => -pi/2 < dec - theta_a < pi/2 => -pi/2 + theta_a < dec < pi/2 + theta
      let r = xyz.z / (1.0 - xyz.z.pow2()).sqrt();
      if r.is_finite() {
        let r = (r - self.tan_ta) / (1.0 + r * self.tan_ta);  // to avoid computing tan(dec - ta)
        let r = self.cos_nu * (self.cotan_ta - r);
        let (sinc, cosc) = (self.c * lon).sin_cos();
        Some(ProjXY::new(r * sinc, self.y0 - r * cosc))
      } else { // r == NaN | +Infinity | -Infinity
        debug_assert!(is_close_to_one(xyz.z));
        Some(ProjXY::new(0.0, self.y0))
      }
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    const EPS: f64 = 1.0e-14;
    let x2d = pos.x;
    let y2d = self.y0 - pos.y;
    let r2 = x2d.pow2() + y2d.pow2();
    let r = if self.conic.negative_ta { -(r2.sqrt()) } else { r2.sqrt() };
    let lon =  y2d.atan2(x2d) / self.c; // no need to divide both y2d and x2d by r
    if (-PI - EPS..PI + EPS).contains(&lon) {
      // dec = this.ta + atan(this.cotanta - r / this.cosnu);
      // sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
      // cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
      // cos(atan(x)) = 1 / sqrt(1+x^2)
      // sin(atan(x)) = x / sqrt(1+x^2)
      let r = self.cotan_ta - r / self.cos_nu;
      let d = (1.0 + r.pow2()).sqrt();
      let (sinl, cosl) = lon.sin_cos();
      let cosb = (self.cos_ta - r * self.sin_ta) / d;
      let sinb = (self.sin_ta + r * self.cos_ta) / d;
      Some(XYZ::new(cosb * cosl, cosb * sinl, sinb))
    } else {
      None
    }
  }
}

fn is_close_to_one(x: f64) -> bool {
  const EPS: f64 = 1.0e-15;
  (1.0 - x).abs() < EPS
}
