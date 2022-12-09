//! Fisheye projection.

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

const D_MAX: f64 = 1.6580627893946132; // 95.0_f64.to_radians();

/// Fisheye projection.
pub struct Feye;

impl Default for Feye {
  fn default() -> Self {
    Self::new()
  }
}

impl Feye {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Feye {

  const NAME: &'static str = "Fisheye";
  const WCS_NAME: &'static str = "FEYE";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-D_MAX..=D_MAX),
      Some(-D_MAX..=D_MAX)
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // Distance in the Euclidean plane (yz).
    // Angular distance is acos(x), but for small separation, asin(r) is more accurate.
    let r = (xyz.y.pow2() + xyz.z.pow2()).sqrt();
    let r = if xyz.x > 0.0 { // Angular distance < PI/2, angular distance = asin(r)
      r.asinc()
    } else { // Angular distance > PI/2, angular distance = acos(x)
      xyz.x.acos() / r
    };
    if r <= D_MAX {
      Some(ProjXY::new(xyz.y * r, xyz.z * r))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    // r <= pi
    let r = (pos.x.pow2() + pos.y.pow2()).sqrt();
    if r <= D_MAX {
      let x = r.cos();
      let r = r.sinc();
      Some(XYZ::new(x, pos.x * r, pos.y * r))
    } else {
      None
    }
  }
}
