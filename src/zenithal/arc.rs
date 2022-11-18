//! Zenithal (or azimuthal) equidistant projection.

use std::f64::consts::PI;
use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

/// Zenithal (or azimuthal) equidistant projection.
pub struct Arc;

impl Default for Arc {
  fn default() -> Self {
    Self::new()
  }
}

impl Arc {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Arc {

  const NAME: &'static str = "Zenithal (or azimuthal) equidistant";
  const WCS_NAME: &'static str = "ARC";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if xyz.x > -1.0 {
      // Distance in the Euclidean plane (yz)
      // Angular distance is acos(x), but for small separation, asin(r) is more accurate.
      let r = (xyz.y.pow2() + xyz.z.pow2()).sqrt();
      let r = if xyz.x > 0.0 { // Angular distance < PI/2, angular distance = asin(r)
        r.asinc()
      } else { // Angular distance > PI/2, angular distance = acos(x)
        xyz.x.acos() / r
      };
      Some(ProjXY::new(xyz.y * r, xyz.z * r))
    } else {
      Some(ProjXY::new(PI, 0.0))
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    // r <= pi
    let r = (pos.x.pow2() + pos.y.pow2()).sqrt();
    if r <= PI {
      let x = r.cos();
      let r = r.sinc();
      Some(XYZ::new(x, pos.x * r, pos.y * r))
    } else {
      None
    }
  }
}

