//! Gnomonic projection.

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

/// Gnomonic projection.
pub struct Tan;

impl Default for Tan {
  fn default() -> Self {
    Self::new()
  }
}

impl Tan {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Tan {
  
  const NAME: &'static str = "Gnomonic";
  const WCS_NAME: &'static str = "TAN";
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if xyz.x > 0.0 { // EPSILON ??
      Some(ProjXY::new(xyz.y / xyz.x, xyz.z / xyz.x))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let x = 1_f64 / (1_f64 + pos.x.pow2() + pos.y.pow2()).sqrt();
    Some(XYZ::new(x, pos.x * x, pos.y * x))
  }
}
