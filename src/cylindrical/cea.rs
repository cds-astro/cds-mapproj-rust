//! Cylindrical equal area projection.

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

/// Cylindrical equal area projection.
/// With default value `lambda = 1`, this projection is a 
/// Lambert's Cylindrical or Lambert's Equal Area projection.
pub struct Cea {
  lambda: f64
}

impl Default for Cea {
  fn default() -> Self {
    Self::new()
  }
}

impl Cea {
  
  pub fn new() -> Self {
    Self::from_param(1.0)
  }
  
  /// # Params
  /// * `lambda`: correspond to the WCS `PVi_1a` parameter
  /// # Panics 
  /// * if `lambda = 0` or is not finite.
  pub fn from_param(lambda: f64) -> Self {
    assert!(lambda != 0.0 && lambda.is_finite());
    Self { lambda }
  }
  
}

impl CanonicalProjection for Cea {

  const NAME: &'static str = "Cylindrical equal area";
  const WCS_NAME: &'static str = "CEA";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    Some(ProjXY::new(xyz.y.atan2(xyz.x), xyz.z / self.lambda))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let (sinl, cosl) = pos.x.sin_cos();
    let sinb = self.lambda * pos.y; // = z
    let cosb = (1.0 - sinb.pow2()).sqrt(); // = sqrt(1 - z^2) = sqrt(x^2 + y^2) = sinb.asin().cos();
    Some(XYZ::new(cosb * cosl, cosb * sinl, sinb))
  }
}
