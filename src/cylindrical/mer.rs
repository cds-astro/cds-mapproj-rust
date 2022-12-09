//! Mercator projection.

use std::f64::consts::PI;
use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Mercator projection.
pub struct Mer;

impl Default for Mer {
  fn default() -> Self {
    Self::new()
  }
}

impl Mer {
  pub fn new() -> Self { Self }
}

impl CanonicalProjection for Mer {

  const NAME: &'static str = "Mercator";
  const WCS_NAME: &'static str = "MER";
  
  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-PI..=PI),
      None
    );
    &PROJ_BOUNDS
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if -1.0 < xyz.z && xyz.z < 1.0 {
      Some(ProjXY::new(xyz.y.atan2(xyz.x), xyz.z.atanh()))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    if (-PI..PI).contains(&pos.x) {
      let (sinl, cosl) = pos.x.sin_cos();
      let z = pos.y.tanh();
      let r = (1.0 - z.pow2()).sqrt(); // = cos(asin(z))
      Some(XYZ::new(r * cosl, r * sinl, z)) 
    } else {
      None
    }
  }
}