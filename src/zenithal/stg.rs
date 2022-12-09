//! Stereographic projection.

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Stereographic projection.
pub struct Stg;

impl Default for Stg {
  fn default() -> Self {
    Self::new()
  }
}

impl Stg {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Stg {

  const NAME: &'static str = "Stereographic";
  const WCS_NAME: &'static str = "STG";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      None,
      None
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // All positions are valid, but diverges at lat = -PI/2
    let w = (1.0 + xyz.x).half();
    if w > 0.0 {
      Some(ProjXY::new(xyz.y / w, xyz.z / w))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    // All positions valid, just opposite pole
    let r = 0.25 * (pos.x.pow2() + pos.y.pow2());
    let w = 1.0 + r;
    Some(XYZ::new((1.0 - r)  /w,  pos.x / w,  pos.y / w))
  }
}
