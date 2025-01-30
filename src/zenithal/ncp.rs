//! North Celestial Pole orthographic projection.
use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// North Celestial Pole orthographic projection.
pub struct Ncp;

impl Default for Ncp {
  fn default() -> Self {
    Self::new()
  }
}

impl Ncp {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Ncp {

  const NAME: &'static str = "North Celestial Pole orthographic";
  const WCS_NAME: &'static str = "NCP";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-1.0..=1.0),
      Some(-1.0..=1.0)
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if xyz.x >= 0.0 { // Front hemisphere
      Some(ProjXY::new(xyz.y, xyz.z))
    } else {
      None // Back hemisphere
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let n2 = pos.x.pow2() + pos.y.pow2(); // = 1 - (x^2 + y^2) = z^2
    if n2 <= 1.0 {
      Some(XYZ::new((1.0 - n2).sqrt(), pos.x, pos.y))
    } else {
      None
    }
  }
}
