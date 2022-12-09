//! Lambert's zenithal (or azimuthal) equal area projection.

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Lambert's zenithal (or azimuthal) equal area projection.
pub struct Zea;

impl Default for Zea {
  fn default() -> Self {
    Self::new()
  }
}

impl Zea {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Zea {

  const NAME: &'static str = "Lambert's zenithal (or azimuthal) equal area";
  const WCS_NAME: &'static str = "ZEA";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-2.0..=2.0),
      Some(-2.0..=2.0)
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // Whole sphere, r <= 2 (equal area)
    let w = (0.5 + 0.5 * xyz.x).sqrt(); // <=> sqrt[(1 + x) / 2]
    if w > 0.0 {
      Some(ProjXY::new(xyz.y / w, xyz.z / w))
    } else {
      Some(ProjXY::new(2.0, 0.0))
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    // Whole sphere, r <= 2 (equal area)
    let r = 0.25 * (pos.x.pow2() + pos.y.pow2());
    if r <= 1.0 {
      let w = (1.0 - r).sqrt();
      Some(XYZ::new_renorming_if_necessary(1.0 - r.twice(),  pos.x * w,  pos.y * w))
    } else {
      None
    }
  }
}
