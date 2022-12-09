//! Samson-Flamsteed projection.
use std::f64::consts::PI;
use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};
use crate::math::HALF_PI;

/// Samson-Flamsteed projection.
pub struct Sfl;

impl Default for Sfl {
  fn default() -> Self {
    Self::new()
  }
}

impl Sfl {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Sfl {
  
  // Also called:
  // * Global Sinusoidal
  // * Mercator equal-area
  // * Mercator-Samson
  // * Sanson's
  const NAME: &'static str = "Samson-Flamsteed";
  const WCS_NAME: &'static str = "SFL";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-PI..=PI),
      Some(-HALF_PI..=HALF_PI)
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // We use x^2 + y^2 instead of 1 - z^2 for numerical precision on small distances
    let r = (xyz.x.pow2() + xyz.y.pow2()).sqrt();
    let lat = xyz.z.atan2(r); // could have use z.asin(), but atan2 is more accurate
    Some(ProjXY::new(
      xyz.y.atan2(xyz.x) * r, 
      lat
    ))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    if (-HALF_PI..=HALF_PI).contains(&pos.y) {
      let (z, r) = pos.y.sin_cos();
      // let z = pos.y.sin();
      // let r = (1.0 - z.pow2()).sqrt(); // = cos^2(Y);
      let l = if r == 0.0 { 0.0 } else { pos.x / r };
      if (-PI..=PI).contains(&l) {
        let (sinl, cosl) = l.sin_cos();
        Some(XYZ::new(r * cosl, r * sinl, z))
      } else {
        None
      }
    } else {
      None
    }
  }
}
