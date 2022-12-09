//! Parabolic projection.

use std::f64::consts::PI;
use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Parabolic projection.
pub struct Par;

impl Par {
  pub fn new() -> Self {
    Self
  }
}

impl Default for Par {
  fn default() -> Self {
    Self::new()
  }
}

impl CanonicalProjection for Par {

  const NAME: &'static str = "Parabolic";
  const WCS_NAME: &'static str = "PAR";
  
  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-PI..=PI),          // +- pi
      Some(-0.5_f64..=0.5_f64) // +-sin(pi/6)
    );
    &PROJ_BOUNDS
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // more computations, but more precise than lat = xyz.z.asin()
    let r2 = xyz.x.pow2() + xyz.y.pow2();
    let lat = xyz.z.atan2(r2.sqrt());
    // let lat = xyz.z.asin();
    Some(ProjXY::new(
      xyz.y.atan2(xyz.x) * ((lat / 1.5).cos().twice() - 1.0), 
      (lat / 3.0).sin()
    ))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    if (-0.5_f64..=0.5_f64).contains(&pos.y) {
      let (z, r) = (3.0 * pos.y.asin()).sin_cos(); // sin(b)
      let m = 1.0 - pos.y.twice().pow2(); // r = 0 => sin(b / 3) = 1/2 => b = +-90
      let l = if m == 0.0 { 0.0 } else { pos.x / m };
      if (-PI..=PI).contains(&l) {
        let (sinl, cosl) = l.sin_cos();
        // let r = (1.0 - z.pow2()).sqrt(); // cos(b), =0 if b = +-90 => z = +-1
        Some(XYZ::new(r * cosl, r * sinl, z))
      } else {
        None
      }
    } else {
      None
    }
  }
}
