//! Hammer-Aitoff (equal area) projection.

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

/// Hammer-Aitoff (equal area) projection.
pub struct Ait;

impl Default for Ait {
  fn default() -> Self {
    Self::new()
  }
}

impl Ait {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Ait {

  const NAME: &'static str = "Hammer-Aitoff (equal area)";
  const WCS_NAME: &'static str = "AIT";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let r = (xyz.x.pow2() + xyz.y.pow2()).sqrt();
    let w = (r * (r + xyz.x)).half().sqrt(); // = cos(b) cos(l/2)
    let w = (1.0 + w).half().sqrt();       // = 1 / gamma
    let y2d = xyz.z / w;
    let w = (r * (r - xyz.x)).twice().sqrt() / w; // = 2 * gamma * cos(b) sin(l/2)
    let x2d = if xyz.y < 0.0 { -w } else { w };
    Some(ProjXY::new(x2d, y2d))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    // Ellipse, dimensions sqrt(2) x 2.sqrt(2)
    let r = 0.125 * pos.x.pow2() + pos.y.pow2().half(); //  = 1 - cos(b) cos(l/2)
    if r > 1.0  {
      None
    } else {
      let mut x = 1.0 - r; // cos(b) cos(l/2)
      let mut w = (1.0 - r.half()).sqrt(); // sqrt(HALF * (1 + x)) ;  //  = Z = sqrt[ (1 + cos(b) cos(l/2)) / 2]
      let mut y = pos.x.half() * w; // cos(b) sin(l/2)
      let z = pos.y * w; // z
      // Convert from Cartesian (l/2, b) to Cartesian (l, b) 
      let r = (x.pow2() + y.pow2()).sqrt();  // cos(b)
      if r > 0.0 {
        w = x;
        x = (w.pow2() - y.pow2()) / r; // cos(b) cos(l)
        y = (w * y / r).twice(); // cos(b) sin(l)
      }
      Some(XYZ::new_renorming_if_necessary(x, y, z))
      // Some(XYZ::new(x, y, z))
    }
  }
}