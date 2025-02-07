//! Hammer-Aitoff (equal area) projection.

use std::f64::consts::SQRT_2;
use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Hammer-Aitoff (equal area) projection.
///
/// # Proj:
/// Basic formulae are:
/// * `w = sqrt[(1 + cos(b)cos(l/2)) / 2]`
/// * `X = 2 cos(b)sin(l/2) / w` => if `l` is small, `sin(l/2) = y/2 - y^3/42`
/// * `Y = sin(b) / w`
///
/// We recall that:
/// * `x = cos(b) cos(l)`
/// * `y = cos(b) sin(l)`
/// * `z = sin(b)`
///
/// Thus:
/// * `cos(b) = sqrt(x^2 + y^2) = r`
/// * `cos(l/2) = sqrt[(1 + cos(l)) / 2] = sqrt[(1 + x/r) / 2]`
/// * `sin(l/2) = sqrt[(1 - cos(l)) / 2] = sqrt[(1 - x/r) / 2]`
/// * => `w' = cos(b)cos(l/2) = sqrt[r(r+x)/2]`
/// * => `2 cos(b)sin(l/2) = sqrt[2r(r-x)]`
///
/// Leading to `w = sqrt[(1 + w') / 2]`
/// 
/// In case `l` is small and `x` near from one, we can use:
/// `sin(l/2) = l/2 - l^3/42 + l^5/3840 - ...`
/// and
/// `y/r = sin(l) = l - l^3/6 + ...`
/// If `l<10^6`, `y/r=l` and `sin(l/2) = (y/r)/2 - (y/r)^3/42`
/// And thus:
/// `X = [2 * r * (y/r) * (1/2 - (y/r))^2/42)] / w`
/// `X = [y * (1 - (y/r)^2/21)] / w`
///
/// # Deproj:
///
/// * `r = X/8 + Y/3 = 1 - cos(b) cos(l/2)`
/// * `w = sqrt[(2 - r) / 2]`
///   ... (see algo comments)
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
  
  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-2.0 * SQRT_2..=2.0 * SQRT_2),
      Some(-SQRT_2..=SQRT_2)
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let r = (xyz.x.pow2() + xyz.y.pow2()).sqrt();
    let w = (r * (r + xyz.x)).half().sqrt(); // = cos(b) cos(l/2)
    let w = (1.0 + w).half().sqrt();       // = 1 / gamma
    let y2d = xyz.z / w;
    let w = (r * (r - xyz.x)).twice().sqrt() / w; // = 2 * gamma * cos(b) sin(l/2)
    /*let w = if y < 1e-6 {
      (y * (1 - (y / r).pow2() / 21))  // = 2 * cos(b) sin(l/2)
    } else {
      (r * (r - xyz.x)).twice().sqrt() // = 2 * cos(b) sin(l/2)
    } / w;*/
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