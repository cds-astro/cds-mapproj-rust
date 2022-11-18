//! Plate Carre projection.
use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

/// Plate Carre projection.
pub struct Car;

impl Default for Car {
  fn default() -> Self {
    Self::new()
  }
}

impl Car {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Car {

  const NAME: &'static str = "Plate Carre";
  const WCS_NAME: &'static str = "CAR";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // We do not use lat = asin(z) for precision purpose
    let r2 = xyz.x.pow2() + xyz.y.pow2();
    let lat = xyz.z.atan2(r2.sqrt());
    Some(ProjXY::new(xyz.y.atan2(xyz.x), lat))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let (slon, clon) = pos.x.sin_cos();
    let (slat, clat) = pos.y.sin_cos();
    Some(XYZ::new(clat * clon, clat * slon, slat))
  }
}