//! Orthographic projections.

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Orthographic projection.
pub struct Sin;

impl Default for Sin {
  fn default() -> Self {
    Self::new()
  }
}

impl Sin {
  pub fn new() -> Self {
    Self
  }
}

impl CanonicalProjection for Sin {
  
  const NAME: &'static str = "Orthographic";
  const WCS_NAME: &'static str = "SIN";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      Some(-1.0..=1.0),
      Some(-1.0..=1.0)
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if xyz.x >= 0.0 {
      Some(ProjXY::new(xyz.y, xyz.z))
    } else { // Back hemisphere
      None
    }
  }

  fn unproj(&self, xy: &ProjXY) -> Option<XYZ> {
    let r2 = xy.x.pow2() + xy.y.pow2();
    if r2 <= 1.0 {
      let x = (1.0 - r2).sqrt(); //  since x^2 + y^2 + z^2 = 1
      Some(XYZ::new(x, xy.x, xy.y))
    } else {
      None
    }
  }
}

/// Slant Orthographic projection.
pub struct SinSlant {
  xi: f64,
  eta: f64,
  xp: f64,
  yp: f64,
  zp: f64,
  tg2: f64,
  proj_bounds: ProjBounds,
}

impl SinSlant {
  
  /// # Params
  /// * `xi`: corresponds to `-PV1` in WCS
  /// * `eta`: corresponds to `PV2` in WCS
  /// # Remark
  /// if `xi = eta = 0`, use `Sin` instead of `SinSlant` 
  pub fn new(xi: f64, eta: f64) -> Self {
    let tg2 = xi.pow2() + eta.pow2();
    let tmp = (1.0 + tg2).sqrt();
    Self {
      xi, eta,
      xp: -1.0 / tmp,
      yp: -xi / tmp,
      zp: -eta / tmp,
      tg2,
      proj_bounds: ProjBounds::new(
        Some(-1.0..=1.0 + xi * 2.0),
        Some(-1.0..=1.0 + eta * 2.0)
      )
    }
  }
}

impl CanonicalProjection for SinSlant {

  const NAME: &'static str = "Slant orthographic";
  const WCS_NAME: &'static str = "SIN";

  fn bounds(&self) -> &ProjBounds {
    &self.proj_bounds
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let s = xyz.x * self.xp + xyz.y * self.yp + xyz.z * self.zp;
    if s <= 0.0 {
      let s = 1.0 - xyz.x;
      Some(ProjXY::new(xyz.y + self.xi * s, xyz.z + self.eta * s))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let x2d = pos.x;
    let y2d = pos.y;
    let r2 = x2d.pow2() + y2d.pow2();
    let s = self.xp + x2d * self.yp + y2d * self.yp;
    let v = (1.0 - self.xp * s).pow2()
      + (x2d - self.yp * s).pow2()
      + (y2d - self.zp * s).pow2();
    if v < 1.0 {
      let rp = self.xi * x2d + self.eta * y2d; // R'
      let a = 1.0 + self.tg2;
      let b = 2.0 * (rp - self.tg2);
      let c = r2 - 2.0 * rp + self.tg2 - 1.0;
      let x = (-b + (b.pow2() - 4.0 * a * c).sqrt()) / (2.0 * a);
      let tmp = 1.0 - x;
      Some(XYZ::new(x, x2d - self.xi * tmp, y2d - self.eta * tmp))
    } else {
      None
    }
  }
}
