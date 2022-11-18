//! Cylindrical perspective projection.

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

/// Cylindrical perspective projection.
pub struct Cyp {
  /// `mu` parameter.
  mu: f64,
  /// `lambda` parameter.
  lambda: f64,
  /// re-computed value `mu + lambda`
  lpm: f64,
}

impl Default for Cyp {
  fn default() -> Self {
    Self::new()
  }
}

impl Cyp {
  
  pub fn new() -> Self {
    Self::from_params(1.0, 0.5 * 2.0_f64.sqrt())
  }
  
  /// # Params
  /// * `mu`: keyword `PVi_1a`.
  /// * `lambda`: keyword `PVi_2a`.
  pub fn from_params(mu: f64, lambda: f64) -> Self {
    Self {
      mu, lambda, lpm: mu + lambda
    }
  }
  
}

impl CanonicalProjection for Cyp {

  const NAME: &'static str = "Cylindrical perspective";
  const WCS_NAME: &'static str = "CYP";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let r = (xyz.x.pow2() + xyz.y.pow2()).sqrt(); // more accurate than sqrt(1 - z^2)
    Some(ProjXY::new(
      self.lambda * xyz.y.atan2(xyz.x),
      xyz.z * self.lpm / (self.mu + r)
    ))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let (sinl, cosl) = (pos.x / self.lambda).sin_cos();
    let nu = pos.y / self.lpm;
    let t = (self.mu * nu) / (1.0 + nu.pow2()).sqrt();
    let sqrt_1_p_nu2 = (1.0 + nu.pow2()).sqrt();
    let sqrt_1_m_t2 = (1.0 - t.pow2()).sqrt();
    let cosb = (sqrt_1_m_t2 - nu * t) / sqrt_1_p_nu2;
    let sinb = (nu * sqrt_1_m_t2 + t) / sqrt_1_p_nu2;
    Some(XYZ::new(cosb * cosl, cosb * sinl, sinb))
  }
}
