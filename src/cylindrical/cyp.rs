//! Cylindrical perspective projection.

use std::f64::consts::PI;
use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

/// Cylindrical perspective projection.
pub struct Cyp {
  /// `mu` parameter.
  mu: f64,
  /// `lambda` parameter.
  lambda: f64,
  /// pre-computed value `mu + lambda`
  lpm: f64,
  /// pre-compute `PI * lambda`
  pi_x_lambda: f64,
  /// pre-computed `1 + lambda / mu`
  lpm_over_mu: f64,
  /// pre-computed projection bounds
  proj_bounds: ProjBounds,
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
    let pi_x_lambda = PI * lambda;
    let lpm_over_mu = (mu + lambda) / mu;
    Self {
      mu, 
      lambda, 
      lpm: mu + lambda,
      pi_x_lambda,
      lpm_over_mu,
      proj_bounds: ProjBounds::new(
        Some(-pi_x_lambda..=pi_x_lambda),
        Some(-lpm_over_mu..=lpm_over_mu)
      )
    }
  }
  
}

impl CanonicalProjection for Cyp {

  const NAME: &'static str = "Cylindrical perspective";
  const WCS_NAME: &'static str = "CYP";

  fn bounds(&self) -> &ProjBounds {
    &self.proj_bounds
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let r = (xyz.x.pow2() + xyz.y.pow2()).sqrt(); // more accurate than sqrt(1 - z^2)
    Some(ProjXY::new(
      self.lambda * xyz.y.atan2(xyz.x),
      xyz.z * self.lpm / (self.mu + r)
    ))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    if (-self.pi_x_lambda..=self.pi_x_lambda).contains(&pos.x) 
      && (-self.lpm_over_mu..=self.lpm_over_mu).contains(&pos.y) {
      let (sinl, cosl) = (pos.x / self.lambda).sin_cos();
      let nu = pos.y / self.lpm;
      let t = (self.mu * nu) / (1.0 + nu.pow2()).sqrt();
      let sqrt_1_p_nu2 = (1.0 + nu.pow2()).sqrt();
      let sqrt_1_m_t2 = (1.0 - t.pow2()).sqrt();
      let cosb = (sqrt_1_m_t2 - nu * t) / sqrt_1_p_nu2;
      let sinb = (nu * sqrt_1_m_t2 + t) / sqrt_1_p_nu2;
      Some(XYZ::new(cosb * cosl, cosb * sinl, sinb))
    } else {
      None
    }
  }
}
