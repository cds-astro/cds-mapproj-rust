//! Zenithal perspective projection.

use std::f64::consts::PI;

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

static HALF_PI: f64 = 0.5 * PI;

/// Zenithal perspective projection.
pub struct Azp {
  /// WCS keyword PVi_1a.
  mu: f64,
  /// WCS keyword PVi_2a (but converted in radians).
  gamma: f64,
  // Pre computed quantities
  tg: f64, // tan(gamma)
  cg: f64, // cos(gamma)
  sg: f64, // sin(gamma)
  abs_mu: f64, // |mu|
  m_p_1: f64,  // mu + 1;
  sqrt_mu2_m_1: f64, // sqrt(mu^2 - 1);
  x_min: f64,        // mu == 0 ? 0 : -1 / this.mu;
}

impl Default for Azp {
  fn default() -> Self {
    Self::new()
  }
}

impl Azp {

  /// New AZP projection with default parameters:
  /// * `mu = 1.35`
  /// * `gamma = 0`
  pub fn new() -> Self {
    // Gnomonic             if mu = 0
    // Near-side perspetive if mu = 1.35
    // Clarke's (first)     if mu = 1.35
    // Clarke's (second)    if mu = 1.65
    // James'               if mu = 1.367
    // La Hire's            if mu = 1.71
    // Approx equidist zenith if mu = 1.7519
    // Aprrox equ-area zenith if mu = 2.4142
    Self::from_params(1.35, 0.0)
  }

  /// New AZP projection with custom parameters:
  /// # Paras
  /// * `mu`, WCS parameter PVi_1a
  /// * `gamma`, WCS parameter PVi_2a (but converted in radians)
  /// # Pancis
  /// * if `gamma` not in `[-pi/2, pi/2]`
  pub fn from_params(mu: f64, gamma: f64) -> Self {
    assert!((-HALF_PI..=HALF_PI).contains(&gamma));
    Self {
      mu,
      gamma,
      tg: gamma.tan(),
      cg: gamma.cos(),
      sg: gamma.sin(),
      abs_mu: mu.abs(),
      m_p_1: mu + 1.0,
      sqrt_mu2_m_1: (mu.pow2() - 1.0).sqrt(),
      x_min: if mu == 0.0 { 0.0 } else { -1.0 / mu },
    }
  }

  /// Get the value of the `mu` parameter.
  pub fn mu(&self) -> f64 {
    self.mu
  }

  /// Get the value of the `gamma` parameter.
  pub fn gamma(&self) -> f64 {
    self.gamma
  }
  
  
}

impl CanonicalProjection for Azp {

  const NAME: &'static str = "Zenithal perspective";
  const WCS_NAME: &'static str = "AZP";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let wx = self.mu + xyz.x;
    let wy = self.cg * wx - xyz.z * self.sg;
    let wx = wx - xyz.z * self.tg;
    if xyz.x < self.x_min || wx == 0.0 || wy == 0.0
      || (self.abs_mu < 1.0 && (xyz.x + self.mu) * self.cg <= xyz.z * self.sg) {
      None
    } else {
      Some(ProjXY::new((xyz.y * self.m_p_1) / wx, (xyz.z * self.m_p_1) / wy))
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let y2d_cg = pos.y * self.cg;
    let big_r = (pos.x.pow2() + y2d_cg.pow2()).sqrt();
    let w = self.m_p_1 + (pos.y * self.sg);
    if big_r == 0.0 {
      // Center of the projection
      Some(XYZ::new(1.0, 0.0, 0.0))
    } else if w == 0.0 {
      // Special case, happens only if |mu| < 1
      debug_assert!(self.abs_mu <= 1.0);
      let w = (1.0 - self.mu.pow2()).sqrt() / big_r;
      Some(XYZ::new(-self.mu, pos.x * w, y2d_cg * w))
    } else if (self.abs_mu > 1.0 && (big_r * self.sqrt_mu2_m_1) > w)
     || (self.abs_mu < 1.0 && big_r > 1.0e9) {
      None
    } else {
      let w = big_r / w;
      let w1 = (w.pow2() + 1.0).sqrt();
      let w2 = (self.mu * w) / w1;
      let mut w3 = (1.0 - w2.pow2()).sqrt();
      let mut sina = (w2 + w3 * w) / w1;
      if self.abs_mu < 1.0 && sina < 0.0 {
        w3 = -w3;
        sina = (w2 + w3 * w) / w1;
      }
      let cosa = (w3 - w2 * w) / w1;
      let w = sina / big_r;
      Some(XYZ::new_renorming_if_necessary(cosa, pos.x * w, y2d_cg * w))
    }
  }
}
