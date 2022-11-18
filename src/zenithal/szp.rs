//! Slant zenithal perspective projection.

use std::f64::consts::PI;

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

static HALF_PI: f64 = 0.5 * PI;
static EPSILON: f64 = 1e-15;

/// Slant zenithal perspective projection.
pub struct Szp {
  /// Keyword PVi_1a.
  mu: f64, // 0.0
  /// Keyword PVi_2a (converted in radians)
  phi: f64, // 0.0
  /// Keyword PVi_3a (converted in radians)
  theta: f64, // pi/2
  /// Euclidean coordinates of the projection point from the center of the unit
  /// sphere (u) and from the center of the projection plane (p).
  xp: f64,
  yp: f64,
  zp: f64,
  /// Pre-computed variable. Set to `true` if `xp < 0`
  neg_xp: bool,
  /// Pre-computed variable. `(1 - xp)` i.e. distance from point P to the proj plane
  o_m_xp: f64,
  /// Pre-computed variable. `(1 - xp)^2`
  o_m_xp_2: f64,
  /// Pre-computed variable. `|mu|`
  abs_mu: f64,
  /// Pre-computed variable.  `mu^2 - 1`
  mu2_m_1: f64,
}

impl Default for Szp {
  fn default() -> Self {
    Self::new()
  }
}

impl Szp {

  pub fn new() -> Self {
    // Self::from_params(2.0, PI, 2.0 * HALF_PI / 3.0)
    Self::from_params(0.0, 0.0, HALF_PI)
  }

  /// Get the value of the `mu` parameter.
  pub fn mu(&self) -> f64 {
    self.mu
  }

  /// Get the value of the `phi` parameter.
  pub fn phi(&self) -> f64 {
    self.phi
  }

  /// Get the value of the `theta` parameter.
  pub fn theta(&self) -> f64 {
    self.theta
  }
  
  /// # Params
  /// * `mu`: WCS `PVi_1a` parameter.
  /// * `phi`: WCS `PVi_1a` parameter converted in radians.
  /// * `theta`: WCS `PVi_1a` parameter converted in radians.
  /// # Panics
  /// * if `mu` not a finite number
  /// * if `theta` not in `[-PI/2, PI/2]`
  pub fn from_params(mu: f64, phi: f64, theta: f64) -> Self {
    assert!(mu.is_finite());
    assert!((-HALF_PI..=HALF_PI).contains(&theta));
    // Convert from FITS WCS standards to "my" conventions
    let theta_p = HALF_PI - phi;
    let (rho_p, abs_mu) = if mu < 0.0 {
      (HALF_PI - theta, -mu)
    } else {
      (HALF_PI + theta, mu)
    };
    // Compute Cartesian coordinates of the projection point P
    let (sin_rho_p, cos_rho_p) = rho_p.sin_cos();
    let (sin_theta_p, cos_theta_p) = theta_p.sin_cos();
    let r_sr = abs_mu * sin_rho_p;
    let xp = abs_mu * cos_rho_p;
    let yp = r_sr * cos_theta_p;
    let zp = r_sr * sin_theta_p;
    // Compute other useful quantities
    let neg_xp = xp < 0.0;
    let o_m_xp = 1.0 - xp;
    let o_m_xp_2 = o_m_xp.pow2();
    let mu2_m_1 = mu.pow2() - 1.0;
    Self {
      mu, phi, theta,
      xp, yp, zp,
      neg_xp, o_m_xp, o_m_xp_2, abs_mu, mu2_m_1
    }
  }

  pub fn is_in_proj_bounds(&self, xyz: &XYZ) -> bool {
    // We use epsilon instead of 0 to avoid points close to the divergence
    if self.abs_mu <= 1.0 {
      xyz.x - self.xp > EPSILON // point between proj. plane and P
    } else {
      debug_assert!(self.abs_mu > 1.0);
      let sm1 = self.xp * xyz.x + self.yp * xyz.y + self.zp * xyz.z - 1.0;
      if self.neg_xp { // xp < 0 => ap> pi/2 => opposite hemisphere
        sm1 < -EPSILON
      } else { // xp >= 0 => ap <= pi/2 => planewards hemisphere
        sm1 >  EPSILON
      }
    }
  }

}

impl CanonicalProjection for Szp {

  const NAME: &'static str = "Slant zenithal perspective";
  const WCS_NAME: &'static str = "SZP";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    if self.is_in_proj_bounds(xyz) {
      let o_m_x = 1.0 - xyz.x;
      let d = xyz.x - self.xp;
      Some(ProjXY::new(
        (self.o_m_xp * xyz.y - self.yp * o_m_x) / d,
        (self.o_m_xp * xyz.z - self.zp * o_m_x) / d
      ))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let tx = pos.x - self.yp; // X - yp
    let ty = pos.y - self.zp; // Y - zp
    let tx2 = tx.pow2();
    let ty2 = ty.pow2();
    let t = self.xp * self.o_m_xp + self.yp * tx + self.zp * ty;
    // Check if X and Y are in the bounds of the projection area
    if t.pow2() - (self.o_m_xp_2 + tx2 + ty2) * self.mu2_m_1 <= EPSILON {
      None
    } else {
      // Start the deprojection
      let tx = tx / self.o_m_xp;
      let ty = ty / self.o_m_xp;
      let txp = tx * self.xp - self.yp;
      let typ = ty * self.xp - self.zp;
      let a = tx2 + ty2 + 1.0;
      let b = -(tx * txp + ty * typ);
      let c = txp.pow2() + typ.pow2() - 1.0;
      let x = (-b + (b.pow2() - a * c).sqrt()) / a;
      // Double check (should be useless!)
      if self.abs_mu > 1.0 && (!x.is_finite() || !(-1.0..=1.0).contains(&x)) {
        None
      } else {
        let x_m_xp = x - self.xp;
        Some(XYZ::new_renorming_if_necessary(x, tx * x_m_xp + self.yp, ty * x_m_xp + self.zp)) // new_renorming_if_necessary
      }
    }
  }
}
