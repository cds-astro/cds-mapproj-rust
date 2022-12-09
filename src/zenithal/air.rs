//! Airy projection.

use std::f64::consts::PI;

use crate::{CanonicalProjection, CustomFloat, ProjBounds, ProjXY, XYZ};

static HALF_PI: f64 = 0.5 * PI;

/// Airy projection.
pub struct Air {
  /// Angular distance (in radians) from the proj center at which error is minimized
  rho_b: f64,
  /// Precomputed constant deriving from `rho_b`
  cte_b: f64,
  /// Max number of iteration for the Newton-Raphson iterative method.
  n_iter: u8,
  /// Precision for the Newton-Raphson iterative method
  eps: f64,
}

impl Default for Air {
  fn default() -> Self {
    Self::new()
  }
}

impl Air {
  
  pub fn new() -> Self {
    Self::from_param(HALF_PI)
  }
  
  /// # Params
  /// * `rho_b`: WCS parameter `PVi_2a` (but converted in radians);
  ///            angular distance from the proj center at which error is minimized
  /// # Panics
  /// * if `rho_b` no in `]0, pi[`
  pub fn from_param(rho_b: f64) -> Self {
    let xb = rho_b.cos();
    assert!(xb.abs() < 1.0, "In AIR, angle must be in ]0, pi[");
    let xb_p_1 = xb + 1.0;
    let cte_b = xb_p_1 * (0.5 * xb_p_1).ln() / (1.0 - xb);
    Self {
      rho_b,
      cte_b,
      n_iter: 100,
      eps: (1.0_f64 / (2.0_f64 * 60.0 * 60.0 * 1000.0 * 2000.0)).to_radians() / PI // <=> max error of 1 mas (at b = PI)
    }
  }
  
  /// Get the value of the `rho_b` parameter.
  pub fn rho_b(&self) -> f64 {
    self.rho_b
  }
  
  /// Set the max number of iteration for the Newton-Raphson iterative method.
  pub fn set_n_iter(&mut self, n_iter: u8) {
    self.n_iter = n_iter;
  }

  /// Set precision for the Newton-Raphson iterative method
  pub fn set_eps(&mut self, eps: f64) {
    self.eps = eps;
  }

  fn newton_solve(&self, big_r: f64) -> f64 {
    let mut x = 0.0_f64;
    let mut x_p_1 = 1.0 + x;
    let mut x_m_1 = 1.0 - x;
    let mut ln = (0.5 * x_p_1).ln();
    let mut sqrt_1_p_x = x_p_1.sqrt();
    let mut sqrt_1_m_x = x_m_1.sqrt();
    let mut sqrt_1p_1m = -sqrt_1_p_x / sqrt_1_m_x;
    let mut i = 0_u8;
    let mut f = sqrt_1p_1m * ln + (self.cte_b / sqrt_1p_1m) - big_r;
    while i < self.n_iter && f.abs() > self.eps {
      x -= -sqrt_1_p_x * sqrt_1_m_x * f / (1.0 + ln / x_m_1 - self.cte_b / x_p_1);
      if x <= -1.0 {
        x = -1.0 + 1e-14;
      } else if x >= 1.0 {
        x =  1.0 - 1e-14;
      }
      x_p_1 = 1.0 + x;
      x_m_1 = 1.0 - x;
      ln = (0.5 * x_p_1).ln();
      sqrt_1_p_x = x_p_1.sqrt();
      sqrt_1_m_x = x_m_1.sqrt();
      sqrt_1p_1m = -sqrt_1_p_x / sqrt_1_m_x;
      i += 1;
      f = sqrt_1p_1m * ln + (self.cte_b / sqrt_1p_1m) - big_r
    }
    x
  }
}



impl CanonicalProjection for Air {

  const NAME: &'static str = "Airy";
  const WCS_NAME: &'static str = "AIR";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(
      None,
      None
    );
    &PROJ_BOUNDS
  }
  
  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let x_p_1 = xyz.x + 1.0;
    if xyz.x == 1.0 { // rho = 0
      Some(ProjXY::new(0.0, 0.0))
    } else if x_p_1 == 0.0 { // rho = PI
      debug_assert!(xyz.x == -1.0);
      None
    } else {
      // Compute R/sin(rho): r over sin rho
      debug_assert!(1.0 - xyz.x > 0.0, "1 - x must be > 0");
      let r_o_sr = (0.5 * x_p_1).ln() / (xyz.x - 1.0) - (self.cte_b / x_p_1);
      debug_assert!(r_o_sr >= 0.0, "r: {} < 0", r_o_sr);
      Some(ProjXY::new(xyz.y * r_o_sr, xyz.z * r_o_sr))
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let big_r = (pos.x.pow2() + pos.y.pow2()).sqrt();
    if big_r == 0.0 {
      Some(XYZ::new(1.0, 0.0, 0.0))
    } else {
      let x = self.newton_solve(big_r);
      let w = (1.0 - x.pow2()).sqrt() / big_r;
      Some(XYZ::new(x, pos.x * w, pos.y * w))
    }
  }
}

   

