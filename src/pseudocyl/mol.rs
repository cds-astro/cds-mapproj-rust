//! Mollweide projection.

use std::f64::consts::{PI, SQRT_2};

use crate::{CanonicalProjection, CustomFloat, ProjXY, XYZ};

static HALF_PI: f64 = 0.5 * PI;

/// Mollweide projection.
pub struct Mol {
  /// Max number of iteration for the Newton-Raphson iterative method.
  n_iter: u8,
  /// Precision for the Newton-Raphson iterative method.
  eps: f64,
}

impl Default for Mol {
  fn default() -> Self {
    Self::new()
  }
}

impl Mol {
  
  pub fn new() -> Self {
    Self {
      n_iter: 100,
      eps: (1.0 / (2.0_f64 * 60.0 * 60.0 * 1000.0 * 2000.0)).to_radians() / PI, // <=> max error (at b = PI) of 1 microarcsec
    }
  }
  
  /// Set the max number of iteration for the Newton method.
  pub fn set_n_iter(&mut self, n_iter: u8) {
    self.n_iter = n_iter;
  }
  
  /// Set the precision (in radians) to stop the Newton algo.
  pub fn set_epsilon(&mut self, eps: f64) {
    self.eps = eps;
  }


  fn newton_solve(&self, z: f64) -> f64 {
    let cte = PI * z;
    // Initial guess so that for z ~= 1, gamma ~= PI/2.
    // Smooth function for small |z|, so no big deal having a bad init value.
    let mut x = 2.0 * z.asin();
    let mut f = x + x.sin() - cte; 
    let mut i = 0_u8;
    while f.abs() > self.eps && i < self.n_iter {
      x -= f / (1.0 + x.cos()); // VERIFIER COEFF 2!!
      f = x + x.sin() - cte;
      i += 1;
    }
    0.5 * x
  }
  
}

impl CanonicalProjection for Mol {

  const NAME: &'static str = "Mollweide";
  const WCS_NAME: &'static str = "MOL";

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // find g iteratively using Newton-Raphson
    // Solve f(g) = 0, with f(g) = 2g + sin(2g) - pi*z
    //                 and f'(g) = 2 [1 + cos(2g)]
    //                           = 4 cos^2(g)
    // x=2g  f(x) = x + sin(x) - pi*z
    //      f'(x) = 1 + cos(x)
    // => f(x) / f'(x) = (x + sin(x) - pi*z) / (1 + cos(x))
    // So we have a problem if x near from +-pi, since (cos(+-pi) = -1)
    // Determine experimentaly possible range of z (problematic if
    // z near from +-1
    let g = self.newton_solve(xyz.z);
    if (-HALF_PI..=HALF_PI).contains(&g) {
      let (sing, cosg) = g.sin_cos();
      Some(ProjXY::new((SQRT_2 * xyz.y.atan2(xyz.x) * cosg).twice() / PI, SQRT_2 * sing))
    } else {
      None
    }
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let sqrt_2_m_y2 = 2.0 - pos.y.pow2();
    if sqrt_2_m_y2 <= 0.0 {
      let z = if pos.y > 0.0 { 1.0 } else { -1.0 };
      Some(XYZ::new(0.0, 0.0, z))
    } else {
      let sqrt_2_m_y2 = sqrt_2_m_y2.sqrt();
      let z = ((pos.y / SQRT_2).asin().twice() + pos.y * sqrt_2_m_y2) / PI;
      if (-1.0..=1.0).contains(&z) {
        let (sinl, cosl) = ((pos.x * HALF_PI) / sqrt_2_m_y2).sin_cos();
        let r = (1.0 - z.pow2()).sqrt();
        Some(XYZ::new(r * cosl, r * sinl, z))
      } else {
        // Should not happen!
        None
      }
    }
  }
}
