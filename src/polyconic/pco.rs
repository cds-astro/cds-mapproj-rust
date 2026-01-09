//! Polyconic projection.

use crate::{CanonicalProjection, ProjBounds, ProjXY, XYZ};
use std::f64::consts::PI;

const D2R: f64 = PI / 180.0;
const R2D: f64 = 180.0 / PI;
const TOL: f64 = 1.0e-12;

/// Polyconic projection.
///
/// The American polyconic projection.
#[derive(Debug, Clone, Copy)]
pub struct Pco {
  r0: f64,
  w0: f64,
  w1: f64,
  w2: f64,
  w3: f64,
}

impl Default for Pco {
  fn default() -> Self {
    Self::new()
  }
}

impl Pco {
  /// Construct new polyconic projection
  #[must_use]
  pub fn new() -> Self {
    let r0 = R2D;
    let w0 = 1.0;
    let w1 = 1.0;
    let w2 = 360.0 / PI;
    let w3 = D2R / w2;

    Self { r0, w0, w1, w2, w3 }
  }
}

impl CanonicalProjection for Pco {
  const NAME: &'static str = "Polyconic";
  const WCS_NAME: &'static str = "PCO";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(None, None);
    &PROJ_BOUNDS
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // Convert from unit vector to spherical coordinates
    let lon = xyz.y.atan2(xyz.x);
    let lat = xyz.z.asin();

    let phi = lon;
    let theta = lat;

    if theta == 0.0 {
      let x = self.w0 * phi;
      let y = 0.0;
      return Some(ProjXY::new(x, y));
    }

    if theta.abs() < 1.0e-4 {
      // Small angle approximation to avoid cot(theta) blowing up
      let x = self.w0 * phi * theta.cos();
      let y = (self.w0 + self.w3 * phi * phi) * theta;
      return Some(ProjXY::new(x, y));
    }

    let (sin_theta, cos_theta) = theta.sin_cos();
    let (sin_psi, cos_psi) = (phi * sin_theta).sin_cos();
    let cot_theta = cos_theta / sin_theta;

    let x = self.r0 * cot_theta * sin_psi;
    let y = self.r0 * (cot_theta * (1.0 - cos_psi) + theta);

    Some(ProjXY::new(x, y))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let xj = pos.x;
    let yj = pos.y;

    let w = (yj * self.w1).abs();

    let (phi, theta) = if w < TOL {
      // Near the equator
      (xj * self.w1, 0.0)
    } else if (w - 90.0).abs() < TOL {
      // Near the pole
      (0.0, yj.signum() * 90.0 * D2R)
    } else {
      // Iterative solution using weighted division of the interval
      let (the, ymthe, tan_the) = if w < 1.0e-4 {
        // To avoid cot(theta) blowing up near theta == 0
        let the = yj / (self.w0 + self.w3 * xj * xj);
        let ymthe = yj - self.w0 * the;
        let tan_the = the.tan();
        (the, ymthe, tan_the)
      } else {
        // Iterative solution
        let mut the_pos = yj / self.w0;
        let mut the_neg = 0.0;

        let xx = xj * xj;
        let mut f_pos = xx;
        let mut f_neg = -xx;

        let mut the = the_pos;
        let mut ymthe = 0.0;
        let mut tan_the = 0.0;

        for _ in 0..64 {
          // Weighted division of the interval
          let lambda = (f_pos / (f_pos - f_neg)).clamp(0.1, 0.9);
          the = the_pos - lambda * (the_pos - the_neg);

          // Compute the residue
          ymthe = yj - self.w0 * the;
          tan_the = the.tan();
          let f = xx + ymthe * (ymthe - self.w2 / tan_the);

          // Check for convergence
          if f.abs() < TOL {
            break;
          }
          if (the_pos - the_neg).abs() < TOL {
            break;
          }

          // Redefine the interval
          if f > 0.0 {
            the_pos = the;
            f_pos = f;
          } else {
            the_neg = the;
            f_neg = f;
          }
        }
        (the, ymthe, tan_the)
      };

      let x1 = self.r0 - ymthe * tan_the;
      let y1 = xj * tan_the;

      let phi = if x1 == 0.0 && y1 == 0.0 {
        0.0
      } else {
        y1.atan2(x1) / the.sin()
      };

      (phi, the)
    };

    // Convert spherical to unit vector
    let (sin_theta, cos_theta) = theta.sin_cos();
    let (sin_phi, cos_phi) = phi.sin_cos();

    let x = cos_theta * cos_phi;
    let y = cos_theta * sin_phi;
    let z = sin_theta;

    Some(XYZ::new(x, y, z))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_pco_roundtrip() {
    let proj = Pco::new();
    let xyz = XYZ::new(1.0, 0.0, 0.0);

    if let Some(proj_xy) = proj.proj(&xyz) {
      if let Some(xyz_back) = proj.unproj(&proj_xy) {
        assert!((xyz.x - xyz_back.x).abs() < 1e-8);
        assert!((xyz.y - xyz_back.y).abs() < 1e-8);
        assert!((xyz.z - xyz_back.z).abs() < 1e-8);
      }
    }
  }
}
