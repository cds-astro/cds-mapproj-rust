//! Bonne's projection.

use crate::{CanonicalProjection, ProjBounds, ProjXY, XYZ};
use std::f64::consts::PI;

const D2R: f64 = PI / 180.0;
const R2D: f64 = 180.0 / PI;

/// Bonne's projection (equal-area polyconic projection).
///
/// This projection requires a parameter `theta1` (PV i_1), the standard parallel.
/// If theta1 = 0, this reduces to the Sanson-Flamsteed projection (SFL).
#[derive(Debug, Clone, Copy)]
pub struct Bon {
  theta1: f64, // Standard parallel in radians
  r0: f64,
  w1: f64,
  w2: f64,
}

impl Default for Bon {
  fn default() -> Self {
    Self::new()
  }
}

impl Bon {
  /// Construct with default theta1 = 60 degrees
  #[must_use]
  pub fn new() -> Self {
    Self::from_params(60.0 * D2R)
  }

  /// Construct from provided `theta1` parameter (in radians).
  /// theta1 is the standard parallel.
  #[must_use]
  pub fn from_params(theta1: f64) -> Self {
    let r0 = R2D;
    let w1 = 1.0;
    let w2 = r0 * theta1.cos() / theta1.sin() + theta1;

    Self { theta1, r0, w1, w2 }
  }

  /// Construct from provided `theta1` parameter in degrees.
  #[must_use]
  pub fn from_params_deg(theta1_deg: f64) -> Self {
    Self::from_params(theta1_deg * D2R)
  }
}

impl CanonicalProjection for Bon {
  const NAME: &'static str = "Bonne's";
  const WCS_NAME: &'static str = "BON";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(None, None);
    &PROJ_BOUNDS
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    // Convert from unit vector to spherical coordinates
    let lon = xyz.y.atan2(xyz.x);
    let lat = xyz.z.asin();

    // Bonne's projection
    let theta = lat;
    let phi = lon;

    let r = self.w2 - self.w1 * theta;
    let cos_theta = theta.cos();

    if cos_theta == 0.0 {
      // At the pole
      return Some(ProjXY::new(0.0, self.w2));
    }

    let alpha = cos_theta / r;
    let s = alpha * phi;

    let (sin_s, cos_s) = s.sin_cos();
    let x = r * sin_s;
    let y = -r * cos_s - (-self.w2);

    Some(ProjXY::new(x, y))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let xj = pos.x;
    let yj = pos.y;

    let dy = self.w2 - yj;
    let dy2 = dy * dy;

    let r = (xj * xj + dy2).sqrt();
    let r = if self.theta1 < 0.0 { -r } else { r };

    let alpha = if r == 0.0 {
      0.0
    } else {
      (xj / r).atan2(dy / r)
    };

    let theta = (self.w2 - r) / self.w1;
    let cos_theta = theta.cos();

    let phi = if cos_theta == 0.0 {
      0.0
    } else {
      alpha * (r / self.r0) / cos_theta
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
  fn test_bon_roundtrip() {
    let proj = Bon::new();
    let xyz = XYZ::new(1.0, 0.0, 0.0);

    if let Some(proj_xy) = proj.proj(&xyz) {
      if let Some(xyz_back) = proj.unproj(&proj_xy) {
        assert!((xyz.x - xyz_back.x).abs() < 1e-10);
        assert!((xyz.y - xyz_back.y).abs() < 1e-10);
        assert!((xyz.z - xyz_back.z).abs() < 1e-10);
      }
    }
  }
}
