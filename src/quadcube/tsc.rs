//! Tangential Spherical Cube projection.

use crate::{CanonicalProjection, ProjBounds, ProjXY, XYZ};

const TOL: f64 = 1.0e-12;

/// Tangential Spherical Cube projection.
/// 
/// Maps the sphere onto the six faces of a cube.
#[derive(Debug, Clone, Copy)]
pub struct Tsc {
  w0: f64,
  w1: f64,
}

impl Default for Tsc {
  fn default() -> Self {
    Self::new()
  }
}

impl Tsc {
  /// Construct new TSC projection
  #[must_use]
  pub fn new() -> Self {
    let w0 = 45.0;
    let w1 = 1.0 / 45.0;
    
    Self { w0, w1 }
  }
}

impl CanonicalProjection for Tsc {
  const NAME: &'static str = "Tangential Spherical Cube";
  const WCS_NAME: &'static str = "TSC";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(None, None);
    &PROJ_BOUNDS
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let (sin_phi, cos_phi) = xyz.y.atan2(xyz.x).sin_cos();
    let (sin_theta, cos_theta) = xyz.z.asin().sin_cos();
    
    let l = cos_theta * cos_phi;
    let m = cos_theta * sin_phi;
    let n = sin_theta;
    
    // Determine which face
    let mut face = 0;
    let mut zeta = n;
    if l > zeta {
      face = 1;
      zeta = l;
    }
    if m > zeta {
      face = 2;
      zeta = m;
    }
    if -l > zeta {
      face = 3;
      zeta = -l;
    }
    if -m > zeta {
      face = 4;
      zeta = -m;
    }
    if -n > zeta {
      face = 5;
      zeta = -n;
    }
    
    let (xf, yf, x0, y0) = match face {
      1 => (m / zeta, n / zeta, 0.0, 0.0),
      2 => (-l / zeta, n / zeta, 2.0, 0.0),
      3 => (-m / zeta, n / zeta, 4.0, 0.0),
      4 => (l / zeta, n / zeta, 6.0, 0.0),
      5 => (m / zeta, l / zeta, 0.0, -2.0),
      _ => (m / zeta, -l / zeta, 0.0, 2.0), // face 0
    };
    
    // Apply bounds tolerance
    let xf = if xf.abs() > 1.0 {
      if xf.abs() > 1.0 + TOL {
        return None;
      }
      xf.signum()
    } else {
      xf
    };
    
    let yf = if yf.abs() > 1.0 {
      if yf.abs() > 1.0 + TOL {
        return None;
      }
      yf.signum()
    } else {
      yf
    };
    
    let x = self.w0 * (xf + x0);
    let y = self.w0 * (yf + y0);
    
    Some(ProjXY::new(x, y))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let mut xf = pos.x * self.w1;
    let mut yf = pos.y * self.w1;
    
    // Bounds checking
    if xf.abs() <= 1.0 {
      if yf.abs() > 3.0 {
        return None;
      }
    } else if xf.abs() > 7.0 || yf.abs() > 1.0 {
      return None;
    }
    
    // Map negative faces to the other side
    if xf < -1.0 {
      xf += 8.0;
    }
    
    // Determine the face
    let (l, m, n) = if xf > 5.0 {
      // face = 4
      xf -= 6.0;
      let denom = (1.0 + xf * xf + yf * yf).sqrt();
      let m = -1.0 / denom;
      (-m * xf, m, -m * yf)
    } else if xf > 3.0 {
      // face = 3
      xf -= 4.0;
      let denom = (1.0 + xf * xf + yf * yf).sqrt();
      let l = -1.0 / denom;
      (l, l * xf, -l * yf)
    } else if xf > 1.0 {
      // face = 2
      xf -= 2.0;
      let denom = (1.0 + xf * xf + yf * yf).sqrt();
      let m = 1.0 / denom;
      (-m * xf, m, m * yf)
    } else if yf > 1.0 {
      // face = 0
      yf -= 2.0;
      let denom = (1.0 + xf * xf + yf * yf).sqrt();
      let n = 1.0 / denom;
      (-n * yf, n * xf, n)
    } else if yf < -1.0 {
      // face = 5
      yf += 2.0;
      let denom = (1.0 + xf * xf + yf * yf).sqrt();
      let n = -1.0 / denom;
      (-n * yf, -n * xf, n)
    } else {
      // face = 1
      let denom = (1.0 + xf * xf + yf * yf).sqrt();
      let l = 1.0 / denom;
      (l, l * xf, l * yf)
    };
    
    Some(XYZ::new(l, m, n))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_tsc_roundtrip() {
    let proj = Tsc::new();
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
