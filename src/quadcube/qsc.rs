//! Quadrilateralized spherical cube projection.

use crate::{CanonicalProjection, ProjBounds, ProjXY, XYZ};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

const TOL: f64 = 1.0e-12;

/// Quadrilateralized Spherical Cube projection.
///
/// This is an equal-area projection that maps the sphere onto
/// the six faces of a cube.
#[derive(Debug, Clone, Copy)]
pub struct Qsc {
  w0: f64,
  w1: f64,
}

impl Default for Qsc {
  fn default() -> Self {
    Self::new()
  }
}

impl Qsc {
  /// Construct new QSC projection
  #[must_use]
  pub fn new() -> Self {
    let w0 = 45.0;
    let w1 = 1.0 / 45.0;

    Self { w0, w1 }
  }
}

impl CanonicalProjection for Qsc {
  const NAME: &'static str = "Quadrilateralized Spherical Cube";
  const WCS_NAME: &'static str = "QSC";

  fn bounds(&self) -> &ProjBounds {
    const PROJ_BOUNDS: ProjBounds = ProjBounds::new(None, None);
    &PROJ_BOUNDS
  }

  fn proj(&self, xyz: &XYZ) -> Option<ProjXY> {
    let lat = xyz.z.asin();

    // Handle poles specially
    if lat.abs() == PI / 2.0 {
      let x = 0.0;
      let y = lat.signum() * 2.0 * self.w0;
      return Some(ProjXY::new(x, y));
    }

    let (sin_phi, cos_phi) = xyz.y.atan2(xyz.x).sin_cos();
    let (sin_theta, cos_theta) = lat.sin_cos();

    let l = cos_theta * cos_phi;
    let m = cos_theta * sin_phi;
    let n = sin_theta;

    // Determine which face has maximum component
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

    let zeco = 1.0 - zeta;

    let (xi, eta, x0, y0, zeco_adj) = match face {
      1 => {
        let xi = m;
        let eta = n;
        let zeco_adj = if zeco < 1.0e-8 {
          // Small angle formula
          let t = lat;
          let p = xyz.y.atan2(xyz.x);
          (p.mul_add(p, t * t)) / 2.0
        } else {
          zeco
        };
        (xi, eta, 0.0, 0.0, zeco_adj)
      }
      2 => {
        let xi = -l;
        let eta = n;
        let zeco_adj = if zeco < 1.0e-8 {
          let t = lat;
          let p = xyz.y.atan2(xyz.x) - PI / 2.0;
          (p.mul_add(p, t * t)) / 2.0
        } else {
          zeco
        };
        (xi, eta, 2.0, 0.0, zeco_adj)
      }
      3 => {
        let xi = -m;
        let eta = n;
        let zeco_adj = if zeco < 1.0e-8 {
          let t = lat;
          let mut p = xyz.y.atan2(xyz.x);
          p -= p.signum() * PI;
          (p.mul_add(p, t * t)) / 2.0
        } else {
          zeco
        };
        (xi, eta, 4.0, 0.0, zeco_adj)
      }
      4 => {
        let xi = l;
        let eta = n;
        let zeco_adj = if zeco < 1.0e-8 {
          let t = lat;
          let p = xyz.y.atan2(xyz.x) + PI / 2.0;
          (p.mul_add(p, t * t)) / 2.0
        } else {
          zeco
        };
        (xi, eta, 6.0, 0.0, zeco_adj)
      }
      5 => {
        let xi = m;
        let eta = l;
        let zeco_adj = if zeco < 1.0e-8 {
          let t = lat + PI / 2.0;
          t * t / 2.0
        } else {
          zeco
        };
        (xi, eta, 0.0, -2.0, zeco_adj)
      }
      _ => {
        // face 0
        let xi = m;
        let eta = -l;
        let zeco_adj = if zeco < 1.0e-8 {
          let t = PI / 2.0 - lat;
          t * t / 2.0
        } else {
          zeco
        };
        (xi, eta, 0.0, 2.0, zeco_adj)
      }
    };

    let (xf, yf) = if xi == 0.0 && eta == 0.0 {
      (0.0, 0.0)
    } else if -xi > eta.abs() {
      let omega = eta / xi;
      let tau = 1.0 + omega * omega;
      let xf = -(zeco_adj / (1.0 - 1.0 / (1.0 + tau).sqrt())).sqrt();
      let yf = (xf / 15.0)
        * (omega.atan().to_degrees() - (omega / (tau + tau).sqrt()).asin().to_degrees());
      (xf, yf)
    } else if xi > eta.abs() {
      let omega = eta / xi;
      let tau = 1.0 + omega * omega;
      let xf = (zeco_adj / (1.0 - 1.0 / (1.0 + tau).sqrt())).sqrt();
      let yf = (xf / 15.0)
        * (omega.atan().to_degrees() - (omega / (tau + tau).sqrt()).asin().to_degrees());
      (xf, yf)
    } else if -eta >= xi.abs() {
      let omega = xi / eta;
      let tau = 1.0 + omega * omega;
      let yf = -(zeco_adj / (1.0 - 1.0 / (1.0 + tau).sqrt())).sqrt();
      let xf = (yf / 15.0)
        * (omega.atan().to_degrees() - (omega / (tau + tau).sqrt()).asin().to_degrees());
      (xf, yf)
    } else {
      // eta >= xi.abs()
      let omega = xi / eta;
      let tau = 1.0 + omega * omega;
      let yf = (zeco_adj / (1.0 - 1.0 / (1.0 + tau).sqrt())).sqrt();
      let xf = (yf / 15.0)
        * (omega.atan().to_degrees() - (omega / (tau + tau).sqrt()).asin().to_degrees());
      (xf, yf)
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
    let face = if xf > 5.0 {
      xf -= 6.0;
      4
    } else if xf > 3.0 {
      xf -= 4.0;
      3
    } else if xf > 1.0 {
      xf -= 2.0;
      2
    } else if yf > 1.0 {
      yf -= 2.0;
      0
    } else if yf < -1.0 {
      yf += 2.0;
      5
    } else {
      1
    };

    let direct = xf.abs() > yf.abs();

    let (omega, tau, zeta) = if direct && xf != 0.0 {
      let w = 15.0 * yf / xf;
      let (sin_w, cos_w) = (w.to_radians().sin(), w.to_radians().cos());
      let omega = sin_w / (cos_w - FRAC_1_SQRT_2);
      let tau = 1.0 + omega * omega;
      let zeco = xf * xf * (1.0 - 1.0 / (1.0 + tau).sqrt());
      let zeta = 1.0 - zeco;
      (omega, tau, zeta)
    } else if !direct && yf != 0.0 {
      let w = 15.0 * xf / yf;
      let (sin_w, cos_w) = (w.to_radians().sin(), w.to_radians().cos());
      let omega = sin_w / (cos_w - FRAC_1_SQRT_2);
      let tau = 1.0 + omega * omega;
      let zeco = yf * yf * (1.0 - 1.0 / (1.0 + tau).sqrt());
      let zeta = 1.0 - zeco;
      (omega, tau, zeta)
    } else {
      (0.0, 1.0, 1.0)
    };

    if zeta < -1.0 - TOL {
      return None;
    }

    let zeta = zeta.max(-1.0);
    let zeco = 1.0 - zeta;

    let w = if zeta >= -1.0 {
      (zeco * (2.0 - zeco) / tau).sqrt()
    } else {
      0.0
    };

    let (l, m, n) = match face {
      1 => {
        let l = zeta;
        if direct {
          let m = if xf < 0.0 { -w } else { w };
          let n = m * omega;
          (l, m, n)
        } else {
          let n = if yf < 0.0 { -w } else { w };
          let m = n * omega;
          (l, m, n)
        }
      }
      2 => {
        let m = zeta;
        if direct {
          let l = if xf > 0.0 { -w } else { w };
          let n = -l * omega;
          (l, m, n)
        } else {
          let n = if yf < 0.0 { -w } else { w };
          let l = -n * omega;
          (l, m, n)
        }
      }
      3 => {
        let l = -zeta;
        if direct {
          let m = if xf > 0.0 { -w } else { w };
          let n = -m * omega;
          (l, m, n)
        } else {
          let n = if yf < 0.0 { -w } else { w };
          let m = -n * omega;
          (l, m, n)
        }
      }
      4 => {
        let m = -zeta;
        if direct {
          let l = if xf < 0.0 { -w } else { w };
          let n = l * omega;
          (l, m, n)
        } else {
          let n = if yf < 0.0 { -w } else { w };
          let l = n * omega;
          (l, m, n)
        }
      }
      5 => {
        let n = -zeta;
        if direct {
          let m = if xf < 0.0 { -w } else { w };
          let l = m * omega;
          (l, m, n)
        } else {
          let l = if yf < 0.0 { -w } else { w };
          let m = l * omega;
          (l, m, n)
        }
      }
      _ => {
        // face 0
        let n = zeta;
        if direct {
          let m = if xf < 0.0 { -w } else { w };
          let l = -m * omega;
          (l, m, n)
        } else {
          let l = if yf > 0.0 { -w } else { w };
          let m = -l * omega;
          (l, m, n)
        }
      }
    };

    Some(XYZ::new(l, m, n))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_qsc_roundtrip() {
    let proj = Qsc::new();
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
