//! COBE quadrilateralized spherical cube projection.

use crate::{CanonicalProjection, ProjBounds, ProjXY, XYZ};

const TOL: f64 = 1.0e-7;

// Polynomial coefficients for CSC deprojection

const P00: f32 = -0.27292696;
const P10: f32 = -0.07629969;
const P20: f32 = -0.22797056;
const P30: f32 = 0.54852384;
const P40: f32 = -0.62930065;
const P50: f32 = 0.25795794;
const P60: f32 = 0.02584375;
const P01: f32 = -0.02819452;
const P11: f32 = -0.01471565;
const P21: f32 = 0.48051512;
const P31: f32 = -1.7411445;
const P41: f32 = 1.7154751;
const P51: f32 = -0.53022337;
const P02: f32 = 0.2705816;
const P12: f32 = -0.5680094;
const P22: f32 = 0.30803317;
const P32: f32 = 0.989381;
const P42: f32 = -0.8318047;
const P03: f32 = -0.6044156;
const P13: f32 = 1.5088009;
const P23: f32 = -0.93678576;
const P33: f32 = 0.08693841;
const P04: f32 = 0.9341208;
const P14: f32 = -1.4160192;
const P24: f32 = 0.33887446;
const P05: f32 = -0.63915306;
const P15: f32 = 0.5203224;
const P06: f32 = 0.14381585;

// Polynomial coefficients for CSC projection

const GSTAR: f32 = 1.3748485;
const MM: f32 = 0.004869492;
const GAMMA: f32 = -0.13161671;
const OMEGA1: f32 = -0.15959623;
const D0: f32 = 0.07591962;
const D1: f32 = -0.02177625;
const C00: f32 = 0.14118963;
const C10: f32 = 0.08097013;
const C01: f32 = -0.28152853;
const C11: f32 = 0.15384113;
const C20: f32 = -0.1782512;
const C02: f32 = 0.10695947;

/// COBE quadrilateralized spherical cube projection.
///
/// Used by the COBE satellite, this projection maps the sphere onto
/// the six faces of a cube with a specific polynomial transformation.
#[derive(Debug, Clone, Copy)]
pub struct Csc {
  w0: f64,
  w1: f64,
}

impl Default for Csc {
  fn default() -> Self {
    Self::new()
  }
}

impl Csc {
  /// Construct new CSC projection
  #[must_use]
  pub fn new() -> Self {
    let w0 = 45.0;
    let w1 = 1.0 / 45.0;

    Self { w0, w1 }
  }
}

impl CanonicalProjection for Csc {
  const NAME: &'static str = "COBE Quadrilateralized Spherical Cube";
  const WCS_NAME: &'static str = "CSC";

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

    let (xi, eta, x0, y0) = match face {
      1 => (m, n, 0.0, 0.0),
      2 => (-l, n, 2.0, 0.0),
      3 => (-m, n, 4.0, 0.0),
      4 => (l, n, 6.0, 0.0),
      5 => (m, l, 0.0, -2.0),
      _ => (m, -l, 0.0, 2.0), // face 0
    };

    let chi = (xi / zeta) as f32;
    let psi = (eta / zeta) as f32;

    let chi2 = chi * chi;
    let psi2 = psi * psi;
    let chi2co = 1.0 - chi2;
    let psi2co = 1.0 - psi2;

    // Avoid floating underflows
    let chipsi = (chi * psi).abs();
    let chi4 = if chi2 > 1.0e-16 { chi2 * chi2 } else { 0.0 };
    let psi4 = if psi2 > 1.0e-16 { psi2 * psi2 } else { 0.0 };
    let chi2psi2 = if chipsi > 1.0e-16 { chi2 * psi2 } else { 0.0 };

    let xf = chi
      * (chi2
        + chi2co
          * (GSTAR
            + psi2
              * (GAMMA * chi2co
                + MM * chi2
                + psi2co
                  * (C00 + C10 * chi2 + C01 * psi2 + C11 * chi2psi2 + C20 * chi4 + C02 * psi4))
            + chi2 * (OMEGA1 - chi2co * (D0 + D1 * chi2))));
    let yf = psi
      * (psi2
        + psi2co
          * (GSTAR
            + chi2
              * (GAMMA * psi2co
                + MM * psi2
                + chi2co
                  * (C00 + C10 * psi2 + C01 * chi2 + C11 * chi2psi2 + C20 * psi4 + C02 * chi4))
            + psi2 * (OMEGA1 - psi2co * (D0 + D1 * psi2))));

    // Apply bounds tolerance
    let xf = if xf.abs() > 1.0 {
      if f64::from(xf.abs()) > 1.0 + TOL {
        return None;
      }
      xf.signum()
    } else {
      xf
    };

    let yf = if yf.abs() > 1.0 {
      if f64::from(yf.abs()) > 1.0 + TOL {
        return None;
      }
      yf.signum()
    } else {
      yf
    };

    let x = self.w0 * (f64::from(xf) + x0);
    let y = self.w0 * (f64::from(yf) + y0);

    Some(ProjXY::new(x, y))
  }

  fn unproj(&self, pos: &ProjXY) -> Option<XYZ> {
    let mut xf = (pos.x * self.w1) as f32;
    let mut yf = (pos.y * self.w1) as f32;

    // Bounds checking
    if f64::from(xf.abs()) <= 1.0 {
      if f64::from(yf.abs()) > 3.0 {
        return None;
      }
    } else if f64::from(xf.abs()) > 7.0 || f64::from(yf.abs()) > 1.0 {
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

    let xx = xf * xf;
    let yy = yf * yf;

    // Compute chi using polynomial in x
    let z0 = P00 + xx * (P10 + xx * (P20 + xx * (P30 + xx * (P40 + xx * (P50 + xx * P60)))));
    let z1 = P01 + xx * (P11 + xx * (P21 + xx * (P31 + xx * (P41 + xx * P51))));
    let z2 = P02 + xx * (P12 + xx * (P22 + xx * (P32 + xx * P42)));
    let z3 = P03 + xx * (P13 + xx * (P23 + xx * P33));
    let z4 = P04 + xx * (P14 + xx * P24);
    let z5 = P05 + xx * P15;
    let z6 = P06;

    let chi = z0 + yy * (z1 + yy * (z2 + yy * (z3 + yy * (z4 + yy * (z5 + yy * z6)))));
    let chi = xf + xf * (1.0 - xx) * chi;

    // Compute psi using polynomial in y
    let z0 = P00 + yy * (P10 + yy * (P20 + yy * (P30 + yy * (P40 + yy * (P50 + yy * P60)))));
    let z1 = P01 + yy * (P11 + yy * (P21 + yy * (P31 + yy * (P41 + yy * P51))));
    let z2 = P02 + yy * (P12 + yy * (P22 + yy * (P32 + yy * P42)));
    let z3 = P03 + yy * (P13 + yy * (P23 + yy * P33));
    let z4 = P04 + yy * (P14 + yy * P24);
    let z5 = P05 + yy * P15;
    let z6 = P06;

    let psi = z0 + xx * (z1 + xx * (z2 + xx * (z3 + xx * (z4 + xx * (z5 + xx * z6)))));
    let psi = yf + yf * (1.0 - yy) * psi;

    // Convert to unit vector
    let t = 1.0 / (f64::from(chi.mul_add(chi, psi * psi)) + 1.0).sqrt();
    let (l, m, n) = match face {
      1 => {
        let l = t;
        (l, f64::from(chi) * l, f64::from(psi) * l)
      }
      2 => {
        let m = t;
        (-f64::from(chi) * m, m, f64::from(psi) * m)
      }
      3 => {
        let l = -t;
        (l, f64::from(chi) * l, -f64::from(psi) * l)
      }
      4 => {
        let m = -t;
        (-f64::from(chi) * m, m, -f64::from(psi) * m)
      }
      5 => {
        let n = -t;
        (-f64::from(psi) * n, -f64::from(chi) * n, n)
      }
      _ => {
        // face 0
        let n = t;
        (-f64::from(psi) * n, f64::from(chi) * n, n)
      }
    };

    Some(XYZ::new(l, m, n))
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_csc_roundtrip() {
    let proj = Csc::new();
    let xyz = XYZ::new(1.0, 0.0, 0.0);

    if let Some(proj_xy) = proj.proj(&xyz) {
      if let Some(xyz_back) = proj.unproj(&proj_xy) {
        assert!((xyz.x - xyz_back.x).abs() < 1e-6);
        assert!((xyz.y - xyz_back.y).abs() < 1e-6);
        assert!((xyz.z - xyz_back.z).abs() < 1e-6);
      }
    }
  }
}
