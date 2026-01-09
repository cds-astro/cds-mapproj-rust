//! TPV (Tangent Plane + Polynomial) distortion convention.
//!
//! TPV extends the TAN projection with polynomial distortion correction,
//! widely used in astronomy to model telescope optical distortions.
//!
//! ## Formulas
//!
//! ```text
//! xi'  = PV1_0 + PV1_1*xi + PV1_2*eta + PV1_3*r + ... + PV1_39*r^7
//! eta' = PV2_0 + PV2_1*eta + PV2_2*xi + PV2_3*r + ... + PV2_39*r^7
//! ```
//!
//! Where r = sqrt(xi^2 + eta^2). Note: Axis 2 swaps xi/eta in linear terms.
//! **IMPORTANT**: xi, eta, xi', eta' are all in DEGREES (intermediate world coords).
//! Per TPV specification: "the units of xi, eta, f_xi, and f_eta are also degrees."
//!
//! ## Example
//!
//! ```rust
//! use mapproj::tpv::{Tpv, TpvCoeff, TpvPV};
//!
//! // Radial distortion: xi' = xi + 0.01*r, eta' = eta + 0.01*r
//! let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.01]);
//! let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.01]);
//! let tpv = Tpv::new(TpvPV::new(pv1, pv2));
//!
//! // Forward and inverse (values in degrees)
//! let (xi_p, eta_p) = tpv.project(10.0, 5.0);
//! if let Some(result) = tpv.inverse(xi_p, eta_p) {
//!     println!("Round-trip: ({}, {})", result.x(), result.y());
//! }
//! ```
//!
//! ## Reference
//!
//! <https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html>

use crate::{CustomFloat, ImgXY};

/// TPV polynomial coefficients for a single axis.
///
/// Stores 40 coefficients (PV0-PV39) including Cartesian (xi^p*eta^q) and
/// radial (r, r^3, r^5, r^7) terms. Defaults: PV1_1=1.0, PV2_1=1.0, others=0.
///
/// # Examples
///
/// ```rust
/// use mapproj::tpv::TpvCoeff;
///
/// let identity = TpvCoeff::new_axis1(&[]);           // xi' = xi
/// let scale = TpvCoeff::new_axis1(&[0.0, 2.0]);      // xi' = 2*xi
/// let radial = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.01]); // xi' = xi + 0.01*r
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TpvCoeff {
  /// Polynomial coefficients (PV0 through PV39)
  /// Index i corresponds to PVn_i coefficient
  coefficients: [f64; 40],

  /// Whether this is for axis 1 or axis 2 (affects evaluation formula)
  axis: TpvAxis,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum TpvAxis {
  Axis1,
  Axis2,
}

impl TpvCoeff {
  /// Create TPV coefficients for axis 1 (xi -> xi').
  ///
  /// Missing coefficients default to 0, except PV1_1 defaults to 1.0.
  #[must_use]
  pub fn new_axis1(coeffs: &[f64]) -> Self {
    Self::new_axis(coeffs, TpvAxis::Axis1)
  }

  /// Create TPV coefficients for axis 2 (eta -> eta').
  ///
  /// **Note:** Axis 2 swaps xi/eta in linear terms (PV2_1 multiplies eta, PV2_2 multiplies xi).
  /// Missing coefficients default to 0, except PV2_1 defaults to 1.0.
  #[must_use]
  pub fn new_axis2(coeffs: &[f64]) -> Self {
    Self::new_axis(coeffs, TpvAxis::Axis2)
  }

  /// Internal helper to create TPV coefficients for a given axis.
  fn new_axis(coeffs: &[f64], axis: TpvAxis) -> Self {
    let mut coefficients = [0.0; 40];
    // Default PV1_1 or PV2_1 = 1
    coefficients[1] = 1.0;

    for (i, &val) in coeffs.iter().enumerate().take(40) {
      coefficients[i] = val;
    }

    Self { coefficients, axis }
  }

  /// Evaluate the polynomial at (xi, eta).
  ///
  /// Returns xi' for axis 1, eta' for axis 2.
  ///
  /// # Units
  /// * Input: `xi`, `eta` in DEGREES (intermediate world coordinates after CD matrix)
  /// * Output: corrected intermediate world coordinate (also in DEGREES)
  ///
  /// TPV is a SEQUENT distortion per WCSLIB - it operates on intermediate world
  /// coordinates (DEGREES) that have already been transformed by the CD matrix.
  /// Unlike SIP (which is a PRIOR distortion operating on pixel offsets), TPV
  /// computes corrected coordinates directly, not an additive correction.
  /// Per TPV specification: "the units of xi, eta, f_xi, and f_eta are also degrees."
  #[must_use]
  pub fn eval(&self, xi: f64, eta: f64) -> f64 {
    let xi2 = xi * xi;
    let eta2 = eta * eta;
    let xi3 = xi2 * xi;
    let eta3 = eta2 * eta;
    let xi4 = xi2 * xi2;
    let eta4 = eta2 * eta2;
    let xi5 = xi4 * xi;
    let eta5 = eta4 * eta;
    let xi6 = xi3 * xi3;
    let eta6 = eta3 * eta3;
    let xi7 = xi6 * xi;
    let eta7 = eta6 * eta;

    let r = (xi2 + eta2).sqrt();
    let r2 = r * r;
    let r3 = r2 * r;
    let r5 = r3 * r2;
    let r7 = r5 * r2;

    let c = &self.coefficients;

    match self.axis {
      TpvAxis::Axis1 => {
        c[0]
          + c[1] * xi
          + c[2] * eta
          + c[3] * r
          + c[4] * xi2
          + c[5] * xi * eta
          + c[6] * eta2
          + c[7] * xi3
          + c[8] * xi2 * eta
          + c[9] * xi * eta2
          + c[10] * eta3
          + c[11] * r3
          + c[12] * xi4
          + c[13] * xi3 * eta
          + c[14] * xi2 * eta2
          + c[15] * xi * eta3
          + c[16] * eta4
          + c[17] * xi5
          + c[18] * xi4 * eta
          + c[19] * xi3 * eta2
          + c[20] * xi2 * eta3
          + c[21] * xi * eta4
          + c[22] * eta5
          + c[23] * r5
          + c[24] * xi6
          + c[25] * xi5 * eta
          + c[26] * xi4 * eta2
          + c[27] * xi3 * eta3
          + c[28] * xi2 * eta4
          + c[29] * xi * eta5
          + c[30] * eta6
          + c[31] * xi7
          + c[32] * xi6 * eta
          + c[33] * xi5 * eta2
          + c[34] * xi4 * eta3
          + c[35] * xi3 * eta4
          + c[36] * xi2 * eta5
          + c[37] * xi * eta6
          + c[38] * eta7
          + c[39] * r7
      }
      TpvAxis::Axis2 => {
        c[0]
          + c[1] * eta
          + c[2] * xi
          + c[3] * r
          + c[4] * eta2
          + c[5] * eta * xi
          + c[6] * xi2
          + c[7] * eta3
          + c[8] * eta2 * xi
          + c[9] * eta * xi2
          + c[10] * xi3
          + c[11] * r3
          + c[12] * eta4
          + c[13] * eta3 * xi
          + c[14] * eta2 * xi2
          + c[15] * eta * xi3
          + c[16] * xi4
          + c[17] * eta5
          + c[18] * eta4 * xi
          + c[19] * eta3 * xi2
          + c[20] * eta2 * xi3
          + c[21] * eta * xi4
          + c[22] * xi5
          + c[23] * r5
          + c[24] * eta6
          + c[25] * eta5 * xi
          + c[26] * eta4 * xi2
          + c[27] * eta3 * xi3
          + c[28] * eta2 * xi4
          + c[29] * eta * xi5
          + c[30] * xi6
          + c[31] * eta7
          + c[32] * eta6 * xi
          + c[33] * eta5 * xi2
          + c[34] * eta4 * xi3
          + c[35] * eta3 * xi4
          + c[36] * eta2 * xi5
          + c[37] * eta * xi6
          + c[38] * xi7
          + c[39] * r7
      }
    }
  }

  /// Partial derivative df/dxi (used for Newton-Raphson inverse).
  #[must_use]
  pub fn dxi(&self, xi: f64, eta: f64) -> f64 {
    let xi2 = xi * xi;
    let eta2 = eta * eta;
    let xi3 = xi2 * xi;
    let eta3 = eta2 * eta;
    let xi4 = xi2 * xi2;
    let eta4 = eta2 * eta2;
    let xi5 = xi4 * xi;
    let eta5 = eta4 * eta;
    let xi6 = xi3 * xi3;
    let eta6 = eta3 * eta3;

    let r2 = xi2 + eta2;
    let r = r2.sqrt();

    // dr/dxi = xi/r (when r > 0)
    let drxi = if r > 1e-15 { xi / r } else { 0.0 };
    let c = &self.coefficients;

    match self.axis {
      TpvAxis::Axis1 => {
        c[1]
          + c[3] * drxi
          + 2.0 * c[4] * xi
          + c[5] * eta
          + 3.0 * c[7] * xi2
          + 2.0 * c[8] * xi * eta
          + c[9] * eta2
          + 3.0 * c[11] * r2 * drxi
          + 4.0 * c[12] * xi3
          + 3.0 * c[13] * xi2 * eta
          + 2.0 * c[14] * xi * eta2
          + c[15] * eta3
          + 5.0 * c[17] * xi4
          + 4.0 * c[18] * xi3 * eta
          + 3.0 * c[19] * xi2 * eta2
          + 2.0 * c[20] * xi * eta3
          + c[21] * eta4
          + 5.0 * c[23] * r2 * r2 * drxi
          + 6.0 * c[24] * xi5
          + 5.0 * c[25] * xi4 * eta
          + 4.0 * c[26] * xi3 * eta2
          + 3.0 * c[27] * xi2 * eta3
          + 2.0 * c[28] * xi * eta4
          + c[29] * eta5
          + 7.0 * c[31] * xi6
          + 6.0 * c[32] * xi5 * eta
          + 5.0 * c[33] * xi4 * eta2
          + 4.0 * c[34] * xi3 * eta3
          + 3.0 * c[35] * xi2 * eta4
          + 2.0 * c[36] * xi * eta5
          + c[37] * eta6
          + 7.0 * c[39] * r2 * r2 * r2 * drxi
      }
      TpvAxis::Axis2 => {
        c[2]
          + c[3] * drxi
          + c[5] * eta
          + 2.0 * c[6] * xi
          + c[8] * eta2
          + 2.0 * c[9] * eta * xi
          + 3.0 * c[10] * xi2
          + 3.0 * c[11] * r2 * drxi
          + c[13] * eta3
          + 2.0 * c[14] * eta2 * xi
          + 3.0 * c[15] * eta * xi2
          + 4.0 * c[16] * xi3
          + c[18] * eta4
          + 2.0 * c[19] * eta3 * xi
          + 3.0 * c[20] * eta2 * xi2
          + 4.0 * c[21] * eta * xi3
          + 5.0 * c[22] * xi4
          + 5.0 * c[23] * r2 * r2 * drxi
          + c[25] * eta5
          + 2.0 * c[26] * eta4 * xi
          + 3.0 * c[27] * eta3 * xi2
          + 4.0 * c[28] * eta2 * xi3
          + 5.0 * c[29] * eta * xi4
          + 6.0 * c[30] * xi5
          + c[32] * eta6
          + 2.0 * c[33] * eta5 * xi
          + 3.0 * c[34] * eta4 * xi2
          + 4.0 * c[35] * eta3 * xi3
          + 5.0 * c[36] * eta2 * xi4
          + 6.0 * c[37] * eta * xi5
          + 7.0 * c[38] * xi6
          + 7.0 * c[39] * r2 * r2 * r2 * drxi
      }
    }
  }

  /// Partial derivative df/deta (used for Newton-Raphson inverse).
  #[must_use]
  pub fn deta(&self, xi: f64, eta: f64) -> f64 {
    let xi2 = xi * xi;
    let eta2 = eta * eta;
    let xi3 = xi2 * xi;
    let eta3 = eta2 * eta;
    let xi4 = xi2 * xi2;
    let eta4 = eta2 * eta2;
    let xi5 = xi4 * xi;
    let eta5 = eta4 * eta;
    let xi6 = xi3 * xi3;
    let eta6 = eta3 * eta3;

    let r2 = xi2 + eta2;
    let r = r2.sqrt();

    // dr/deta = eta/r (when r > 0)
    let dreta = if r > 1e-15 { eta / r } else { 0.0 };
    let c = &self.coefficients;

    match self.axis {
      TpvAxis::Axis1 => {
        c[2]
          + c[3] * dreta
          + c[5] * xi
          + 2.0 * c[6] * eta
          + c[8] * xi2
          + 2.0 * c[9] * xi * eta
          + 3.0 * c[10] * eta2
          + 3.0 * c[11] * r2 * dreta
          + c[13] * xi3
          + 2.0 * c[14] * xi2 * eta
          + 3.0 * c[15] * xi * eta2
          + 4.0 * c[16] * eta3
          + c[18] * xi4
          + 2.0 * c[19] * xi3 * eta
          + 3.0 * c[20] * xi2 * eta2
          + 4.0 * c[21] * xi * eta3
          + 5.0 * c[22] * eta4
          + 5.0 * c[23] * r2 * r2 * dreta
          + c[25] * xi5
          + 2.0 * c[26] * xi4 * eta
          + 3.0 * c[27] * xi3 * eta2
          + 4.0 * c[28] * xi2 * eta3
          + 5.0 * c[29] * xi * eta4
          + 6.0 * c[30] * eta5
          + c[32] * xi6
          + 2.0 * c[33] * xi5 * eta
          + 3.0 * c[34] * xi4 * eta2
          + 4.0 * c[35] * xi3 * eta3
          + 5.0 * c[36] * xi2 * eta4
          + 6.0 * c[37] * xi * eta5
          + 7.0 * c[38] * eta6
          + 7.0 * c[39] * r2 * r2 * r2 * dreta
      }
      TpvAxis::Axis2 => {
        c[1]
          + c[3] * dreta
          + 2.0 * c[4] * eta
          + c[5] * xi
          + 3.0 * c[7] * eta2
          + 2.0 * c[8] * eta * xi
          + c[9] * xi2
          + 3.0 * c[11] * r2 * dreta
          + 4.0 * c[12] * eta3
          + 3.0 * c[13] * eta2 * xi
          + 2.0 * c[14] * eta * xi2
          + c[15] * xi3
          + 5.0 * c[17] * eta4
          + 4.0 * c[18] * eta3 * xi
          + 3.0 * c[19] * eta2 * xi2
          + 2.0 * c[20] * eta * xi3
          + c[21] * xi4
          + 5.0 * c[23] * r2 * r2 * dreta
          + 6.0 * c[24] * eta5
          + 5.0 * c[25] * eta4 * xi
          + 4.0 * c[26] * eta3 * xi2
          + 3.0 * c[27] * eta2 * xi3
          + 2.0 * c[28] * eta * xi4
          + c[29] * xi5
          + 7.0 * c[31] * eta6
          + 6.0 * c[32] * eta5 * xi
          + 5.0 * c[33] * eta4 * xi2
          + 4.0 * c[34] * eta3 * xi3
          + 3.0 * c[35] * eta2 * xi4
          + 2.0 * c[36] * eta * xi5
          + c[37] * xi6
          + 7.0 * c[39] * r2 * r2 * r2 * dreta
      }
    }
  }
}

/// TPV coefficients for both axes
///
/// Combines the polynomial coefficients for both the xi and eta directions
/// into a single structure for applying the complete TPV transformation.
///
/// # Example
///
/// ```rust
/// use mapproj::tpv::{TpvCoeff, TpvPV};
///
/// // Create axis-specific coefficients
/// let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.01]); // xi' = xi + 0.01*r
/// let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.01]); // eta' = eta + 0.01*r
///
/// // Combine into single TPV parameter set
/// let tpv_params = TpvPV::new(pv1, pv2);
/// ```
#[derive(Debug, Clone)]
pub struct TpvPV {
  /// First axis coefficients (PV1_n)
  pv1: TpvCoeff,

  /// Second axis coefficients (PV2_n)
  pv2: TpvCoeff,
}

impl TpvPV {
  /// Create TPV coefficients for both axes.
  #[must_use]
  pub fn new(pv1: TpvCoeff, pv2: TpvCoeff) -> Self {
    Self { pv1, pv2 }
  }
}

/// TPV distortion transformer.
///
/// Applies polynomial distortion with forward `project()` and inverse `inverse()` methods.
/// Inverse uses Newton-Raphson iteration.
#[derive(Debug, Clone)]
pub struct Tpv {
  /// TPV coefficients for both axes
  pv: TpvPV,

  /// Maximum iterations for Newton-Raphson
  n_iter: u8,

  /// Precision for Newton-Raphson convergence
  eps: f64,
}

impl Tpv {
  /// Create a new TPV distortion handler.
  #[must_use]
  pub fn new(pv: TpvPV) -> Self {
    Self {
      pv,
      n_iter: 100,
      eps: 1.0e-9,
    }
  }

  /// Apply forward TPV distortion: (xi, eta) -> (xi', eta').
  #[must_use]
  pub fn project(&self, xi: f64, eta: f64) -> (f64, f64) {
    let xi_prime = self.pv.pv1.eval(xi, eta);
    let eta_prime = self.pv.pv2.eval(xi, eta);

    (xi_prime, eta_prime)
  }

  /// Apply inverse TPV distortion: (xi', eta') -> (xi, eta) via Newton-Raphson.
  #[must_use]
  pub fn inverse(&self, xi_prime: f64, eta_prime: f64) -> Option<ImgXY> {
    self.bivariate_newton(xi_prime, eta_prime)
  }

  /// Newton-Raphson solver for inverse transformation.
  #[must_use]
  fn bivariate_newton(&self, xi_prime: f64, eta_prime: f64) -> Option<ImgXY> {
    // Try multiple initial guesses
    let attempts = [
      // Often a good guess for small distortions
      (xi_prime, eta_prime),
      (0.0_f64, 0.0_f64),
      (1e-3, 0.0),
      (-1e-3, 0.0),
      (0.0, 1e-3),
      (0.0, -1e-3),
      (1e-2, 1e-2),
      (-1e-2, -1e-2),
    ];

    let eps2 = self.eps.pow2();

    for (dxi, deta) in attempts.iter().copied() {
      let mut xi = xi_prime + dxi;
      let mut eta = eta_prime + deta;

      // Initial residuals
      let mut f = self.pv.pv1.eval(xi, eta) - xi_prime;
      let mut g = self.pv.pv2.eval(xi, eta) - eta_prime;

      let mut norm2 = f.pow2() + g.pow2();
      let mut i = 0;

      while i < self.n_iter && norm2 > eps2 {
        // Jacobian: [df_xi/dxi, df_xi/deta]
        //           [df_eta/dxi, df_eta/deta]
        let a = self.pv.pv1.dxi(xi, eta);
        let b = self.pv.pv1.deta(xi, eta);
        let c = self.pv.pv2.dxi(xi, eta);
        let d = self.pv.pv2.deta(xi, eta);

        let jac = a * d - b * c;

        // Check for singular Jacobian
        let jac_tol = f64::EPSILON * (a.abs() + d.abs()).max(1.0);
        if jac.is_nan() || jac.is_infinite() || jac.abs() <= jac_tol {
          break;
        }

        let inv_jac = 1.0 / jac;

        // Newton step
        let step_xi = inv_jac * (f * d - g * b);
        let step_eta = inv_jac * (g * a - f * c);

        // Damping for stability
        let max_step = 1e6_f64;
        let step_norm = step_xi.abs().max(step_eta.abs());
        let (dxi_step, deta_step) = if step_norm > max_step {
          (
            step_xi * (max_step / step_norm),
            step_eta * (max_step / step_norm),
          )
        } else {
          (step_xi, step_eta)
        };

        xi -= dxi_step;
        eta -= deta_step;

        f = self.pv.pv1.eval(xi, eta) - xi_prime;
        g = self.pv.pv2.eval(xi, eta) - eta_prime;
        norm2 = f.pow2() + g.pow2();
        i += 1;
      }

      // Check convergence
      if norm2 <= eps2 {
        return Some(ImgXY::new(xi, eta));
      }
    }

    None
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_tpv_identity() {
    // With default coefficients (PV1_1=1, PV2_1=1, all others 0),
    // TPV should be identity transformation
    let pv1 = TpvCoeff::new_axis1(&[]);
    let pv2 = TpvCoeff::new_axis2(&[]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 3.5;
    let eta = -2.7;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    // For identity: xi' = xi, eta' = eta
    assert!((xi_p - xi).abs() < 1e-12);
    assert!((eta_p - eta).abs() < 1e-12);
  }

  #[test]
  fn test_tpv_linear() {
    // Simple linear transformation: xi' = 2*xi, eta' = 3*eta
    // PV1_0=0, PV1_1=2
    let pv1 = TpvCoeff::new_axis1(&[0.0, 2.0]);
    // PV2_0=0, PV2_1=3
    let pv2 = TpvCoeff::new_axis2(&[0.0, 3.0]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 2.0;
    let eta = 1.5;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    assert!((xi_p - 2.0 * xi).abs() < 1e-12);
    assert!((eta_p - 3.0 * eta).abs() < 1e-12);
  }

  #[test]
  fn test_tpv_radial() {
    // Test radial term: xi' = xi + 0.1*r, eta' = eta + 0.1*r
    // PV1_3 = 0.1
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.1]);
    // PV2_3 = 0.1
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.1]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 3.0_f64;
    let eta = 4.0_f64;
    let r = (xi * xi + eta * eta).sqrt();
    let (xi_p, eta_p) = tpv.project(xi, eta);

    assert!((xi_p - (xi + 0.1 * r)).abs() < 1e-12);
    assert!((eta_p - (eta + 0.1 * r)).abs() < 1e-12);
  }

  #[test]
  fn test_tpv_derivatives() {
    // Test derivative calculation with a quadratic polynomial
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.5, 0.0, 0.1, 0.2, 0.3]);
    let _pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.3, 0.0, 0.2, 0.1, 0.4]);

    let xi = 2.0;
    let eta = 1.5;

    // Numerical derivatives
    let eps = 1e-8;
    let f0 = pv1.eval(xi, eta);
    let f_xi = pv1.eval(xi + eps, eta);
    let f_eta = pv1.eval(xi, eta + eps);
    let dxi_numerical = (f_xi - f0) / eps;
    let deta_numerical = (f_eta - f0) / eps;

    // Analytical derivatives
    let dxi_analytical = pv1.dxi(xi, eta);
    let deta_analytical = pv1.deta(xi, eta);

    assert!((dxi_analytical - dxi_numerical).abs() < 1e-6);
    assert!((deta_analytical - deta_numerical).abs() < 1e-6);
  }

  #[test]
  fn test_tpv_inverse_linear() {
    // Test inverse for linear transformation
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.5, 0.0]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 2.0, 0.0]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 3.0;
    let eta = -2.0;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    // Invert
    let result = tpv.inverse(xi_p, eta_p).unwrap();

    assert!((result.x() - xi).abs() < 1e-8);
    assert!((result.y() - eta).abs() < 1e-8);
  }

  #[test]
  fn test_tpv_inverse_quadratic() {
    // Test inverse with quadratic distortion
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.0, 1e-3, 0.0, 1e-3]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.0, 1e-3, 0.0, 1e-3]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 10.0;
    let eta = 15.0;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    let result = tpv.inverse(xi_p, eta_p).unwrap();

    assert!((result.x() - xi).abs() < 1e-7);
    assert!((result.y() - eta).abs() < 1e-7);
  }

  #[test]
  fn test_tpv_no_bounds() {
    // Test that TPV works without bounds checking
    let pv1 = TpvCoeff::new_axis1(&[]);
    let pv2 = TpvCoeff::new_axis2(&[]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    // Test that coordinates outside the provided range still work (identity transform)
    let (xi_p, eta_p) = tpv.project(15.0, 0.0);
    assert!((xi_p - 15.0).abs() < 1e-12);
    assert!((eta_p - 0.0).abs() < 1e-12);

    // Test that inverse also works without bounds
    let result = tpv.inverse(100.0, 100.0);
    assert!(result.is_some());
    let xy = result.unwrap();
    assert!((xy.x() - 100.0).abs() < 1e-6);
    assert!((xy.y() - 100.0).abs() < 1e-6);
  }

  #[test]
  fn test_tpv_cross_terms_axis1() {
    // Verify that PV1_2 multiplies eta (not xi) in axis 1 formula
    // According to spec: xi' = PV1_0 + PV1_1*xi + PV1_2*eta + ...
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.5]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 2.0;
    let eta = 3.0;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    // xi' should be: 0.0 + 1.0*xi + 0.5*eta = 2.0 + 1.5 = 3.5
    // eta' should be: 0.0 + 1.0*eta + 0.0*xi = 3.0
    assert!((xi_p - 3.5).abs() < 1e-14, "xi' = {xi_p}, expected 3.5");
    assert!((eta_p - 3.0).abs() < 1e-14, "eta' = {eta_p}, expected 3.0");
  }

  #[test]
  fn test_tpv_cross_terms_axis2() {
    // Verify that PV2_2 multiplies xi (not eta) in axis 2 formula
    // According to spec: eta' = PV2_0 + PV2_1*eta + PV2_2*xi + ...
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.7]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 2.0;
    let eta = 3.0;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    // xi' should be: 0.0 + 1.0*xi + 0.0*eta = 2.0
    // eta' should be: 0.0 + 1.0*eta + 0.7*xi = 3.0 + 1.4 = 4.4
    assert!((xi_p - 2.0).abs() < 1e-14, "xi' = {xi_p}, expected 2.0");
    assert!((eta_p - 4.4).abs() < 1e-14, "eta' = {eta_p}, expected 4.4");
  }

  #[test]
  fn test_tpv_quadratic_comprehensive() {
    // Test from Python: PV1=[0,1,0,0,0.001,0,0.001], PV2=[0,1,0,0,0.001,0,0.001]
    // For xi=10, eta=15:
    // xi' = xi + 0.001*xi^2 + 0.001*eta^2 = 10 + 0.1 + 0.225 = 10.325
    // eta' = eta + 0.001*eta^2 + 0.001*xi^2 = 15 + 0.225 + 0.1 = 15.325
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.0, 1e-3, 0.0, 1e-3]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.0, 1e-3, 0.0, 1e-3]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 10.0;
    let eta = 15.0;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    let expected_xi_p = 10.325;
    let expected_eta_p = 15.325;

    assert!(
      (xi_p - expected_xi_p).abs() < 1e-14,
      "xi' = {xi_p}, expected {expected_xi_p}"
    );
    assert!(
      (eta_p - expected_eta_p).abs() < 1e-14,
      "eta' = {eta_p}, expected {expected_eta_p}"
    );
  }

  #[test]
  fn test_tpv_radial_higher_order() {
    // Test radial distortion: r + r^3 + r^5
    // PV1_3 = 0.01, PV1_11 = 0.0001, PV1_23 = 0.000001
    let mut pv1_coeffs = vec![0.0; 40];
    pv1_coeffs[1] = 1.0;
    pv1_coeffs[3] = 0.01;
    pv1_coeffs[11] = 0.0001;
    pv1_coeffs[23] = 0.000001;

    let mut pv2_coeffs = vec![0.0; 40];
    pv2_coeffs[1] = 1.0;
    pv2_coeffs[3] = 0.01;
    pv2_coeffs[11] = 0.0001;
    pv2_coeffs[23] = 0.000001;

    let pv1 = TpvCoeff::new_axis1(&pv1_coeffs);
    let pv2 = TpvCoeff::new_axis2(&pv2_coeffs);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    // Test point: xi=3, eta=4, r=5
    let xi = 3.0_f64;
    let eta = 4.0_f64;
    let r = (xi * xi + eta * eta).sqrt();
    assert!((r - 5.0).abs() < 1e-10);

    let (xi_p, eta_p) = tpv.project(xi, eta);

    // xi' = xi + 0.01*r + 0.0001*r^3 + 0.000001*r^5
    // r=5: xi' = 3 + 0.05 + 0.0125 + 0.003125 = 3.065625
    // eta' = 4 + 0.05 + 0.0125 + 0.003125 = 4.065625
    let expected_xi_p = 3.065625;
    let expected_eta_p = 4.065625;

    assert!(
      (xi_p - expected_xi_p).abs() < 1e-14,
      "xi' = {xi_p}, expected {expected_xi_p}"
    );
    assert!(
      (eta_p - expected_eta_p).abs() < 1e-14,
      "eta' = {eta_p}, expected {expected_eta_p}"
    );
  }

  #[test]
  fn test_tpv_asymmetric() {
    // Test asymmetric distortion from Python verification
    // PV1: [0, 1, 0.05, 0, 0.001, 0, 0.0005]
    // PV2: [0, 1, -0.03, 0, 0.0008, 0, 0.0012]
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.05, 0.0, 0.001, 0.0, 0.0005]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, -0.03, 0.0, 0.0008, 0.0, 0.0012]);
    let pv = TpvPV::new(pv1, pv2);
    let tpv = Tpv::new(pv);

    let xi = 2.0;
    let eta = 3.0;
    let (xi_p, eta_p) = tpv.project(xi, eta);

    // xi' = 0 + 1*xi + 0.05*eta + 0.001*xi^2 + 0.0005*eta^2
    //     = 2 + 0.15 + 0.004 + 0.0045 = 2.1585
    // eta' = 0 + 1*eta + (-0.03)*xi + 0.0008*eta^2 + 0.0012*xi^2
    //      = 3 - 0.06 + 0.0072 + 0.0048 = 2.952
    let expected_xi_p = 2.1585;
    let expected_eta_p = 2.952;

    assert!(
      (xi_p - expected_xi_p).abs() < 1e-14,
      "xi' = {xi_p}, expected {expected_xi_p}"
    );
    assert!(
      (eta_p - expected_eta_p).abs() < 1e-14,
      "eta' = {eta_p}, expected {expected_eta_p}"
    );
  }
}
