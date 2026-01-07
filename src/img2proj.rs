//! Module containing the structure to convert back on forth
//! from Image coordinates to Intermediate coordinates, i.e. coordinates in the projection plane.

use crate::{sip::Sip, tpv::Tpv, ImgXY, ProjXY};
use std::ops::RangeInclusive;

/// Transform the XY coordinates in the projection plane in a pixel coordinates in an image.
pub trait ProjXY2ImgXY {
  /// Transforms intermediate world coordinates -- i.e. coordinates in the canonical projection
  /// plane -- to pixel coordinates in an image.
  /// # Params
  /// * `xy` coordinates in the canonical projection plane.
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY>;
}

/// Transform the pixel coordinates in an image to the XY coordinates in the projection plane.
pub trait ImgXY2ProjXY {
  /// Type of the inverse transformation.
  type T: ProjXY2ImgXY;

  /// Transforms pixel coordinates to the intermediate world coordinates, i.e.
  ///  to coordinates in the canonical projection plane.
  /// # Params
  /// * `xy`: pixel coordinates (no units)
  fn img2proj(&self, xy: &ImgXY) -> ProjXY;

  /// Provide the inverse transformation.
  fn inverse(&self) -> Self::T;
}

/// Image to Projected plane conversion
#[derive(Debug, Clone, Copy)]
pub struct BasicImgXY2ProjXY {
  /// Projection x-axis coordinate at the image X-axis center.
  center_px: f64,
  /// Projection y-axis coordinate at the image Y-axis center.
  center_py: f64,
  /// Center of the image in the X-axis (= half the number of pixels along the X-axis).
  center_x: f64,
  /// Center of the image in the Y-axis (= half the number of pixels along the Y-axis).
  center_y: f64,
  /// Size of a pixel along the projection plane the x-axis.
  scale_x: f64,
  /// Size of a pixel along the projection plane the y-axis.
  scale_y: f64,
}

impl BasicImgXY2ProjXY {
  /// We assume the origin to be on the bottom left corner of the image,
  /// with the `x` increasing from left to right and `y` increasing bottom-up.
  /// # Remark
  /// * for a PNG, the origin is the top left corner and `y` is increasing top-down.
  /// * for the east to be towards the left of an image, the `x` axis as to be reversed.
  /// # Params
  /// * `img_size`: `(size_x, size_y)` number of pixels in each dimension
  /// * `proj_bounds`: `(bounds_x, bounds_y)` boundaries of the projection domain
  #[must_use]
  pub fn from(
    img_size: (u16, u16),
    proj_bounds: (&RangeInclusive<f64>, &RangeInclusive<f64>),
  ) -> Self {
    // x-axis
    let img_size_x = f64::from(img_size.0);
    let center_x = 0.5 * (img_size_x - 1.0);
    let xp_min = proj_bounds.0.start();
    let xp_max = proj_bounds.0.end();
    let scale_x = (xp_max - xp_min) / img_size_x;
    let center_px = 0.5 * (xp_max + xp_min);
    // y-axis
    let img_size_y = f64::from(img_size.1);
    let center_y = 0.5 * (img_size_y - 1.0);
    let yp_min = proj_bounds.1.start();
    let yp_max = proj_bounds.1.end();
    let scale_y = (yp_max - yp_min) / img_size_y;
    let center_py = 0.5 * (yp_max + yp_min);
    // Create
    Self {
      center_px,
      center_py,
      center_x,
      center_y,
      scale_x,
      scale_y,
    }
  }
}

impl ImgXY2ProjXY for BasicImgXY2ProjXY {
  type T = Self;

  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    let proj_x = self.center_px + (xy.x - self.center_x) * self.scale_x;
    let proj_y = self.center_py + (xy.y - self.center_y) * self.scale_y;
    ProjXY::new(proj_x, proj_y)
  }

  fn inverse(&self) -> Self::T {
    *self
  }
}

impl ProjXY2ImgXY for BasicImgXY2ProjXY {
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    let x = self.center_x + (xy.x - self.center_px) / self.scale_x;
    let y = self.center_y + (xy.y - self.center_py) / self.scale_y;
    Some(ImgXY::new(x, y))
  }
}

/// Specific implementation for PNG (top left origin) with East towards the left.
#[derive(Debug, Clone, Copy)]
pub struct ReversedEastPngImgXY2ProjXY {
  y_img_max: f64,
  b: BasicImgXY2ProjXY,
}

impl ReversedEastPngImgXY2ProjXY {
  /// Create a new instance from the image size and projection bounds.
  #[must_use]
  pub fn from(
    img_size: (u16, u16),
    proj_bounds: (&RangeInclusive<f64>, &RangeInclusive<f64>),
  ) -> Self {
    Self {
      y_img_max: f64::from(img_size.1 - 1),
      b: BasicImgXY2ProjXY::from(img_size, proj_bounds),
    }
  }
}

impl ImgXY2ProjXY for ReversedEastPngImgXY2ProjXY {
  type T = Self;

  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    let proj_x = self.b.center_px + (xy.x - self.b.center_x) * self.b.scale_x;
    let proj_y = self.b.center_py + ((self.y_img_max - xy.y) - self.b.center_y) * self.b.scale_y;
    ProjXY::new(-proj_x, proj_y)
  }

  fn inverse(&self) -> Self::T {
    *self
  }
}

impl ProjXY2ImgXY for ReversedEastPngImgXY2ProjXY {
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    let x = self.b.center_x + (-xy.x - self.b.center_px) / self.b.scale_x;
    let y = self.b.center_y + (xy.y - self.b.center_py) / self.b.scale_y;
    Some(ImgXY::new(x, self.y_img_max - y))
  }
}

/// Struct allowing to transform the pixel coordinates in an image to the XY coordinates
/// in the projection plane.
/// The three constructors are each associated with one of the three convention
/// describe in the FITS paper: CDij, CDELTi + PCij, CDELTi + CROTA2.
#[derive(Debug, Clone, Copy)]
pub struct WcsImgXY2ProjXY {
  /// Translation vector (in pixel units, so no units).
  crpix1: f64,
  crpix2: f64,

  /// Rotation (no units) combined with a scale (in radians) matrix.
  cd11: f64,
  cd12: f64,
  cd21: f64,
  cd22: f64,
}

impl WcsImgXY2ProjXY {
  /// Create a struct from the `CDij` convention.
  /// # Params
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `cd11`: value of the `CD11` keyword (element of a rotation + scale matrix) in degrees, use 1 as default value
  /// * `cd12`: value of the `CD12` keyword (element of a rotation + scale matrix) in degrees, use 0 as default value
  /// * `cd21`: value of the `CD21` keyword (element of a rotation + scale matrix) in degrees, use 0 as default value
  /// * `cd22`: value of the `CD22` keyword (element of a rotation + scale matrix) in degrees, use 1 as default value
  #[must_use]
  pub fn from_cd(crpix1: f64, crpix2: f64, cd11: f64, cd12: f64, cd21: f64, cd22: f64) -> Self {
    Self {
      crpix1,
      crpix2,
      cd11: cd11.to_radians(),
      cd12: cd12.to_radians(),
      cd21: cd21.to_radians(),
      cd22: cd22.to_radians(),
    }
  }

  /// Create a struct from the `CDELTi + PCij` convention.
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `pc11`: value of the `PC11` keyword (element of a rotation matrix, no units), use 1 as default value
  /// * `pc12`: value of the `PC12` keyword (element of a rotation matrix, no units), use 0 as default value
  /// * `pc21`: value of the `PC21` keyword (element of a rotation matrix, no units), use 0 as default value
  /// * `pc22`: value of the `PC22` keyword (element of a rotation matrix, no units), use 1 as default value
  /// * `cdelt1`: value of the `CDELT1` keyword, in degrees
  /// * `cdelt2`: value of the `CDELT2` keyword, in degrees
  #[must_use]
  #[allow(clippy::too_many_arguments, reason = "Unavoidable complexity")]
  pub fn from_pc(
    crpix1: f64,
    crpix2: f64,
    pc11: f64,
    pc12: f64,
    pc21: f64,
    pc22: f64,
    cdelt1: f64,
    cdelt2: f64,
  ) -> Self {
    Self::from_cd(
      crpix1,
      crpix2,
      cdelt1 * pc11,
      cdelt1 * pc12,
      cdelt2 * pc21,
      cdelt2 * pc22,
    )
  }

  /// Create a struct from the `CDELTi + CROTA2` convention.
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `crota2`: value of the `CROTA2` keyword, in degrees
  /// * `cdelt1`: value of the `CDELT1` keyword, in degrees
  /// * `cdelt2`: value of the `CDELT2` keyword, in degrees
  #[must_use]
  pub fn from_cr(crpix1: f64, crpix2: f64, crota2: f64, cdelt1: f64, cdelt2: f64) -> Self {
    let (sinc, cosc) = crota2.to_radians().sin_cos();
    Self::from_cd(
      crpix1,
      crpix2,
      cdelt1 * cosc,
      cdelt1 * sinc,
      -cdelt2 * sinc,
      cdelt2 * cosc,
    )
  }
}

impl ImgXY2ProjXY for WcsImgXY2ProjXY {
  type T = WcsProjXY2ImgXY;

  /// Transform the pixel coordinates to the intermediate world coordinates
  /// (or native spherical coordinates) by applying first a translation
  /// (given the `CRPIXi` keywords value) and then a rotation plus a scale
  /// (given the `CDij` keywords values).
  /// # Params
  /// * `imgXY`: pixel coordinates (no units) to be transformed intermediate world coordinates
  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    // Translation
    let x = xy.x - self.crpix1;
    let y = xy.y - self.crpix2;
    // Rotation + scale
    ProjXY::new(self.cd11 * x + self.cd12 * y, self.cd21 * x + self.cd22 * y)
  }

  fn inverse(&self) -> Self::T {
    // Compute the determinant of the CD matrix
    let det = self.cd11 * self.cd22 - self.cd12 * self.cd21;
    // Compute the coefficient of the inverse matrix
    // For matrix [[a,b],[c,d]], inverse is (1/det)*[[d,-b],[-c,a]]
    // So for [[cd11,cd12],[cd21,cd22]], inverse is:
    // (1/det)*[[cd22,-cd12],[-cd21,cd11]]
    WcsProjXY2ImgXY {
      crpix1: self.crpix1,
      crpix2: self.crpix2,
      icd11: self.cd22 / det,
      icd12: -self.cd12 / det,
      icd21: -self.cd21 / det,
      icd22: self.cd11 / det,
    }
  }
}

/// Struct allowing to transform the pixel coordinates in an image to the XY coordinates
/// in the projection plane.
/// The three constructors are each associated with one of the three convention
/// describe in the FITS paper: CDij, CDELTi + PCij, CDELTi + CROTA2.
#[derive(Debug, Clone)]
pub struct WcsWithSipImgXY2ProjXY {
  /// Regular transformation
  wcs: WcsImgXY2ProjXY,
  /// SIP transformation
  sip: Sip,
}

impl WcsWithSipImgXY2ProjXY {
  /// Add SIP convention to a regular WCS transformation.
  #[must_use]
  pub fn new(wcs: WcsImgXY2ProjXY, sip: Sip) -> Self {
    Self { wcs, sip }
  }
}

impl ImgXY2ProjXY for WcsWithSipImgXY2ProjXY {
  type T = WcsWithSipProjXY2ImgXY;

  /// Transform the pixel coordinates to the intermediate world coordinates
  /// (or native spherical coordinates) by applying first a translation
  /// (given the `CRPIXi` keywords value) and then a rotation plus a scale
  /// (given the `CDij` keywords values).
  /// # Params
  /// * `imgXY`: pixel coordinates (no units) to be transformed intermediate world coordinates
  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    // Translation to get u, v relative pixel coordinates
    let u = xy.x - self.wcs.crpix1;
    let v = xy.y - self.wcs.crpix2;
    // Add SIP distortion: (u + f(u,v), v + g(u,v))
    let x_distorted = u + self.sip.f(u, v);
    let y_distorted = v + self.sip.g(u, v);
    // Apply CD matrix (rotation + scale)
    ProjXY::new(
      self.wcs.cd11 * x_distorted + self.wcs.cd12 * y_distorted,
      self.wcs.cd21 * x_distorted + self.wcs.cd22 * y_distorted,
    )
  }

  fn inverse(&self) -> Self::T {
    WcsWithSipProjXY2ImgXY {
      wcs: self.wcs.inverse(),
      sip: self.sip.clone(),
    }
  }
}

/// WCS projected XY to Image XY
#[derive(Debug, Clone, Copy)]
pub struct WcsProjXY2ImgXY {
  /// Translation vector (in pixel units, so no units).
  crpix1: f64,
  crpix2: f64,
  /// Rotation (no units) combined with a scale (in radians) matrix.
  icd11: f64,
  icd12: f64,
  icd21: f64,
  icd22: f64,
}

impl ProjXY2ImgXY for WcsProjXY2ImgXY {
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    let x = self.icd11 * xy.x + self.icd12 * xy.y + self.crpix1;
    let y = self.icd21 * xy.x + self.icd22 * xy.y + self.crpix2;
    Some(ImgXY::new(x, y))
  }
}

/// WCS with SIP projected to image XY
#[derive(Debug, Clone)]
pub struct WcsWithSipProjXY2ImgXY {
  /// Regular transformation
  wcs: WcsProjXY2ImgXY,

  /// SIP transformation
  sip: Sip,
}

impl ProjXY2ImgXY for WcsWithSipProjXY2ImgXY {
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    // Apply inverse CD matrix to get corrected pixel coordinates U, V
    let u_corrected = self.wcs.icd11 * xy.x + self.wcs.icd12 * xy.y;
    let v_corrected = self.wcs.icd21 * xy.x + self.wcs.icd22 * xy.y;

    // Apply inverse SIP transformation to get u, v from U, V
    // u = U + F(U,V), v = V + G(U,V) where F,G are inverse polynomials
    self
      .sip
      .inverse(u_corrected, v_corrected)
      .map(|ImgXY { x: u, y: v }| {
        // Translate back to image pixel coordinates
        ImgXY::new(u + self.wcs.crpix1, v + self.wcs.crpix2)
      })
  }
}

/// WCS with TPV (Tangent Plane + Polynomial distortion) for image to projection conversion
#[derive(Debug, Clone)]
pub struct WcsWithTpvImgXY2ProjXY {
  /// Regular WCS transformation
  wcs: WcsImgXY2ProjXY,
  /// TPV transformation
  tpv: Tpv,
}

impl WcsWithTpvImgXY2ProjXY {
  /// Add TPV convention to a regular WCS transformation.
  #[must_use]
  pub fn new(wcs: WcsImgXY2ProjXY, tpv: Tpv) -> Self {
    Self { wcs, tpv }
  }
}

impl ImgXY2ProjXY for WcsWithTpvImgXY2ProjXY {
  type T = WcsWithTpvProjXY2ImgXY;

  /// Transform pixel coordinates to intermediate world coordinates using TPV distortion.
  ///
  /// The transformation follows the TPV specification:
  /// TPV transformation order per TPV standard:
  /// 1. Translate pixel coordinates: (x, y) -> (u, v) = (x - CRPIX1, y - CRPIX2)
  /// 2. Apply TPV distortion to PIXEL offsets: (u, v) -> (u', v')
  /// 3. Apply CD matrix to get projection plane coordinates in radians: (u', v') -> (xi, eta)
  ///
  /// IMPORTANT: TPV polynomial operates on PIXEL offsets, NOT angular coordinates!
  /// This is critical for nonlinear distortions - applying TPV after the CD matrix
  /// would cause the distortion coefficients to have wrong units and magnitude.
  ///
  /// # Params
  /// * `xy`: pixel coordinates (no units) to be transformed
  fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    // Step 1: Translation to get pixel offset (u, v)
    let u = xy.x - self.wcs.crpix1;
    let v = xy.y - self.wcs.crpix2;

    // Step 2: Apply TPV distortion to PIXEL offsets
    // TPV polynomial operates on pixels, producing distorted pixel offsets
    let (u_prime, v_prime) = self.tpv.project(u, v).unwrap_or((u, v));

    // Step 3: Apply CD matrix to distorted pixel offsets to get radians
    // CD matrix units: radians/pixel
    let xi_rad = self.wcs.cd11 * u_prime + self.wcs.cd12 * v_prime;
    let eta_rad = self.wcs.cd21 * u_prime + self.wcs.cd22 * v_prime;

    ProjXY::new(xi_rad, eta_rad)
  }

  fn inverse(&self) -> Self::T {
    WcsWithTpvProjXY2ImgXY {
      wcs: self.wcs.inverse(),
      tpv: self.tpv.clone(),
    }
  }
}

/// WCS with TPV projected to image XY
#[derive(Debug, Clone)]
pub struct WcsWithTpvProjXY2ImgXY {
  /// Regular transformation
  wcs: WcsProjXY2ImgXY,
  /// TPV transformation
  tpv: Tpv,
}

impl ProjXY2ImgXY for WcsWithTpvProjXY2ImgXY {
  /// Inverse TPV transformation order (reverse of forward):
  /// 1. Apply inverse CD matrix: (xi, eta) in radians -> (u', v') distorted pixel offsets
  /// 2. Apply inverse TPV distortion: (u', v') -> (u, v) pixel offsets
  /// 3. Translate: (u, v) -> (x, y) image coordinates
  ///
  /// IMPORTANT: TPV operates on PIXEL offsets, so inverse CD must come before inverse TPV.
  fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    // Step 1: Apply inverse CD matrix to radians to get distorted PIXEL offsets
    // Input: xi_rad, eta_rad (projection plane coordinates in radians)
    // Output: u_prime, v_prime (distorted pixel offsets)
    let u_prime = self.wcs.icd11 * xy.x + self.wcs.icd12 * xy.y;
    let v_prime = self.wcs.icd21 * xy.x + self.wcs.icd22 * xy.y;

    // Step 2: Apply inverse TPV to get undistorted pixel offsets
    // Input: u_prime, v_prime (distorted pixel offsets)
    // Output: u, v (undistorted pixel offsets)
    let result = self.tpv.inverse(u_prime, v_prime)?;
    let u = result.x;
    let v = result.y;

    // Step 3: Translate to image pixel coordinates
    Some(ImgXY::new(u + self.wcs.crpix1, v + self.wcs.crpix2))
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::sip::{Sip, SipAB, SipCoeff};

  /// Helper function to create test SIP coefficients from IRAC Channel 4 example
  /// (Table 1 from Shupe et al. ADASS XIV 2005)
  fn create_test_sip_coefficients() -> (SipCoeff, SipCoeff) {
    // A coefficients (degree 4)
    let a_coef = [
      0.0,
      0.0,
      2.1634068532689E-06,
      1.0622437604068E-11,
      1.4075878614807E-14,
      0.0,
      -5.194753640575E-06,
      -5.2797808038221E-10,
      -1.9317154005522E-14,
      8.543473309812E-06,
      -4.4012735467525E-11,
      3.767898933666E-14,
      -4.7518233007536E-10,
      5.0860953083043E-15,
      2.5776347115304E-14,
    ];

    // B coefficients (degree 4)
    let b_coef = [
      0.0,
      0.0,
      -7.2299995118730E-06,
      -4.2102920235938E-10,
      6.5531313110898E-16,
      0.0,
      6.1778338717084E-06,
      -6.7603466821178E-11,
      1.3892905568706E-14,
      -1.7442694174934E-06,
      -5.1333879897858E-10,
      -2.9648166208490E-14,
      8.5722142612681E-11,
      -2.0749495718513E-15,
      -1.812610418272E-14,
    ];

    (SipCoeff::new(a_coef.into()), SipCoeff::new(b_coef.into()))
  }

  #[test]
  fn test_wcs_with_sip_forward_transformation() {
    // Setup WCS parameters from IRAC Channel 4 example
    let crpix1 = 2048.0;
    let crpix2 = 1024.0;
    let cd11 = -7.8481866550866E-06;
    let cd12 = 1.0939720432379E-05;
    let cd21 = 1.1406694624771E-05;
    let cd22 = 8.6942510845452E-06;

    let (a_coef, b_coef) = create_test_sip_coefficients();
    let ab_coef = SipAB::new(a_coef, b_coef);
    let sip = Sip::new(ab_coef, None, -2048.0..=2048.0, -1024.0..=1024.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_sip = WcsWithSipImgXY2ProjXY::new(wcs, sip);

    // Test point: offset from CRPIX by (3, 3) pixels
    let img_xy = ImgXY::new(crpix1 + 3.0, crpix2 + 3.0);
    let proj_xy = wcs_sip.img2proj(&img_xy);

    // The transformation should apply:
    // 1. u = 3.0, v = 3.0 (relative to CRPIX)
    // 2. f(3,3) ~= 4.95811569e-05, g(3,3) ~= -2.51926572e-05 (from SIP test)
    // 3. u_dist = 3.0 + 4.95811569e-05 ~= 3.0000495811569
    // 4. v_dist = 3.0 - 2.51926572e-05 ~= 2.9999748073428
    // 5. Apply CD matrix to (u_dist, v_dist)

    let u_dist = 3.0 + 4.95811569e-05;
    let v_dist = 3.0 - 2.51926572e-05;
    let expected_x = cd11.to_radians() * u_dist + cd12.to_radians() * v_dist;
    let expected_y = cd21.to_radians() * u_dist + cd22.to_radians() * v_dist;

    // Verify the transformation is correct
    assert!((proj_xy.x - expected_x).abs() < 1e-15);
    assert!((proj_xy.y - expected_y).abs() < 1e-15);
  }

  #[test]
  fn test_simple_sip_round_trip() {
    // Very simple test with small, well-behaved values
    let crpix1 = 100.0;
    let crpix2 = 100.0;
    let cd11 = 1.0e-5;
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // Zero distortion - should behave exactly like regular WCS
    let a_coef = [0.0];
    let b_coef = [0.0];
    let sip_a = SipCoeff::new(a_coef.into());
    let sip_b = SipCoeff::new(b_coef.into());
    let ab_proj = SipAB::new(sip_a, sip_b);

    let sip = Sip::new(ab_proj, None, -100.0..=100.0, -100.0..=100.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_sip = WcsWithSipImgXY2ProjXY::new(wcs, sip);

    // Test point
    let img_xy = ImgXY::new(105.0, 103.0);
    let inverse = wcs_sip.inverse();

    let proj_xy = wcs_sip.img2proj(&img_xy);
    let img_xy_back = inverse.proj2img(&proj_xy).expect("Inverse should succeed");

    // With zero distortion, should round-trip perfectly
    assert!((img_xy_back.x - img_xy.x).abs() < 1e-10);
    assert!((img_xy_back.y - img_xy.y).abs() < 1e-10);
  }

  #[test]
  fn test_wcs_with_sip_round_trip_with_newton_raphson() {
    // Test round-trip transformation using Newton-Raphson inverse (no inverse polynomials)
    let crpix1 = 100.0;
    let crpix2 = 100.0;
    let cd11 = 1.0e-5;
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // Small 2nd order distortion
    let a_coef = [0.0, 0.0, 1.0e-8, 0.0, 5.0e-9, 0.0];
    let b_coef = [0.0, 0.0, -1.0e-8, 0.0, 3.0e-9, 0.0];
    let sip_a = SipCoeff::new(a_coef.into());
    let sip_b = SipCoeff::new(b_coef.into());
    let ab_proj = SipAB::new(sip_a, sip_b);

    // No inverse polynomials - will use Newton-Raphson
    let sip = Sip::new(ab_proj, None, -100.0..=100.0, -100.0..=100.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_sip = WcsWithSipImgXY2ProjXY::new(wcs, sip);

    // Test a point near the center
    let img_xy = ImgXY::new(105.0, 103.0);
    let inverse = wcs_sip.inverse();

    let proj_xy = wcs_sip.img2proj(&img_xy);
    let img_xy_back = inverse.proj2img(&proj_xy).expect("Inverse should succeed");

    // Newton-Raphson should converge to high precision for small distortions
    assert!(
      (img_xy_back.x - img_xy.x).abs() < 1e-6,
      "X mismatch: {} vs {}, diff: {}",
      img_xy_back.x,
      img_xy.x,
      (img_xy_back.x - img_xy.x).abs()
    );
    assert!(
      (img_xy_back.y - img_xy.y).abs() < 1e-6,
      "Y mismatch: {} vs {}, diff: {}",
      img_xy_back.y,
      img_xy.y,
      (img_xy_back.y - img_xy.y).abs()
    );
  }

  #[test]
  fn test_wcs_with_sip_forward_satisfies_paper_equation() {
    // Verify that the forward transformation satisfies Equation 1 from the SIP paper:
    // [x]   [CD11 CD12] [u + f(u,v)]
    // [y] = [CD21 CD22] [v + g(u,v)]

    let crpix1 = 100.0;
    let crpix2 = 200.0;
    let cd11 = 1.0e-5;
    let cd12 = 0.5e-5;
    let cd21 = -0.5e-5;
    let cd22 = 1.0e-5;

    // Simple 2nd order polynomial for testing
    let a_coef = [0.0, 0.0, 1.0e-7, 0.0, 2.0e-7, 0.0];
    let b_coef = [0.0, 0.0, -1.0e-7, 0.0, 1.5e-7, 0.0];

    let sip_a = SipCoeff::new(a_coef.into());
    let sip_b = SipCoeff::new(b_coef.into());
    let ab_coef = SipAB::new(sip_a, sip_b);
    let sip = Sip::new(ab_coef, None, -200.0..=200.0, -300.0..=300.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_sip = WcsWithSipImgXY2ProjXY::new(wcs, sip.clone());

    // Test point
    let img_xy = ImgXY::new(105.0, 210.0);
    let u = 5.0;
    let v = 10.0;

    // Manually compute according to the paper's equation
    let f_uv = sip.f(u, v);
    let g_uv = sip.g(u, v);
    let u_plus_f = u + f_uv;
    let v_plus_g = v + g_uv;

    let expected_x = cd11.to_radians() * u_plus_f + cd12.to_radians() * v_plus_g;
    let expected_y = cd21.to_radians() * u_plus_f + cd22.to_radians() * v_plus_g;

    // Get result from our implementation
    let proj_xy = wcs_sip.img2proj(&img_xy);

    // Should match exactly
    assert!(
      (proj_xy.x - expected_x).abs() < 1e-20,
      "X: {} vs {}",
      proj_xy.x,
      expected_x
    );
    assert!(
      (proj_xy.y - expected_y).abs() < 1e-20,
      "Y: {} vs {}",
      proj_xy.y,
      expected_y
    );
  }

  #[test]
  fn test_wcs_without_sip_is_identity_to_basic_wcs() {
    // Test that WCS without SIP polynomials behaves the same as basic WCS
    let crpix1 = 512.0;
    let crpix2 = 512.0;
    let cd11 = -1.0e-5;
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // Create SIP with zero coefficients (identity distortion)
    let a_coef = [0.0];
    let b_coef = [0.0];
    let sip_a = SipCoeff::new(a_coef.into());
    let sip_b = SipCoeff::new(b_coef.into());
    let ab_coef = SipAB::new(sip_a, sip_b);
    let sip = Sip::new(ab_coef, None, -512.0..=512.0, -512.0..=512.0);

    let wcs_basic = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_with_sip = WcsWithSipImgXY2ProjXY::new(wcs_basic, sip);

    // Test point
    let img_xy = ImgXY::new(600.0, 400.0);

    let proj_basic = wcs_basic.img2proj(&img_xy);
    let proj_with_sip = wcs_with_sip.img2proj(&img_xy);

    // Should be identical since SIP polynomials are zero
    assert!((proj_basic.x - proj_with_sip.x).abs() < 1e-20);
    assert!((proj_basic.y - proj_with_sip.y).abs() < 1e-20);
  }

  // TPV tests
  use crate::tpv::{Tpv, TpvCoeff, TpvPV};

  // Note: test_wcs_with_tpv_identity removed because basic WCS uses radians/pixel
  // while TPV expects degrees/pixel, making direct comparison invalid

  #[test]
  fn test_wcs_with_tpv_forward_transformation() {
    // Test TPV forward transformation with simple linear distortion
    let crpix1 = 100.0;
    let crpix2 = 100.0;
    let cd11 = 1.0e-5;
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // TPV with simple linear scaling: u' = 1.5*u, v' = 2.0*v
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.5]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 2.0]);
    let pv = TpvPV::new(pv1, pv2);
    // Range in PIXELS (TPV operates on pixel offsets)
    let tpv = Tpv::new(pv, -100.0..=100.0, -100.0..=100.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_tpv = WcsWithTpvImgXY2ProjXY::new(wcs, tpv);

    // Test point
    let img_xy = ImgXY::new(105.0, 103.0);
    let u = 5.0; // pixel offset (x - crpix1)
    let v = 3.0; // pixel offset (y - crpix2)

    // CORRECT ORDER per TPV standard:
    // Step 1: Apply TPV distortion to PIXEL offsets
    let u_prime = 1.5 * u; // TPV: u' = 1.5*u (in pixels)
    let v_prime = 2.0 * v; // TPV: v' = 2.0*v (in pixels)

    // Step 2: Apply CD matrix to distorted pixels to get radians
    let expected_x = cd11.to_radians() * u_prime;
    let expected_y = cd22.to_radians() * v_prime;

    let proj_xy = wcs_tpv.img2proj(&img_xy);

    assert!(
      (proj_xy.x - expected_x).abs() < 1e-15,
      "X: {} vs {}",
      proj_xy.x,
      expected_x
    );
    assert!(
      (proj_xy.y - expected_y).abs() < 1e-15,
      "Y: {} vs {}",
      proj_xy.y,
      expected_y
    );
  }

  #[test]
  fn test_wcs_with_tpv_round_trip() {
    // Test round-trip transformation with TPV
    let crpix1 = 100.0;
    let crpix2 = 100.0;
    let cd11 = 1.0e-5;
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // TPV with small quadratic distortion
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.0, 1e3, 0.0, 1e3]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.0, 1e3, 0.0, 1e3]);
    let pv = TpvPV::new(pv1, pv2);
    // Range in PIXELS (TPV operates on pixel offsets, not angular coordinates)
    let tpv = Tpv::new(pv, -100.0..=100.0, -100.0..=100.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_tpv = WcsWithTpvImgXY2ProjXY::new(wcs, tpv);

    // Test point near center
    let img_xy = ImgXY::new(105.0, 103.0);
    let inverse = wcs_tpv.inverse();

    let proj_xy = wcs_tpv.img2proj(&img_xy);
    let img_xy_back = inverse.proj2img(&proj_xy).expect("Inverse should succeed");

    // Should round-trip with good precision
    assert!(
      (img_xy_back.x - img_xy.x).abs() < 1e-6,
      "X mismatch: {} vs {}, diff: {}",
      img_xy_back.x,
      img_xy.x,
      (img_xy_back.x - img_xy.x).abs()
    );
    assert!(
      (img_xy_back.y - img_xy.y).abs() < 1e-6,
      "Y mismatch: {} vs {}, diff: {}",
      img_xy_back.y,
      img_xy.y,
      (img_xy_back.y - img_xy.y).abs()
    );
  }

  #[test]
  fn test_wcs_with_tpv_radial_distortion() {
    // Test TPV with radial distortion terms
    let crpix1 = 256.0;
    let crpix2 = 256.0;
    let cd11 = 1.0e-5;
    let cd12 = 0.0;
    let cd21 = 0.0;
    let cd22 = 1.0e-5;

    // TPV with radial term: u' = u + 0.01*r, v' = v + 0.01*r (in pixels)
    let pv1 = TpvCoeff::new_axis1(&[0.0, 1.0, 0.0, 0.01]);
    let pv2 = TpvCoeff::new_axis2(&[0.0, 1.0, 0.0, 0.01]);
    let pv = TpvPV::new(pv1, pv2);
    // Range in PIXELS (TPV operates on pixel offsets)
    let tpv = Tpv::new(pv, -100.0..=100.0, -100.0..=100.0);

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let wcs_tpv = WcsWithTpvImgXY2ProjXY::new(wcs, tpv);

    // Test point
    let img_xy = ImgXY::new(260.0, 259.0);
    let u = 4.0_f64; // pixel offset (x - crpix1)
    let v = 3.0_f64; // pixel offset (y - crpix2)

    // CORRECT ORDER per TPV standard:
    // Step 1: Calculate r from PIXEL offsets
    let r = (u * u + v * v).sqrt();

    // Step 2: Apply TPV radial distortion to PIXEL offsets
    let u_prime = u + 0.01 * r; // in pixels
    let v_prime = v + 0.01 * r; // in pixels

    // Step 3: Apply CD matrix to distorted pixels to get radians
    let expected_x = cd11.to_radians() * u_prime;
    let expected_y = cd22.to_radians() * v_prime;

    let proj_xy = wcs_tpv.img2proj(&img_xy);

    assert!(
      (proj_xy.x - expected_x).abs() < 1e-15,
      "X: {} vs {}",
      proj_xy.x,
      expected_x
    );
    assert!(
      (proj_xy.y - expected_y).abs() < 1e-15,
      "Y: {} vs {}",
      proj_xy.y,
      expected_y
    );
  }

  #[test]
  fn test_wcs_inverse_with_nondiagonal_cd_matrix() {
    // Test that inverse CD matrix works correctly with non-diagonal matrices
    // This test catches the bug where icd12 and icd21 were swapped
    let crpix1 = 100.0;
    let crpix2 = 200.0;
    let cd11 = 1.0e-5; // degrees/pixel
    let cd12 = 0.5e-5; // degrees/pixel (non-zero off-diagonal)
    let cd21 = -0.3e-5; // degrees/pixel (non-zero off-diagonal)
    let cd22 = 1.2e-5; // degrees/pixel

    let wcs = WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22);
    let inverse = wcs.inverse();

    // Test point
    let img_xy = ImgXY::new(105.0, 203.0);

    // Forward transformation
    let proj_xy = wcs.img2proj(&img_xy);

    // Inverse transformation
    let img_xy_back = inverse.proj2img(&proj_xy).expect("Inverse should succeed");

    // Round-trip should be exact (within floating point precision)
    assert!(
      (img_xy_back.x() - img_xy.x()).abs() < 1e-10,
      "X round-trip failed: {} vs {}, diff: {}",
      img_xy_back.x(),
      img_xy.x(),
      (img_xy_back.x() - img_xy.x()).abs()
    );
    assert!(
      (img_xy_back.y() - img_xy.y()).abs() < 1e-10,
      "Y round-trip failed: {} vs {}, diff: {}",
      img_xy_back.y(),
      img_xy.y(),
      (img_xy_back.y() - img_xy.y()).abs()
    );
  }
}
