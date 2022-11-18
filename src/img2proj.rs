//! Module containing the structure to convert back on forth 
//! from Image coordinates to Intermediate coordinates, i.e. coordinates in the projection plane.

use crate::{
  ImgXY, ProjXY,
  sip::Sip
};

/// Struct allowing to transform the pixel coordinates in an image to the XY coordinates 
/// in the projection plane.
/// The three constructors are each associated with one of the three convention
/// describe in the FITS paper: CDij, CDELTi + PCij, CDELTi + CROTA2.
pub struct ImgXY2ProjXY {
  /// Translation vector (in pixel units, so no units).
  crpix1: f64, crpix2: f64,
  /// Rotation (no units) combined with a scale (in radians) matrix.
  cd11: f64, cd12: f64,
  cd21: f64, cd22: f64,
  /// Possible SIP transformation
  sip: Option<Sip>,
}

impl ImgXY2ProjXY {
  
  /// Create a struct from the `CDij` convention.
  /// # Params
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `cd11`: value of the `CD11` keyword (element of a rotation + scale matrix) in degrees, use 1 as default value
  /// * `cd12`: value of the `CD12` keyword (element of a rotation + scale matrix) in degrees, use 0 as default value
  /// * `cd21`: value of the `CD21` keyword (element of a rotation + scale matrix) in degrees, use 0 as default value
  /// * `cd22`: value of the `CD22` keyword (element of a rotation + scale matrix) in degrees, use 1 as default value
  pub fn from_cd(
    crpix1: f64, crpix2: f64,
    cd11: f64, cd12: f64, 
    cd21: f64, cd22: f64,
  ) -> Self {
    Self {
      crpix1, crpix2,
      cd11: cd11.to_radians(),
      cd12: cd12.to_radians(),
      cd21: cd21.to_radians(),
      cd22: cd22.to_radians(),
      sip: None
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
  pub fn from_pc(
    crpix1: f64, crpix2: f64,
    pc11: f64,  pc12: f64,
    pc21: f64,  pc22: f64,
    cdelt1: f64, cdelt2: f64
  ) -> Self {
    Self::from_cd(
      crpix1, crpix2,
      cdelt1 * pc11, cdelt1 * pc12,
      cdelt2 * pc21, cdelt2 * pc22,
    )
  }

  /// Create a struct from the `CDELTi + CROTA2` convention.
  /// * `crpix1`: value of the `CRPIX1` keyword in pixel units (so no units), use 0 as default value
  /// * `crpix2`: value of the `CRPIX2` keyword in pixel units (so no units), use 0 as default value
  /// * `crota2`: value of the `CROTA2` keyword, in degrees
  /// * `cdelt1`: value of the `CDELT1` keyword, in degrees
  /// * `cdelt2`: value of the `CDELT2` keyword, in degrees
  pub fn from_cr(
    crpix1: f64, crpix2: f64,
    crota2: f64,
    cdelt1: f64, cdelt2: f64
  ) -> Self {
    let (sinc, cosc) = crota2.to_radians().sin_cos();
    Self::from_cd(
      crpix1, crpix2, 
      cdelt1 * cosc, cdelt1 * sinc,
      -cdelt2 * sinc, cdelt2 * cosc
    )
  }
  
  /// Returnthe previous value, if any.
  pub fn set_sip(&mut self, sip: Sip) -> Option<Sip> {
    self.sip.replace(sip)
  }

  /// Remove SIP component
  pub fn rm_sip(&mut self) {
    self.sip = None
  }
  
  /// Transform the pixel coordinates to the intermediate world coordinates
  /// (or native spherical coordinates) by applying first a translation
  /// (given the `CRPIXi` keywords value) and then a rotation plus a scale
  /// (given the `CDij` keywords values).
  /// # Params
  /// * `imgXY`: pixel coordinates (no units) to be transformed intermediate world coordinates
  pub fn img2proj(&self, xy: &ImgXY) -> ProjXY {
    // Translation
    let mut x = xy.x - self.crpix1;
    let mut y = xy.y - self.crpix2;
    // Possible SIP convention
    if let Some(sip) = &self.sip {
      let tmp = x;
      x += sip.f(x, y);
      y += sip.g(tmp, y);
    } 
    // Rotation + scale
    ProjXY::new(
      self.cd11 * x + self.cd12 * y,
      self.cd21 * x + self.cd22 * y
    )
  }
  
  pub fn inverse(&self) -> ProjXY2ImgXY {
    // Compute the determinant of the CD matrix
    let det = self.cd11 * self.cd22 - self.cd12 * self.cd21;
    // Compute the coefficient of the inverse matrix
    ProjXY2ImgXY {
      crpix1: self.crpix1, 
      crpix2: self.crpix2,
      icd11:  self.cd22 / det,
      icd12: -self.cd21 / det,
      icd21: -self.cd12 / det,
      icd22:  self.cd11 / det,
      sip: self.sip.clone()
    }
  }
  
}

pub struct ProjXY2ImgXY {
  /// Translation vector (in pixel units, so no units).
  crpix1: f64, crpix2: f64,
  /// Rotation (no units) combined with a scale (in radians) matrix.
  icd11: f64, icd12: f64,
  icd21: f64, icd22: f64,
  /// Possible SIP transformation
  sip: Option<Sip>,
}

impl ProjXY2ImgXY {
  
  pub fn proj2img(&self, xy: &ProjXY) -> Option<ImgXY> {
    let x = self.icd11 * xy.x + self.icd12 * xy.y + self.crpix1;
    let y = self.icd21 * xy.x + self.icd22 * xy.y + self.crpix2;
    if let Some(sip) = &self.sip {
      sip.inverse(x, y).map(|ImgXY{x: rx, y: ry}| ImgXY::new(x + rx, y + ry))
    } else {
      Some(ImgXY::new(x, y))
    }
   }
}



