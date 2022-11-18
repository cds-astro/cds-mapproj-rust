//! Module containing the structure to convert back on forth 
//! from Image coordinates to Celestial coordinates.

use crate::img2proj::{ImgXY2ProjXY, ProjXY2ImgXY};
use crate::{CanonicalProjection, CenteredProjection, ImgXY, LonLat, Projection, XYZ};

/// Structure to convert back on forth from Image coordinates to Celestial coordinates.
pub struct Img2Celestial<P: CanonicalProjection> {
  img2proj: ImgXY2ProjXY,
  proj2img: ProjXY2ImgXY,
  proj: CenteredProjection<P>,
}

impl<P: CanonicalProjection> Img2Celestial<P> {
  
  pub fn new(img2proj: ImgXY2ProjXY, proj: CenteredProjection<P>) -> Self {
    let proj2img = img2proj.inverse();
    Self { img2proj, proj2img, proj }
  }

  /// Change the projection center.
  /// # Param
  /// * `lonlat`: new projection center
  pub fn set_proj_center_from_lonlat(&mut self, lonlat: &LonLat) {
    self.proj.set_proj_center_from_lonlat(lonlat)
  }

  /// Change the projection center.
  /// # Param
  /// * `xyz`: new projection center
  pub fn set_proj_center_from_xyz(&mut self, xyz: &XYZ) {
    self.proj.set_proj_center_from_xyz(xyz)
  }
  
  pub fn lonlat2img(&self, lonlat: &LonLat) -> Option<ImgXY> {
    self.proj.proj_lonlat(lonlat).and_then(|xy| self.proj2img.proj2img(&xy))
  }

  pub fn xyz2img(&self, xyz: &XYZ) -> Option<ImgXY> {
    self.proj.proj_xyz(xyz).and_then(|xy| self.proj2img.proj2img(&xy))
  }
  
  pub fn img2lonlat(&self, img_pos: &ImgXY) -> Option<LonLat> {
    self.proj.unproj_lonlat(&self.img2proj.img2proj(img_pos))
  }

  pub fn img2xyz(&self, img_pos: &ImgXY) -> Option<XYZ> {
    self.proj.unproj_xyz(&self.img2proj.img2proj(img_pos))
  }
  
}
